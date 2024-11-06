"""
Script to measure whether confident interface residues between two chains of a
predicted structure are consistently predicted to be close across multiple models
(all with the same input sequences).

For the top ranked model, the script identfies all confident pairs of interface
residues that (i.e., pairs of residues in different chains that have a PAE value below
a specified cutoff) and are closer than a distance cutoff.
Then, for each subsequent model, the script calculates the fraction of these residue
pairs that are still within the distance cutoff (regardless of confidence).

Functions:
get_residue_pairs
find_confident_interface_residues
check_distances_across_models
measure_repeatibility
"""
import re
import numpy as np
from pathlib import Path
from Bio.PDB import Model
from scipy.spatial import cKDTree
from analysis_utility import extract_pae, find_rank_001_files, parse_structure_file, map_chains_and_residues
from phosphorylation_handling import correct_cif_pae

def get_residue_pairs(structure_input, distance_cutoff, abs_res_lookup_dict, all_atom):
    """
    Identify pairs of residues across different chains that are within a specified
    distance cutoff.
    
    Parameters:
        - structure_input (str, Path or Bio.PDB.Model.Model): Either a file path to the PDB file
          or a protein model object.
        - distance_cutoff (float): Maximum distance between CA atoms of residues to
          be considered as interface.
        - abs_res_lookup_dict (dict): Dictionary mapping chain and residue IDs to absolute
          residue IDs.
        - all_atom (bool): True if all atoms should be considered for measuring interface
          distance, False if only CA atoms should be considered.
    
    Returns:
        - residue_pairs (list): List of residue pairs that are within the distance cutoff,
          where each pair is a tuple of residue numbers.
    """
    # Check if structure_input is a file path or a model object
    if isinstance(structure_input, str) or isinstance(structure_input, Path):
        # Treat structure_input as a file path and parse the structure file
        structure_model = parse_structure_file(structure_input)
    elif isinstance(structure_input, Model.Model):
        # Treat structure_input as an already loaded protein model object
        structure_model = structure_input
    else:
        raise TypeError("structure_input must be either a file path (str or Path) or a protein model object (Bio.PDB.Model.Model).")

    chains = list(structure_model.get_chains())
    # Get mapping of residues to absolute residue ids
    def get_abs_res_id(abs_res_lookup_dict, target_chain_id, target_res_id):
        return abs_res_lookup_dict.get((target_chain_id, target_res_id), None)  # Returns None if the key is not found

    # Find residue pairs across chains closer than distance cutoff
    residue_pairs = []
    # Collect residue info: chain_id, res_id, abs_res_id, coord(s)
    chain_residues = {}
    for chain in chains:
        residues = []
        for res in chain.get_residues():
            res_id = get_abs_res_id(abs_res_lookup_dict, chain.id, res.id[1])
            if res_id is None:
                continue
            if not all_atom:
                if 'CA' in res:
                    coord = res['CA'].coord
                    residues.append((res_id, coord))
            else:
                atom_coords = np.array([atom.coord for atom in res.get_atoms()])
                residues.append((res_id, atom_coords))
        chain_residues[chain.id] = residues

    # Now, for each pair of chains
    residue_pairs = []
    chain_ids = list(chain_residues.keys())
    for i, chain1_id in enumerate(chain_ids[:-1]):
        residues1 = chain_residues[chain1_id]
        for chain2_id in chain_ids[i+1:]:
            residues2 = chain_residues[chain2_id]
            if not residues1 or not residues2:
                continue

            if not all_atom:
                # Prepare coordinates and KD-tree
                coords1 = np.array([r[1] for r in residues1])
                coords2 = np.array([r[1] for r in residues2])
                tree = cKDTree(coords2)
                # Query neighbors within distance_cutoff
                indices = tree.query_ball_point(coords1, r=distance_cutoff)
                for idx1, neighbors in enumerate(indices):
                    abs_res_id1 = residues1[idx1][0]
                    for idx2 in neighbors:
                        abs_res_id2 = residues2[idx2][0]
                        residue_pairs.append((abs_res_id1, abs_res_id2))
            else:
                # For all_atom=True, we need to consider all atom coordinates
                # Flatten all atom coordinates for residues in chain2
                all_coords2 = np.concatenate([r[1] for r in residues2])
                res_indices2 = np.concatenate([[r[0]] * len(r[1]) for r in residues2])
                tree = cKDTree(all_coords2)
                # Query for each atom in residues of chain1
                for res_id1, atom_coords1 in residues1:
                    if len(atom_coords1) == 0:
                        continue
                    # Query neighbors within distance_cutoff for each atom
                    neighbors = tree.query_ball_point(atom_coords1, r=distance_cutoff)
                    # Flatten the list of arrays
                    all_neighbor_indices = np.concatenate(neighbors).astype(int)
                    if len(all_neighbor_indices) > 0:
                        # Get the residue IDs corresponding to these atoms
                        abs_res_id2s = res_indices2[all_neighbor_indices]
                        # Add unique residue pairs
                        for abs_res_id2 in np.unique(abs_res_id2s):
                            residue_pairs.append((res_id1, abs_res_id2))
    return residue_pairs



def find_confident_interface_residues(structure_input, pae_input, distance_cutoff, pae_cutoff, abs_res_lookup_dict, is_pdb, all_atom):
    """
    Find pairs of interface residues that are confidently predicted to be close (in the
    top ranked model).
    
    Parameters:
        - structure_input (str, Path or Bio.PDB.Model.Model): Either a file path to the PDB file
          or a protein model object of the top ranked predicted structure.
        - pae_input (str, Path or list of lists): Either a path to the JSON file containing PAE data
          or a pre-loaded PAE matrix (e.g., list of lists or NumPy array).
        - distance_cutoff (float): Maximum distance between CA atoms of residues to
          be considered as interface.
        - pae_cutoff (float): Maximum PAE value between a pair of residues for them
          to be considered as being confidently positioned with relation to each other.
        - abs_res_lookup_dict (dict): Dictionary mapping chain and residue IDs to absolute
          residue IDs.
        - is_pdb (bool): True if the file is in PDB format. False if in CIF format.
        - all_atom (bool): True if all atoms should be considered for measuring interface
          distance, False if only CA atoms should be considered.
        
    Returns:
        - confident_pairs (list): List of residue pairs that are confidently predicted
          to be close in the top ranked model, where each pair is a tuple of residue
          numbers relative to the full sequence (i.e., not just position within the chain).
    """
    # Check if structure_input is a file path or a model object
    if isinstance(structure_input, str) or isinstance(structure_input, Path):
        # Treat structure_input as a file path and parse the structure file
        structure_model = parse_structure_file(structure_input)
    elif isinstance(structure_input, Model.Model):
        # Treat structure_input as an already loaded protein model object
        structure_model = structure_input
    else:
        raise TypeError("structure_input must be either a file path (str or Path) or a protein model object (Bio.PDB.Model.Model).")

    # Check if pae_input is a file path (string) or a pre-loaded matrix
    if isinstance(pae_input, str) or isinstance(pae_input, Path):
        # Extract PAE data from JSON file
        pae_data = extract_pae(pae_input)
    else:
        # Assume it's already a pre-loaded PAE matrix
        pae_data = pae_input

    # Find residue pairs within distance cutoff
    residue_pairs = get_residue_pairs(structure_model, distance_cutoff, abs_res_lookup_dict, all_atom)
    
    # If cif, condense pae data in case of phosphorylated residues  
    # necessary even without modfication, but unsure why
    if not is_pdb:
        pae_data = correct_cif_pae(structure_model, pae_data)

    # Check PAE values for each nearby residue pair, return pairs with PAE < cutoff
    confident_pairs = []
    for res1, res2 in residue_pairs:
        pae_value = pae_data[res1][res2]
        if pae_value < pae_cutoff:
            confident_pairs.append((res1, res2))
    return confident_pairs

def check_distances_across_models(folder_path, confident_pairs, distance_cutoff, abs_res_lookup_dict, is_pdb, all_atom):
    """
    Check if confident interface residues (from the top model) are consistently
    predicted to be close in other structure models.

    Parameters:
        - folder_path (str): Path to the folder containing the structure files
        - confident_pairs (list): List of residue pairs that are confidently predicted
          to be close in the top ranked model, where each pair is a tuple of residue
          numbers (relative to the full sequence, not just the chain position)
        - distance_cutoff (float): Maximum distance between CA atoms of residues (in
          different chains) to be considered as interface
        - abs_res_lookup_dict (dict): Dictionary mapping chain and residue IDs to absolute
          residue IDs
        - is_pdb (bool): True if the files are in PDB format, False if in CIF format
        - all_atom (bool): True if all atoms should be considered for measuring interface
          distance, False if only CA atoms should be considered
    
    Returns:
        - model_consistency_scores (list): List of tuples containing the rank number of
          each model and the proportion of confident pairs from the top model which are
          also close in the current model
        - num_consistent_models (int): Number of models with >25% consistency in interface
          residue positions
        - average_consistency_level (float): Average consistency level of models with >25%
          consistency
        - res_pairs (list): List of residue pairs that are consistently close across models,
          where each pair is a tuple of chain ID and residue number
    """
    # Find all structure files except the rank 1/top model
    # Uses rank_001 and rank_1 for pdb, model_0 for cif
    def filter_and_sort_structure_files(folder_path, file_extension, top_model_patterns, rank_pattern):
        path = Path(folder_path)
        rank_regex = re.compile(rank_pattern)
        files = [f for f in path.glob(f'*{file_extension}') if not any(pattern in f.name for pattern in top_model_patterns)]
        sorted_files = sorted(files, key=lambda x: int(rank_regex.search(x.name).group(1)))
        return sorted_files

    file_extension = '.pdb' if is_pdb else '.cif'
    top_model_patterns = ['rank_001', 'rank_1'] if is_pdb else ['model_0']
    rank_pattern = r'rank_(\d+)' if is_pdb else r'model_(\d+)'
    
    try:
        structure_files = filter_and_sort_structure_files(folder_path, file_extension, top_model_patterns, rank_pattern)
    except Exception as e:
        print(f"Error processing files: {e}")
        return []

    # Identify residue pairs within distance cutoff for each model
    model_consistency_scores = []
    consistent_res_pairs = set(confident_pairs)
    for structure_file in structure_files:
        structure_path = structure_file.as_posix()
        rank_number = int(re.search(rank_pattern, str(structure_file)).group(1))
        structure_model = parse_structure_file(structure_path, is_pdb)
        model_res_pairs = get_residue_pairs(structure_model, distance_cutoff, abs_res_lookup_dict, all_atom)
        # Calculate proportion of confident pairs from top model which are also close in current model
        score = sum(1 for (res1, res2) in confident_pairs if (res1, res2) in model_res_pairs) / len(confident_pairs) if confident_pairs else 0
        if score > 0.25:
            model_res_pairs_set = set(model_res_pairs)
            consistent_res_pairs.intersection_update(model_res_pairs_set)
        model_consistency_scores.append((rank_number, score))

    # Sort scores by model rank number
    sorted_model_consistency_scores = sorted(model_consistency_scores, key=lambda x: x[0])
    # Calculate number of models with >25% consistency and average consistency level within these
    # combination of these gives best indication of how well the interface is conserved across models - allows differentiation between all models having ok consistency vs most models are very good but a few are bad
    num_consistent_models = sum(1 for _, score in sorted_model_consistency_scores if score > 0.25)
    average_consistency_level = sum(score for _, score in sorted_model_consistency_scores if score > 0.25) / num_consistent_models if num_consistent_models else 0
    
    # Sort consistent pairs and convert back to chain specific numbering
    inv_abs_res_lookup_dict = {v: k for k, v in abs_res_lookup_dict.items()}
    res_pairs = []
    for abs_res1, abs_res2 in consistent_res_pairs:
        if abs_res1 in inv_abs_res_lookup_dict and abs_res2 in inv_abs_res_lookup_dict:
            chain1_id, res1_id = inv_abs_res_lookup_dict[abs_res1]
            chain2_id, res2_id = inv_abs_res_lookup_dict[abs_res2]
            res_pairs.append((f"{chain1_id} {res1_id}", f"{chain2_id} {res2_id}"))
        else:
            # Handle cases where a mapping might not be found
            print(f"Mapping not found for abs_res_ids: {abs_res1}, {abs_res2}")
    return sorted_model_consistency_scores, num_consistent_models, average_consistency_level, res_pairs

def measure_repeatability(folder_path, distance_cutoff=10.0, pae_cutoff=15.0, all_atom=False):
    """
    Measure the similarity of confident interfaces across multiple models.

    Outputs the number of models with >25% consistency in interface residue positions
    and the average consistency level of these models - combining these gives the best
    indication of how well the interface is conserved across models, as it allows
    differentiation between all models having some consistency vs most models being very
    good but a few being bad.

    Parameters:
        - folder_path (str): Path to the folder containing the structure files
        - distance_cutoff (float): Maximum distance between CA atoms of residues (in
          different chains) to be considered as interface
        - pae_cutoff (float): Maximum PAE value between a pair of residues for them
          to be considered as being confidently positioned with relation to each other
        - all_atom (bool): True if all atoms should be considered for measuring interface
          distance, False if only CA atoms should be considered

    Returns:
        - num_consistent (int): Number of models with >25% consistency in interface
          residue positions
        - level_consistent (float): Average consistency level of models with >25%
          consistency
    """
    # Find the top ranked predicted structure file and check if it is in PDB or CIF format
    rank1_structure_file, rank1_json, _, _, _ = find_rank_001_files(folder_path)
    is_pdb = True
    if rank1_structure_file.suffix != '.pdb':
        is_pdb = False
    if not is_pdb and rank1_structure_file.suffix != '.cif':
        raise ValueError("Unsupported file format. Structure files must be in CIF or PDB format.")

    rank_1_structure_model = parse_structure_file(rank1_structure_file, is_pdb)
    chain_residue_map = map_chains_and_residues(rank_1_structure_model)
    abs_res_lookup_dict = {(chain_id, res_id): abs_res_id for chain_id, res_id, _, abs_res_id in chain_residue_map}
    
    # Find confident interface residues in the top ranked model
    confident_pairs = find_confident_interface_residues(rank_1_structure_model, rank1_json, distance_cutoff, pae_cutoff, abs_res_lookup_dict, is_pdb, all_atom)
    
    # If confident interface residue pairs are found, check distances between these across other models, otherwise return 0
    if confident_pairs:
        #print(f"Found {len(confident_pairs)} confident interface residue pairs in the top ranked model.")
        all_scores, num_consistent, level_consistent, consistent_pairs = check_distances_across_models(folder_path, confident_pairs, distance_cutoff, abs_res_lookup_dict, is_pdb, all_atom)
        print(f"Individual model scores: {all_scores}")
        #print(f"Num models >25% consistent: {num_consistent} out of {len(all_scores)} total models.")
        #print(f"Avg consistency for consistent models: {level_consistent:.2f}")
        print(f"Consistent pairs: {consistent_pairs}")

    else:
        all_scores, num_consistent, level_consistent = [], 0, 0
        #print("No confident interface residue pairs found in the top ranked model.")
    return num_consistent, level_consistent

####################################################################################################
# Example usage
# pdb files (AF2)
#folder_path = "data/Sak_Sas6/Sak_D3+Sas6_D1"
#num_consistent, level_consistent = measure_repeatability(folder_path, distance_cutoff=5, pae_cutoff=10, all_atom=True)

# cif files (AF3) - run with same command as above, adjustment is internal
#folder_path = "fold_ana2_flp_sak_fl"
#num_consistent, level_consistent = measure_repeatability(folder_path, distance_cutoff=5, pae_cutoff=15, all_atom=True)

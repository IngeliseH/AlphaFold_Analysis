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
FOR GENERAL USE
get_residue_pairs
find_confident_pairs

FOR COMPARING WITH A PARTICULAR MODEL
calculate_rop_score
calculate_percent_rop

FOR USE WITH SINGLE PREDICTION OR RANK 1 MODEL
check_distances_across_models
measure_repeatibility

TODO: Check compatibility with AF3 or multichain predictions
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
        - distance_cutoff (float): Maximum distance between atoms of residues to
          be considered as interface - recommended to use 5 if all_atom
          is True, 10 if all_atom is False
        - abs_res_lookup_dict (dict): Dictionary mapping chain and residue IDs to absolute
          residue IDs.
        - all_atom (bool): True if all atoms should be considered for measuring interface
          distance, False if only CA atoms should be considered.
    
    Returns:
        - residue_pairs (list): List of unique residue pairs that are within the distance cutoff,
          where each pair is a tuple of residue numbers.
    """
    # Load structure
    structure_model = (parse_structure_file(structure_input)
                       if isinstance(structure_input, (str, Path))
                       else structure_input)
    if not isinstance(structure_model, Model.Model):
        raise TypeError("structure_input must be a valid file path or protein model object.")

    # Get residue absolute IDs and coordinates
    chain_residues = {}
    for chain in structure_model.get_chains():
        residues = []
        for res in chain.get_residues():
            abs_res_id = abs_res_lookup_dict.get((chain.id, res.id[1]))
            if abs_res_id is not None:
                if all_atom:
                    atom_coords = np.array([atom.coord for atom in res.get_atoms()])
                    residues.append((abs_res_id, atom_coords))
                elif 'CA' in res:
                    residues.append((abs_res_id, res['CA'].coord))
        chain_residues[chain.id] = residues

    unique_residue_pairs = set()  # Track unique residue pairs
    chain_ids = list(chain_residues.keys())

    # Iterate over chain pairs
    for i, chain1_id in enumerate(chain_ids[:-1]):
        residues1 = chain_residues[chain1_id]
        for chain2_id in chain_ids[i+1:]:
            residues2 = chain_residues[chain2_id]
            if not residues1 or not residues2:
                continue

            if all_atom:
                # Concatenate all atom coordinates for each chain
                coords1 = np.concatenate([r[1] for r in residues1])
                coords2 = np.concatenate([r[1] for r in residues2])
                res_map1 = np.repeat([r[0] for r in residues1], [len(r[1]) for r in residues1])
                res_map2 = np.repeat([r[0] for r in residues2], [len(r[1]) for r in residues2])
            else:
                # Use CA coordinates
                coords1 = np.array([r[1] for r in residues1])
                coords2 = np.array([r[1] for r in residues2])
                res_map1 = np.array([r[0] for r in residues1])
                res_map2 = np.array([r[0] for r in residues2])

            # Use KD-tree for efficient distance calculation
            tree2 = cKDTree(coords2)
            neighbors = tree2.query_ball_point(coords1, r=distance_cutoff)

            for idx1, indices in enumerate(neighbors):
                abs_res_id1 = res_map1[idx1]
                for idx2 in indices:
                    abs_res_id2 = res_map2[idx2]

                    # Add the pair in sorted order to ensure uniqueness
                    pair = tuple(sorted((abs_res_id1, abs_res_id2)))

                    # Skip if the pair has already been identified
                    if pair in unique_residue_pairs:
                        continue

                    # Add the new unique pair
                    unique_residue_pairs.add(pair)

    # Convert the set of unique pairs to a list
    residue_pairs = list(unique_residue_pairs)
    return residue_pairs

def find_confident_pairs(residue_pairs, pae_data, pae_cutoff):
    """
    Find pairs of interface residues that are confidently predicted to be close (in the
    top ranked model).
    
    Parameters:
        - residue_pairs (list): List of residue pairs that are within the distance cutoff,
          where each pair is a tuple of residue numbers.
        - pae_data (dict): Dictionary containing the PAE values for each residue pair.
        - pae_cutoff (float): Maximum PAE value between a pair of residues for them
          to be considered as being confidently positioned with relation to each other.
        
    Returns:
        - confident_pairs (list): List of residue pairs that are confidently predicted
          to be close in the top ranked model, where each pair is a tuple of residue
          numbers relative to the full sequence (i.e., not just position within the chain).
    """
    # Check PAE values for each nearby residue pair, return pairs with PAE < cutoff
    confident_pairs = set()
    for res1, res2 in residue_pairs:
        if pae_data[res1][res2] < pae_cutoff or pae_data[res2][res1] < pae_cutoff:
            confident_pairs.add((res1, res2))
    return confident_pairs

# Functions for calculating rop of a specfic rank model, not assumed to be rank 1
def calculate_rop_score(residue_pairs, other_models_residue_pairs, threshold = 0.25):
    """
    Calculates the ROP score for a set of residue pairs in a model, based on the
    percentage of these pairs that are also present in other models.

    Parameters:
        - residue_pairs (set): Set of residue pairs
        - other_models_residue_pairs (list): List of sets of residue pairs from other models
          (should NOT include the model being compared)
        - threshold (float): Minimum percentage of residue pairs that must be present in
          each other model for the interface to be considered consistent
    Returns:
        - rop_score (int): Number of models where the percentage of residue pairs from the
          model that are also present in the other model is greater than the threshold.
    """
    if not residue_pairs or not other_models_residue_pairs:
        return 0

    if not isinstance(residue_pairs, set):
        residue_pairs = set(residue_pairs)

    if not all(isinstance(other_pairs, set) for other_pairs in other_models_residue_pairs):
        other_models_residue_pairs = [set(other_pairs) for other_pairs in other_models_residue_pairs]

    rop_score = 0
    for other_model in other_models_residue_pairs:
        # Compute the percentage of model's confident_pairs present in other_model's residue_pairs
        shared_pairs = residue_pairs & other_model
        percentage = len(shared_pairs) / len(residue_pairs)
        if percentage >= threshold:
            rop_score += 1
    return rop_score

def calculate_percent_rop(residue_pairs, other_models_residue_pairs, threshold = 0.25):
    """
    Calculates percent_rop for a set of residue pairs in a model, based on the
    percentage of these pairs that are also present in other models.

    Parameters:
        - residue_pairs (set): Set of residue pairs in model of interest
        - other_models_residue_pairs (list): List of sets of residue pairs from other models
          (should NOT include the model being compared)
        - threshold (float): Minimum percentage of residue pairs that must be present in
          each other model for the interface to be considered

    Returns:
        - percent_rop (float): Average percentage of residue pairs from the model that are
          also present in other models surpassing threshold
    """
    if not residue_pairs or not other_models_residue_pairs:
        return 0
    percentages = []
    for other_model in other_models_residue_pairs:
        shared_pairs = residue_pairs & other_model
        percentage = len(shared_pairs) / len(residue_pairs) if residue_pairs else 0
        if percentage >= threshold:
            percentages.append(percentage)
    return np.mean(percentages) if percentages else 0

# Functions for use when only looking at rank 1 model, or for quickly checking a single prediction
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
          different chains) to be considered as interface - recommended to use 5 if all_atom
          is True, 10 if all_atom is False
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
        model_res_pairs = set(get_residue_pairs(structure_model, distance_cutoff, abs_res_lookup_dict, all_atom))
        # Calculate proportion of confident pairs from top model which are also close in current model
        score = sum(1 for (res1, res2) in confident_pairs if (res1, res2) in model_res_pairs) / len(confident_pairs) if confident_pairs else 0
        if score > 0.25:
            consistent_res_pairs.intersection_update(model_res_pairs)
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

def measure_repeatability(folder_path, distance_cutoff=5.0, pae_cutoff=15.0, all_atom=True):
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
          different chains) to be considered as interface - recommended to use 5 if all_atom
          is True, 10 if all_atom is False
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
    
    # Find interface residue pairs in the top ranked model
    # Load PAE data
    pae_data = extract_pae(rank1_json)
    # If cif, condense pae data in case of phosphorylated residues  
    # necessary even without modfication, but unsure why
    if not is_pdb:
        pae_data = correct_cif_pae(rank_1_structure_model, pae_data)

    # Find residue pairs within the distance cutoff in the top ranked model
    residue_pairs = get_residue_pairs(rank_1_structure_model, distance_cutoff, abs_res_lookup_dict, all_atom)
    # Find confident subset of these
    confident_pairs = find_confident_pairs(residue_pairs, pae_data, pae_cutoff)
    
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

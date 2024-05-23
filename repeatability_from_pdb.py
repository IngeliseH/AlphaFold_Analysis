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
from analysis_utility import extract_pae, find_rank_001_files, parse_structure_file, map_chains_and_residues
from phosphorylation_handling import correct_cif_pae

def get_residue_pairs(structure_model, distance_cutoff, abs_res_lookup_dict, all_atom):
    """
    Identify pairs of residues across different chains that are within a specified
    distance cutoff.
    
    Parameters:
        - pdb_path (str): Path to the PDB or CIF file
        - distance_cutoff (float): Maximum distance between CA atoms of residues to
          be considered as interface
        - abs_res_lookup_dict (dict): Dictionary mapping chain and residue IDs to absolute
          residue IDs
        - all_atom (bool): True if all atoms should be considered for measuring interface
          distance, False if only CA atoms should be considered
    
    Returns:
        - residue_pairs (list): List of residue pairs that are within the distance cutoff,
          where each pair is a tuple of residue numbers)
    """
    chains = list(structure_model.get_chains())
    # Get mapping of residues to absolute residue ids
    def get_abs_res_id(abs_res_lookup_dict, target_chain_id, target_res_id):
        return abs_res_lookup_dict.get((target_chain_id, target_res_id), None)  # Returns None if the key is not found

    # Find residue pairs across chains closer than distance cutoff
    residue_pairs = []
    nearby_threshold = distance_cutoff + 10  # Initial broader threshold

    # iterate through each residue in each chain, and compare to all residues in other chains
    for i, chain1 in enumerate(chains[:-1]):
        for chain2 in chains[i+1:]:
            for res1 in chain1.get_residues():
                res1_id = get_abs_res_id(abs_res_lookup_dict, chain1.id, res1.id[1])
                for res2 in chain2.get_residues():
                    # Calculate residue number relative to the full sequence (assuming two chains only)
                    res2_id = get_abs_res_id(abs_res_lookup_dict, chain2.id, res2.id[1])
                    if 'CA' in res1 and 'CA' in res2:
                        distance = np.linalg.norm(res1['CA'].coord - res2['CA'].coord)
                        if not all_atom and distance < distance_cutoff:
                            residue_pairs.append((res1_id, res2_id))
                        elif all_atom and distance < nearby_threshold:
                            found_interaction = False
                            for atom1 in res1:
                                if found_interaction:
                                    break
                                for atom2 in res2:
                                    distance = np.linalg.norm(atom1.coord - atom2.coord)
                                    if distance < distance_cutoff:
                                        residue_pairs.append((res1_id, res2_id))
                                        found_interaction = True
                                        break
    return residue_pairs


def find_confident_interface_residues(structure_model, pae_path, distance_cutoff, pae_cutoff, abs_res_lookup_dict, is_pdb, all_atom):
    """
    Find pairs of interface residues that are confidently predicted to be close (in the
    top ranked model).
    
    Parameters:
        - structure_model (Bio.PDB.Model): Structure model of the top ranked predicted
          structure
        - pae_path (str): Path to the JSON file containing PAE data
        - distance_cutoff (float): Maximum distance between CA atoms of residues to
          be considered as interface
        - pae_cutoff (float): Maximum PAE value between a pair of residues for them
          to be considered as being confidently positioned with relation to each other
        - abs_res_lookup_dict (dict): Dictionary mapping chain and residue IDs to absolute
          residue IDs
        - is_pdb (bool): True if the file is in PDB format, False if in CIF format
        - all_atom (bool): True if all atoms should be considered for measuring interface
          distance, False if only CA atoms should be considered
        
    Returns:
        - confident_pairs (list): List of residue pairs that are confidently predicted
          to be close in the top ranked model, where each pair is a tuple of residue
          numbers relative to the full sequence (ie not just position within the chain)
    """
    # Extract PAE data and find residue pairs within distance cutoff
    pae_data = extract_pae(pae_path)
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
        - num_consistent_models (int): Number of models with >50% consistency in interface
          residue positions
        - average_consistency_level (float): Average consistency level of models with >50%
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
        if score > 0.5:
            model_res_pairs_set = set(model_res_pairs)
            consistent_res_pairs.intersection_update(model_res_pairs_set)
        model_consistency_scores.append((rank_number, score))

    # Sort scores by model rank number
    sorted_model_consistency_scores = sorted(model_consistency_scores, key=lambda x: x[0])
    # Calculate number of models with >50% consistency and average consistency level within these
    # combination of these gives best indication of how well the interface is conserved across models - allows differentiation between all models having ok consistency vs most models are very good but a few are bad
    num_consistent_models = sum(1 for _, score in sorted_model_consistency_scores if score > 0.5)
    average_consistency_level = sum(score for _, score in sorted_model_consistency_scores if score > 0.5) / num_consistent_models if num_consistent_models else 0
    
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

def measure_repeatibility(folder_path, distance_cutoff=10.0, pae_cutoff=15.0, all_atom=False):
    """
    Measure the similarity of confident interfaces across multiple models.

    Outputs the number of models with >50% consistency in interface residue positions
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
        - num_consistent (int): Number of models with >50% consistency in interface
          residue positions
        - level_consistent (float): Average consistency level of models with >50%
          consistency
    """
    # Find the top ranked predicted structure file and check if it is in PDB or CIF format
    rank1_structure_file, rank1_json, _, _ = find_rank_001_files(folder_path)
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
        print(f"Found {len(confident_pairs)} confident interface residue pairs in the top ranked model.")
        all_scores, num_consistent, level_consistent, consistent_pairs = check_distances_across_models(folder_path, confident_pairs, distance_cutoff, abs_res_lookup_dict, is_pdb, all_atom)
        print(f"Individual model scores: {all_scores}")
        print(f"Num models >50% consistent: {num_consistent} out of {len(all_scores)} total models.")
        print(f"Avg consistency for consistent models: {level_consistent:.2f}")
        print(f"Consistent pairs: {consistent_pairs}")

    else:
        all_scores, num_consistent, level_consistent = [], 0, 0
        print("No confident interface residue pairs found in the top ranked model.")
    return num_consistent, level_consistent

####################################################################################################
# Example usage
# pdb files (AF2)
folder_path = "Sak_Sas6/Sak_D3+Sas6_D1"
num_consistent, level_consistent = measure_repeatibility(folder_path, distance_cutoff=5, pae_cutoff=10, all_atom=True)

# cif files (AF3) - run with same command as above, adjustment is internal
folder_path = "fold_ana2_flp_sak_fl"
num_consistent, level_consistent = measure_repeatibility(folder_path, distance_cutoff=5, pae_cutoff=15, all_atom=True)

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
import os
import re
from Bio.PDB import PDBParser
from Bio.PDB import MMCIFParser
import numpy as np
from analysis_utility import extract_pae, find_rank_001_files, parse_cif
from phosphorylation_handling import map_res_to_pos, condense_pae

def get_residue_pairs(pdb_path, distance_cutoff, is_pdb):
    """
    Identify pairs of residues across different chains that are within a specified
    distance cutoff.
    
    Parameters:
        - pdb_path (str): Path to the PDB or CIF file
        - distance_cutoff (float): Maximum distance between CA atoms of residues to
          be considered as interface
        - is_pdb (bool): True if the file is in PDB format, False if in CIF format
    
    Returns:
        - residue_pairs (list): List of residue pairs that are within the distance cutoff,
          where each pair is a tuple of residue numbers)
    """
    # Parse the structure
    parser = PDBParser() if is_pdb else MMCIFParser()
    structure = parser.get_structure('protein', pdb_path)
    model = structure[0]

    # Find chain lengths
    chains = list(model.get_chains())
    chain_lengths = {}
    for chain in chains:
        count = 0
        for _ in chain.get_residues():
            count += 1
        chain_lengths[chain.id] = count
    
    # Find residue pairs across chains closer than distance cutoff
    residue_pairs = []
    for i, chain1 in enumerate(chains[:-1]):
        for chain2 in chains[i+1:]:
            for res1 in chain1.get_residues():
                res1_id = res1.id[1]
                for res2 in chain2.get_residues():
                    try:
                        distance = res1['CA'].coord - res2['CA'].coord
                        if np.linalg.norm(distance) < distance_cutoff:
                            res2_id = res2.id[1] + chain_lengths[chain1.id]
                            residue_pairs.append((res1_id, res2_id))
                    except KeyError:
                        # Handle the absence of 'CA' atom
                        continue
    return residue_pairs

def find_confident_interface_residues(pdb_path, pae_path, distance_cutoff, pae_cutoff, is_pdb):
    """
    Find pairs of interface residues that are confidently predicted to be close (in the
    top ranked model).
    
    Parameters:
        - pdb_path (str): Path to the PDB or CIF file
        - pae_path (str): Path to the JSON file containing PAE data
        - distance_cutoff (float): Maximum distance between CA atoms of residues to
          be considered as interface
        - pae_cutoff (float): Maximum PAE value between a pair of residues for them
          to be considered as being confidently positioned with relation to each other
        - is_pdb (bool): True if the file is in PDB format, False if in CIF format
        
    Returns:
        - confident_pairs (list): List of residue pairs that are confidently predicted
          to be close in the top ranked model, where each pair is a tuple of residue
          numbers relative to the full sequence (ie not just position within the chain)
    """
    # Extract PAE data and find residue pairs within distance cutoff
    pae_data = extract_pae(pae_path)
    residue_pairs = get_residue_pairs(pdb_path, distance_cutoff, is_pdb)
    
    # If cif, condense pae data in case of phosphorylated residues  
    # necessary even without modfication, but unsure why
    if not is_pdb:
        chain_residue_map = parse_cif(pdb_path)
        res_to_pos = map_res_to_pos(chain_residue_map)
        pae_data = condense_pae(pae_data, res_to_pos)

    # Check PAE values for each nearby residue pair, return pairs with PAE < cutoff
    confident_pairs = []
    for res1, res2 in residue_pairs:
        pae_value = pae_data[res1-1][res2-1]
        if pae_value < pae_cutoff:
            confident_pairs.append((res1, res2))
    return confident_pairs

def check_distances_across_models(folder_path, confident_pairs, distance_cutoff, is_pdb):
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
        - is_pdb (bool): True if the files are in PDB format, False if in CIF format
    """
    # Find all structure files except the rank 1/top model
    file_extension = '.pdb' if is_pdb else '.cif'
    top_model_patterns = ['rank_001', 'rank_1'] if is_pdb else ['model_0']
    rank_pattern = r'rank_(\d+)' if is_pdb else r'model_(\d+)'
    files = [f for f in os.listdir(folder_path) if f.endswith(file_extension) and not any(top_model_pattern in f for top_model_pattern in top_model_patterns)]
    structure_files = sorted(files, key=lambda x: int(re.search(rank_pattern, x).group(1)))

    # Identify residue pairs within distance cutoff for each model
    model_scores = []
    for structure_file in structure_files:
        structure_path = os.path.join(folder_path, structure_file)
        rank_number = int(re.search(rank_pattern, structure_file).group(1))
        current_pairs = get_residue_pairs(structure_path, distance_cutoff, is_pdb)
        # Calculate proportion of confident pairs from top model which are also close in current model
        score = sum(1 for (res1, res2) in confident_pairs if (res1, res2) in current_pairs) / len(confident_pairs) if confident_pairs else 0
        model_scores.append((rank_number, score))

    # Sort scores by model rank number
    sorted_scores = sorted(model_scores, key=lambda x: x[0])
    # Calculate number of models with >50% consistency and average consistency level within these
    # combination of these gives best indication of how well the interface is conserved across models - allows differentiation between all models having ok consistency vs most models are very good but a few are bad
    num_consistent = sum(1 for _, score in sorted_scores if score > 0.5)
    level_consistent = sum(score for _, score in sorted_scores if score > 0.5) / num_consistent if num_consistent else 0
    return sorted_scores, num_consistent, level_consistent

def measure_repeatibility(folder_path, distance_cutoff=10.0, pae_cutoff=15.0):
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

    Returns:
        - num_consistent (int): Number of models with >50% consistency in interface
          residue positions
        - level_consistent (float): Average consistency level of models with >50%
          consistency
    """
    # Find the top ranked predicted structure file and check if it is in PDB or CIF format
    rank1_structure, rank1_json, _, _ = find_rank_001_files(folder_path)
    is_pdb = True
    if not rank1_structure.endswith('.pdb'):
        is_pdb = False
    if not is_pdb and not rank1_structure.endswith('.cif'):
        raise ValueError("Unsupported file format. Structure files must be in CIF or PDB format.")
    
    # Find confident interface residues in the top ranked model
    confident_pairs = find_confident_interface_residues(rank1_structure, rank1_json, distance_cutoff, pae_cutoff, is_pdb)
    
    # If confident interface residue pairs are found, check distances between these across other models, otherwise return 0
    if confident_pairs:
        print(f"Found {len(confident_pairs)} confident interface residue pairs in the top ranked model.")
        all_scores, num_consistent, level_consistent = check_distances_across_models(folder_path, confident_pairs, distance_cutoff, is_pdb)
        print(f"Individual model scores: {all_scores}")
        print(f"Number of models with >50% consistency: {num_consistent} out of {len(all_scores)} total models.")
        print(f"Average consistency level: {level_consistent:.2f}")
    else:
        all_scores, num_consistent, level_consistent = [], 0, 0
        print("No confident interface residue pairs found in the top ranked model.")
    return num_consistent, level_consistent

####################################################################################################
# Example usage
# pdb files (AF2)
folder_path = "Sak_Sas6/Sak_D3+Sas6_D1"
num_consistent, level_consistent = measure_repeatibility(folder_path)

# cif files (AF3) - run with same command as above, adjustment is internal
folder_path = "fold_ana2_flp_sak_fl"
num_consistent, level_consistent = measure_repeatibility(folder_path, distance_cutoff=10, pae_cutoff=10)

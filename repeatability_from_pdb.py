import os
import re
from Bio.PDB import PDBParser
import numpy as np
from analysis_utility import extract_pae

def get_residue_pairs(pdb_path, distance_cutoff):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_path)
    model = structure[0]  # Assuming only one model per file

    # Find residue pairs across chains with distance < cutoff
    chains = list(model.get_chains())
    residue_pairs = []
    for i, chain1 in enumerate(chains[:-1]):
        for chain2 in chains[i+1:]:
            for res1 in chain1.get_residues():
                for res2 in chain2.get_residues():
                    try:
                        distance = res1['CA'].coord - res2['CA'].coord
                        if np.linalg.norm(distance) < distance_cutoff:
                            residue_pairs.append((res1, res2))
                    except KeyError:
                        # Handle the absence of 'CA' atom
                        continue
    return residue_pairs

def find_confident_interface_residues(pdb_path, pae_path, distance_cutoff, pae_cutoff):
    residue_pairs = get_residue_pairs(pdb_path, distance_cutoff)
    pae_data = extract_pae(pae_path)
    confident_pairs = []
    for res1, res2 in residue_pairs:
        pae_value = pae_data[res1.id[1]-1][res2.id[1]-1]
        if pae_value < pae_cutoff:
            confident_pairs.append((res1, res2))
    return confident_pairs

def check_distances_across_models(folder_path, confident_pairs, distance_cutoff):
    # Find all pdb files except the rank 1
    pdb_files = sorted([f for f in os.listdir(folder_path) if f.endswith('.pdb') and 'rank_001' not in f],
                       key=lambda x: int(re.search(r'rank_(\d+)', x).group(1)))

    model_scores = []
    for pdb_file in pdb_files:
        rank_number = int(re.search(r'rank_(\d+)', pdb_file).group(1))
        pdb_path = os.path.join(folder_path, pdb_file)
        current_pairs = get_residue_pairs(pdb_path, distance_cutoff)
        # Calculate score for the current model
        score = sum(1 for res1, res2 in confident_pairs if (res1, res2) in current_pairs) / len(confident_pairs) if confident_pairs else 0
        model_scores.append((rank_number, score))

    # Sort scores by rank number and calculate the average
    sorted_scores = sorted(model_scores, key=lambda x: x[0])
    average_score = sum(score for _, score in sorted_scores) / len(sorted_scores) if sorted_scores else 0
    return sorted_scores, average_score

# Example usage
folder_path = "Sak_Sas6/Sak_D3+Sas6_D1"
rank1_pdb_name = next(f for f in os.listdir(folder_path) if 'rank_001' in f and f.endswith('.pdb'))
rank1_pdb = os.path.join(folder_path, rank1_pdb_name)
rank1_pae_name = next(f for f in os.listdir(folder_path) if 'rank_001' in f and f.endswith('.json'))
rank1_pae = os.path.join(folder_path, rank1_pae_name)
distance_cutoff = 10.0  # Example cutoff in Angstroms
pae_cutoff = 15.0       # Example PAE cutoff

confident_pairs = find_confident_interface_residues(rank1_pdb, rank1_pae, distance_cutoff, pae_cutoff)
print(f"Number of confident pairs in rank one model: {len(confident_pairs)}. Checking consistency across models...")
model_scores, average_score = check_distances_across_models(folder_path, confident_pairs, distance_cutoff)
print(f"Individual model scores: {model_scores}")
print(f"Average score across models: {average_score:.2f}")
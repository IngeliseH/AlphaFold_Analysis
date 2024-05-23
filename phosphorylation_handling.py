"""
Convert phosphorylated residue indices to regular residue indices in the PAE matrix.

Functions:
map_res_to_pos(chain_residue_map)
condense_pae(pae_matrix, res_to_pos)
"""
import numpy as np
from analysis_utility import map_chains_and_residues

# create mapping from residue indices to PAE matrix indices
def map_res_to_pos(chain_residue_map):
    """
    Creates a mapping from residue indices to PAE matrix indices.

    Parameters:
        - chain_residue_map (list): List of tuples representing chain, residue ID, residue name, and absolute residue ID.

    Returns:
        - dict: Mapping from residue indices to PAE matrix indices.
    """
    res_to_pos = {}
    prev_pos = 0
    for _, _, res_name, abs_res_id in chain_residue_map:
        if res_name == "TPO":
            num_atoms = 10
        elif res_name == "SEP":
            num_atoms = 9
        else:
            num_atoms = 1
        for _ in range(num_atoms):
            # append prev_pos to list of items associated with abs_res_id
            if abs_res_id in res_to_pos:
                res_to_pos[abs_res_id].append(prev_pos)
            else:
                res_to_pos[abs_res_id] = [prev_pos]
            prev_pos += 1
    return res_to_pos

def condense_pae(pae_matrix, res_to_pos):
    """
    Condenses the PAE matrix by combining data for phosphorylated residues.

    Parameters:
        - pae_matrix (list): The original PAE matrix.
        - res_to_pos (dict): Mapping from residue indices to PAE matrix indices.

    Returns:
        - list: Condensed PAE matrix.
    """
    # Convert the PAE matrix to a numpy array for easier manipulation
    pae_matrix = np.array(pae_matrix)
    
    # Determine the number of residues in the condensed matrix
    num_residues = max(res_to_pos.keys()) + 1

    # Create a new PAE matrix
    new_pae_matrix = np.full((num_residues, num_residues), np.inf)

    for i in range(num_residues):
        for j in range(num_residues):
            i_res_indices = res_to_pos[i]
            j_res_indices = res_to_pos[j]

            min_value = np.min([pae_matrix[ii, jj] for ii in i_res_indices for jj in j_res_indices])
            new_pae_matrix[i, j] = min_value

    return new_pae_matrix

def correct_cif_pae(model, pae_data):
    """
    Overall function to correct the PAE matrix for phosphorylated residues.
    
    Parameters:
        - model (Bio.PDB.Model.Model): The structure model.
        - pae_data (list): The original PAE matrix.
        
    Returns:
        - list: Condensed PAE matrix.
    """
    chain_residue_map = map_chains_and_residues(model)
    res_to_pos = map_res_to_pos(chain_residue_map)
    condensed_pae_matrix = condense_pae(pae_data, res_to_pos)

    return condensed_pae_matrix

####################################################################################################
# Example usage
#from analysis_utility import parse_structure_file, extract_pae, find_rank_001_files
#folder_path = "fold_ana2_flp_sak_fl"

#cif_file, json_file, _, _ = find_rank_001_files(folder_path)
#structure_model = parse_structure_file(cif_file, is_pdb=False)
#pae_data = extract_pae(json_file)
#corrected_pae = correct_cif_pae(structure_model, pae_data)

"""
Convert phosphorylated residue indices to regular residue indices in the PAE matrix.

Functions:
map_res_to_pos(chain_residue_map)
condense_pae(pae_matrix, res_to_pos)
"""
import numpy as np

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
    for chain_id, res_id, res_name, abs_res_id in chain_residue_map:
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

"""
Process PAE matrix

Functions:
find_min_pae
visualize_pae_matrix
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from analysis_utility import determine_chain_lengths

def find_min_pae(structure_model, pae_matrix):
    """
    Identify the minimum interprotein PAE values between all pairs of proteins (chains) in the model.

    Parameters:
        - structure_model (Bio.PDB.Model): The structure model of the protein
        - pae_matrix (np.array): The PAE matrix of the protein

    Returns:
        - tuple or dict: If only two chains, returns a tuple (min PAE, position).
          Otherwise, returns a dictionary of minimum interprotein PAE values and their positions for each chain pair.
    """
    chain_lengths = determine_chain_lengths(structure_model)
    if len(chain_lengths) < 2:
        raise ValueError("This function requires at least two chains for interprotein PAE calculations.")
    
    min_pae_results = {}
    cumulative_lengths = [0] + list(np.cumsum(chain_lengths))
    
    # Iterate over each pair of chains
    for i in range(len(chain_lengths)):
        for j in range(i + 1, len(chain_lengths)):
            start_i = cumulative_lengths[i]
            end_i = cumulative_lengths[i + 1]
            start_j = cumulative_lengths[j]
            end_j = cumulative_lengths[j + 1]

            # Extract PAE between chain i and j
            interprotein_pae1 = pae_matrix[start_i:end_i, start_j:end_j]
            interprotein_pae2 = pae_matrix[start_j:end_j, start_i:end_i]
            min_pae1 = np.min(interprotein_pae1)
            min_pos1 = np.unravel_index(np.argmin(interprotein_pae1), interprotein_pae1.shape)

            min_pae2 = np.min(interprotein_pae2)
            min_pos2 = np.unravel_index(np.argmin(interprotein_pae2), interprotein_pae2.shape)

            if min_pae1 < min_pae2:
                min_pae = min_pae1
                min_pos = min_pos1
            else:
                min_pae = min_pae2
                min_pos = min_pos2

            # Store result in a dictionary
            min_pae_results[f"Chain {i+1} to Chain {j+1}"] = {
                'Min PAE': min_pae,
                'Position': min_pos
            }

    # Simplify output if there are only two chains
    if len(chain_lengths) == 2:
        only_key = next(iter(min_pae_results))  # There will be only one key
        return min_pae_results[only_key]['Min PAE'], min_pae_results[only_key]['Position']

    return min_pae_results


def visualize_pae_matrix(pae_matrix, chain_lengths=None, cmap="Spectral", interval=50, pae_cutoff = 30):
    """
    Visualizes the PAE matrix using a heatmap.
    Recommend cmap Spectral for nice visualization or bwr for match with automatic colormap
    
    Parameters:
        - pae_matrix (np.ndarray): A numpy array representing the PAE matrix.
        - chain_lengths (list, optional): A list of chain lengths.
        - cmap (str): Colormap for the heatmap. Default is "Spectral".
        - interval (int): Interval for marking positions. Default is 50.
        - pae_cutoff (int): Maximum PAE value to display in colour key. Default is 30.

    Returns:
        - None, but displays the heatmap.
    """
    def get_residue_positions(chain_lengths, interval=50):
        """
        For custom axis labels. Returns a list of residue positions within each chain.
        Marks positions at the given interval.

        Parameters:
            - chain_lengths (list): List of integers representing the lengths of each chain.
            - interval (int): Interval for marking positions. Default is 50.

        Returns:
            - list: A list of residue positions for custom axis labels, marked at the given interval.
        """
        residue_positions = []
        current_position = 0

        for chain_length in chain_lengths:
            for i in range(1, chain_length + 1):
                current_position += 1
                if i % interval == 0 or i == 1:
                    residue_positions.append(i)
                else:
                    residue_positions.append("")

        return residue_positions

    # Convert PAE matrix to DataFrame for seaborn heatmap
    pae_df = pd.DataFrame(pae_matrix)
    
    # Plotting
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(pae_df, cmap=cmap, cbar_kws={'label': 'PAE Value'},  vmin=0, vmax=pae_cutoff)
    
    # Add chain boundary lines if chain lengths are provided
    if chain_lengths:
        cumulative_lengths = np.cumsum(chain_lengths)
        for length in cumulative_lengths[:-1]:  # Skip the last one
            ax.axvline(x=length, color='black', linewidth=1)
            ax.axhline(y=length, color='black', linewidth=1)
    
    # Set custom tick labels if chain lengths are provided
    if chain_lengths:
        residue_positions = get_residue_positions(chain_lengths, interval)
        ax.set_xticks(np.arange(len(residue_positions)) + 0.5)
        ax.set_xticklabels(residue_positions, rotation=90)
        ax.set_yticks(np.arange(len(residue_positions)) + 0.5)
        ax.set_yticklabels(residue_positions, rotation=0)
    
    # Plot styling
    plt.title("PAE Matrix Heatmap")
    plt.xlabel("Residue Index")
    plt.ylabel("Residue Index")
    plt.tight_layout()

    plt.show()

####################################################################################################
# Example usage:
#from analysis_utility import extract_pae, find_rank_001_files, parse_structure_file, determine_chain_lengths

# for AF2 files
#folder_path = 'Sak_Sas6/Sak_D3+Sas6_D1'
#pdb_file, json_file, PAE_png, fasta_file = find_rank_001_files(folder_path)
#structure_model = parse_structure_file(pdb_file)
#pae_matrix = extract_pae(json_file)
#min_pae = find_min_pae(structure_model, pae_matrix)
#chain_lengths = determine_chain_lengths(structure_model)
#visualize_pae_matrix(pae_matrix, chain_lengths)

# for AF3 files
#from phosphorylation_handling import correct_cif_pae
#folder_path = 'fold_ana2_flp_sak_fl'
#folder_path = 'fold_ana2p_cterm_dimer_sas6_dimer'
#cif_file, json_file, PAE_png, fasta_file = find_rank_001_files(folder_path)
#structure_model = parse_structure_file(cif_file, is_pdb=False)
#initial_pae_matrix = extract_pae(json_file)
#print(len(initial_pae_matrix))
#pae_matrix = correct_cif_pae(structure_model, initial_pae_matrix)
#print(len(pae_matrix))
#min_pae = find_min_pae(structure_model, pae_matrix)
#chain_lengths = determine_chain_lengths(structure_model)
#visualize_pae_matrix(pae_matrix, chain_lengths)
#visualize_pae_matrix(pae_matrix[400:440, 400:440])

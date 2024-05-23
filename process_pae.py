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
    Identify the minimum interprotein PAE value between two proteins

    Parameters:
        - structure_model (Bio.PDB.Model): The structure model of the protein
        - pae_matrix (np.array): The PAE matrix of the protein

    Returns:
        - float: The minimum interprotein PAE value
    """
    chain_lengths = determine_chain_lengths(structure_model)
    if len(chain_lengths) != 2:
        raise ValueError("This function currently only works for interprotein PAE calculations between 2 proteins")
    len_A, _ = chain_lengths
    interprotein_pae1 = pae_matrix[:len_A, len_A:].T
    interprotein_pae2 = pae_matrix[len_A:, :len_A].T

    min_pae1 = np.min(interprotein_pae1)
    min_pae2 = np.min(interprotein_pae2)
    return min(min_pae1,min_pae2)

def visualize_pae_matrix(pae_matrix, chain_lengths=None, cmap="Spectral", interval=50):
    """
    Visualizes the PAE matrix using a heatmap.
    Recommend cmap Spectral for nice visualization or bwr for match with automatic colormap
    
    Parameters:
    pae_matrix (np.ndarray): A numpy array representing the PAE matrix.
    chain_lengths (list, optional): A list of chain lengths.
    interval (int): Interval for marking positions. Default is 50.
    """
    def get_residue_positions(chain_lengths, interval=50):
        """
        Returns a list of residue positions within each chain for custom axis labels.
        Marks positions at the given interval.

        Parameters:
        chain_lengths (list): List of integers representing the lengths of each chain.
        interval (int): Interval for marking positions. Default is 50.

        Returns:
        list: A list of residue positions for custom axis labels, marked at the given interval.
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
    ax = sns.heatmap(pae_df, cmap=cmap, cbar_kws={'label': 'PAE Value'})
    
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

#folder_path = 'Sak_Sas6/Sak_D3+Sas6_D1'
#pdb_file, json_file, PAE_png, fasta_file = find_rank_001_files(folder_path)
#structure_model = parse_structure_file(pdb_file)
#pae_matrix = extract_pae(json_file)
#min_pae = find_min_pae(structure_model, pae_matrix)
#chain_lengths = determine_chain_lengths(structure_model)
#visualize_pae_matrix(pae_matrix, chain_lengths)

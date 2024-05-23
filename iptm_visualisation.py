"""
Functions to extract and visualize the highest prediction iptm values for each domain pair of a protein pair

Only for use with predictions outputting log files - doesn't work with AlphaFold3

Functions:
extract_highest_iptm_value
sort_domains
create_iptm_matrix
visualize_iptm_matrix
"""
import os
import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns

def extract_highest_iptm_value(log_file_path):
    """
    Extract the highest iptm value from a log.txt file.
    
    Parameters:
        - log_file_path (str): The path to the log.txt file.
        
    Returns:
        - float: The highest iptm value.
    """
    with open(log_file_path, 'r') as file:
        log_content = file.read()
    # for colabfold < 1.5.5
    iptm_values = re.findall(r'iptm (\d*\.?\d+)', log_content)
    # for colabfold 1.5.5
    if not iptm_values:
        iptm_values = re.findall(r'ipTM=(\d*\.?\d+)\n', log_content)
    iptm_values = [float(value) for value in iptm_values]
    return max(iptm_values) if iptm_values else 0

def sort_domains(domain_list):
    """
    Sort domains based on the numerical suffix.
    
    Parameters:
        - domain_list (list): A list of domain names.
    
    Returns:
        - list: The sorted list of domain names.
    """
    return sorted(domain_list, key=lambda x: int(re.search(r'(\d+)$', x).group()))

def create_iptm_matrix(base_folder):
    """
    Create a matrix of the highest iptm values for each domain pair.
    
    Parameters:
        - base_folder (str): The path to the folder containing the log.txt files.
        
    Returns:
        - pd.DataFrame: A DataFrame containing the iptm values for each domain pair.
    """
    domain_pairs = {}
    protein1_domains = set()
    protein2_domains = set()

    for root, _, files in os.walk(base_folder):
        if 'log.txt' in files:
            folder_name = os.path.basename(root)
            protein1_domain, protein2_domain = folder_name.split('+')
            protein1_domains.add(protein1_domain)
            protein2_domains.add(protein2_domain)
            highest_iptm = extract_highest_iptm_value(os.path.join(root, 'log.txt'))
            domain_pairs[(protein1_domain, protein2_domain)] = highest_iptm

    # Sort the domain lists
    protein1_domains = sort_domains(list(protein1_domains))
    protein2_domains = sort_domains(list(protein2_domains))

    matrix = pd.DataFrame(index=protein2_domains, columns=protein1_domains).fillna(0)
    for (protein1_domain, protein2_domain), iptm_value in domain_pairs.items():
        matrix.at[protein2_domain, protein1_domain] = iptm_value

    return matrix

def visualize_iptm_matrix(matrix, output_png_path):
    """
    Visualize the iptm matrix and save as a PNG file.
    
    Parameters:
        - matrix (pd.DataFrame): The iptm matrix.
        - output_png_path (str): The path to save the PNG file.
    
    Returns:
        - None, but saves the PNG file.
    """
    plt.figure(figsize=(10, 8))
    sns.heatmap(matrix, annot=True, fmt=".2f", cmap='coolwarm', vmin=0, vmax=1)
    plt.title('iptm Matrix')
    plt.savefig(output_png_path, dpi=300)

####################################################################################################
# Example usage
#full_path = '/path/to/your/colabfold/output'
#iptm_matrix = create_iptm_matrix(full_path)
#png_file_path = os.path.join(full_path, 'iptm_matrix.png')
#visualize_iptm_matrix(iptm_matrix, png_file_path)

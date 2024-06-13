"""
Functions to extract and visualize top model iptm values

All functons except extract_iptm are desiged only for use with AF2 - as these work
with protein pair folders containing domain pair folders, and AF3 is expected to
only use FL protein pairs

Functions:
extract_iptm
sort_domains
create_iptm_matrix
visualize_iptm_matrix
"""
import os
import json
import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns


def extract_iptm(file_path):
    """
    Extract the highest IPTM value from a log file (.txt) or the ipTM value from a .json
    file based on the file type.

    Parameters:
        - file_path (str): Path to the log or json file from which the ipTM value
          is to be extracted.

    Returns:
        - float: The highest ipTM value found in the file, or None if no ipTM values are
          found.
    """
    # If the file is a txt file, attempt to extract the IPTM value
    if file_path.suffix == '.txt':
        # Regex for extracting ipTM values from log files
        iptm_pattern = re.compile(r'(?:iptm |ipTM=)(\d*\.?\d+)\n')
        max_iptm = None
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    match = iptm_pattern.search(line)
                    if match:
                        iptm_value = float(match.group(1))
                        # Check if it's the highest ipTM value found so far (or the first)
                        if max_iptm is None or iptm_value > max_iptm:
                            max_iptm = iptm_value

        except FileNotFoundError:
            print(f"File not found: {file_path}")
        
        return max_iptm

    # If the file is a JSON file, attempt to extract the IPTM value
    elif file_path.suffix == '.json':
        try:
            with open(file_path, 'r') as file:
                data = json.load(file)
                return float(data["iptm"]) if "iptm" in data else 0
        except (FileNotFoundError, json.JSONDecodeError) as e:
            print(f"Error reading JSON file: {file_path} with error {e}")

    return None

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
#full_path = 'Ana2_mus101'
#iptm_matrix = create_iptm_matrix(full_path)
#png_file_path = os.path.join(full_path, 'iptm_matrix.png')
#visualize_iptm_matrix(iptm_matrix, png_file_path)

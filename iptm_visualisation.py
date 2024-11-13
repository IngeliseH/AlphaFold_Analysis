"""
Functions to extract and visualize top model iptm values

All functions except extract_iptm are desiged only for use with AF2 - as these work
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
from pathlib import Path
import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns
from analysis_utility import find_rank_001_files


def extract_iptm(log_file_path, model_rank='001'):
    """
    Extract the ipTM value for a specific model rank from a log file.

    Parameters:
        - log_file_path (str or Path): Path to the log file from which the ipTM value is to be extracted.
        - model_rank (str): The model rank to extract ipTM value for (e.g., '001').

    Returns:
        - float: The ipTM value for the specified model, or None if not found.
    """
    log_file_path = Path(log_file_path)
    if log_file_path.suffix == '.txt':
        iptm_pattern = re.compile(
            r'rank_(\d+)_.*?pLDDT=\s*(\d*\.\d+)\s*pTM=\s*(\d*\.\d+)\s*ipTM=\s*(\d*\.\d+)'
        )
        try:
            with open(log_file_path, 'r') as file:
                lines = file.readlines()
            iptm_values = {}
            for line in reversed(lines):
                match = iptm_pattern.search(line.strip())
                if match:
                    line_model_rank = match.group(1)
                    iptm_value = float(match.group(4))
                    iptm_values[line_model_rank] = iptm_value
                elif iptm_values:
                    # Break once we've passed the ranking section
                    break
            # Return the ipTM value for the specified model rank
            return iptm_values.get(model_rank, None)
        except FileNotFoundError:
            print(f"File not found: {log_file_path}")
            return None
    
    # If the file is a JSON file, attempt to extract the ipTM value
    elif log_file_path.suffix == '.json':
        try:
            with open(log_file_path, 'r') as file:
                data = json.load(file)
                return float(data["iptm"]) if "iptm" in data else 0
        except (FileNotFoundError, json.JSONDecodeError) as e:
            print(f"Error reading JSON file: {log_file_path} with error {e}")

    else:
        print(f"Unsupported file type: {log_file_path.suffix}")
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
        - base_folder (str): The path to the folder containing the fragment folders.

    Returns:
        - pd.DataFrame: A DataFrame containing the iptm values for each domain pair.
    """
    domain_pairs = {}
    protein1_domains = set()
    protein2_domains = set()

    # Iterate through each folder directly within the base folder
    for fragment_folder in os.listdir(base_folder):
        fragment_folder_path = os.path.join(base_folder, fragment_folder)
        if os.path.isdir(fragment_folder_path):
            # Use find_rank_001_files to locate the log.txt file
            _, _, log_file, _, _ = find_rank_001_files(fragment_folder_path)
            if log_file:
                folder_name = os.path.basename(fragment_folder_path)
                protein1_domain, protein2_domain = folder_name.split('+')
                protein1_domains.add(protein1_domain)
                protein2_domains.add(protein2_domain)
                print(extract_iptm(log_file))
                highest_iptm = extract_iptm(log_file)
                domain_pairs[(protein1_domain, protein2_domain)] = highest_iptm

    # Sort the domain lists
    protein1_domains = sort_domains(list(protein1_domains))
    protein2_domains = sort_domains(list(protein2_domains))

    # Initialize the DataFrame with zeroes
    matrix = pd.DataFrame(index=protein2_domains, columns=protein1_domains).fillna(0)
    
    # Populate the DataFrame with iptm values
    for (protein1_domain, protein2_domain), iptm_value in domain_pairs.items():
        matrix.at[protein2_domain, protein1_domain] = float(iptm_value)

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
full_path = '/Users/poppy/Dropbox/BUB1/BUB1_PLK1'
iptm_matrix = create_iptm_matrix(full_path)
png_file_path = os.path.join(full_path, 'iptm_matrix.png')
visualize_iptm_matrix(iptm_matrix, png_file_path)

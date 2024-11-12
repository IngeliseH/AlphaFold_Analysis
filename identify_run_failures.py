"""
Script to identify missing structure files and analyses in a directory containing protein pair folders.

Functions:
    - find_missing_structure_files
    - find_missing_analyses
"""
import os
import csv
import pandas as pd
from analysis_utility import find_rank_001_files

def find_missing_structure_files(base_folder, **kwargs):
    """
    Identify all protein and domain pair folders where the structure_file is not found.

    Parameters:
        - base_folder (str): Path to the base folder containing all protein pair folders.
        - **kwargs: Additional function options can be specified as keyword arguments.

    Keyword Arguments (kwargs):
        - save_to_csv (bool): If True, saves missing structure files to a CSV file.
        - output_file (str): Name of the CSV file to save missing structure files. Default is "missing_structure_files.csv".

    Returns:
        - missing_files (list of dict): List of dictionaries with details of missing structure files.
    """
    missing_files = []
    save_to_csv = kwargs.get('save_to_csv', False)
    output_file = kwargs.get('output_file', "missing_structure_files.csv")
    headers = ['Protein1', 'Protein2', 'Protein1_Domain', 'Protein2_Domain']

    if save_to_csv:
        output_path = os.path.join(base_folder, output_file)
        file = open(output_path, mode='w', newline='')
        writer = csv.writer(file)
        writer.writerow(headers)  # Write the header row

    # Walk through the base folder containing all protein pair folders
    for protein_pair_folder in os.listdir(base_folder):
        protein_pair_path = os.path.join(base_folder, protein_pair_folder)
        if os.path.isdir(protein_pair_path):
            try:
                # Determine the separator to split the folder name
                if '+' in protein_pair_folder:
                    protein1, protein2 = protein_pair_folder.split('+')
                else:
                    protein1, protein2 = protein_pair_folder.split('_', 1)
                print(f"Checking {protein1} and {protein2}...")

                # Look for domain pair folders either directly within the protein pair folder or within a 'Results' subfolder
                domain_folders = []
                for item in os.listdir(protein_pair_path):
                    item_path = os.path.join(protein_pair_path, item)
                    if os.path.isdir(item_path) and '+' in item:
                        domain_folders.append(item_path)
                    elif item == 'Results' and os.path.isdir(item_path):
                        for sub_item in os.listdir(item_path):
                            sub_item_path = os.path.join(item_path, sub_item)
                            if os.path.isdir(sub_item_path) and '+' in sub_item:
                                domain_folders.append(sub_item_path)

                # Check each domain pair folder identified
                for domain_path in domain_folders:
                    try:
                        domain_pair = os.path.basename(domain_path)
                        protein1_domain, protein2_domain = domain_pair.split('+')
                        
                        # Use the provided function to find the necessary files
                        structure_file, json_file, log_file, PAE_png, fasta_file = find_rank_001_files(domain_path)
                        
                        # If the structure_file is missing, record the information
                        if not structure_file:
                            missing_info = {
                                'Protein1': protein1,
                                'Protein2': protein2,
                                'Protein1_Domain': protein1_domain,
                                'Protein2_Domain': protein2_domain
                            }
                            missing_files.append(missing_info)

                            # Optionally write to the CSV file
                            if save_to_csv:
                                writer.writerow([protein1, protein2, protein1_domain, protein2_domain])

                    except Exception as e:
                        print(f"Error checking domain pair {domain_pair}: {e}")
                
            except Exception as e:
                print(f"Error processing protein pair {protein1} and {protein2}: {e}")
    
    if save_to_csv:
        file.close()

    return missing_files

def find_missing_analyses(csv_path, base_folder):
    """
    Identify predictions that are present in the folder but missing from the CSV file.

    Parameters:
        - csv_path (str): Path to the CSV file containing analyzed predictions.
        - base_directory (str): Path to the base directory containing the folders for predictions.

    Returns:
        - missing_in_csv (set): A set of pairs present in folders but missing from the CSV file.
    """
    # Load the CSV file
    csv_data = pd.read_csv(csv_path)

    # Extract protein pairs and domain pairs from the CSV file
    csv_pairs = set(
        csv_data.apply(lambda row: f"{row[0]}_{row[1]}/{row[2]}+{row[3]}", axis=1)
    )

    # Walk through the base directory to find existing domain pairs
    folder_pairs = set()
    for root, dirs, files in os.walk(base_folder):
        for folder in dirs:
            # Check if the path format matches protein/domain pairs structure
            if '+' in folder and '_' in root:
                protein_pair = os.path.basename(root)
                domain_pair = folder
                full_pair = f"{protein_pair}/{domain_pair}"
                folder_pairs.add(full_pair)

    # Identify missing pairs
    missing_in_csv = folder_pairs - csv_pairs

    # Output results
    print("Pairs in folders but missing from CSV:")
    for pair in missing_in_csv:
        print(pair)
    
    return missing_in_csv

# Example usage
#base_folder = '/Users/poppy/Dropbox/2024.08.14_centriole_core_fragmented_proteins_fastas'
#base_folder = '/Users/poppy/Dropbox/2024.10.11_PCM_run_failures'
#csv_path = '/Users/poppy/Dropbox/PCM/PCM_alphafold_predictions_results_ROP25.csv'
#missing_files = find_missing_structure_files(base_folder, save_to_csv=True)
#print(f"Missing files: {missing_files}")
#missing_analyses = find_missing_analyses(csv_path, base_folder)
#print(f"Missing analyses: {missing_analyses}")
import os
import csv
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

# Example usage
#base_folder = '/Users/poppy/Dropbox/2024.08.14_centriole_core_fragmented_proteins_fastas'
#missing_files = find_missing_structure_files(base_folder, save_to_csv=True)
#print(f"Missing files: {missing_files}")
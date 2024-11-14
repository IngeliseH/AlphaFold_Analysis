"""
Script to run an analysis function and collect data in a csv on a large number of AlphaFold
predictions arranged in a specific folder structure. Expects a base folder containing protein pair
folders, each of which contains domain pair folders. For each domain pair folder, a specified
analysis function is called to process the data.

Functions:
    - process_all_predictions
"""
import os
import pandas as pd
from tqdm import tqdm  # Progress bar
import traceback
from iptm_visualisation import create_iptm_matrix, visualize_iptm_matrix

def process_all_predictions(
    base_folder,
    analysis_function,
    output_file="alphafold_predictions_results.csv",
    ipTM_graphic=True,
):
    """
    Process all AlphaFold predictions in the given base folder and write results to a CSV file.
    Expects base folder to contain protein pair folders, each of which contains domain pair folders.

    Parameters:
        - base_folder (str): Path to the base folder containing all protein pair folders.
        - output_file (str): Output CSV file name. Default is "alphafold_predictions_results.csv".
        - ipTM_graphic (bool): Whether to generate an ipTM matrix graphic. Default is True.
        - analysis_function (func): Function to process each domain pair folder. Should return a
          dictionary for each csv line to be written for a domain pair.
    """

    # Collect all data in a list of dictionaries
    all_data = []

    # Get list of protein pair folders
    protein_pair_folders = [d for d in os.listdir(base_folder) if os.path.isdir(os.path.join(base_folder, d))]

    total_pairs = len(protein_pair_folders)

    # Use tqdm for progress bar
    with tqdm(total=total_pairs, desc='Processing Protein Pairs') as pbar:
        for protein_pair_folder in protein_pair_folders:
            protein_pair_path = os.path.join(base_folder, protein_pair_folder)
            try:
                # Determine the separator to split the folder name
                if '+' in protein_pair_folder:
                    protein1, protein2 = protein_pair_folder.split('+')
                else:
                    protein1, protein2 = protein_pair_folder.split('_', 1)

                # Look for domain pair folders
                domain_folders = [
                    os.path.join(root, dir_name) for root, dirs, _ in os.walk(protein_pair_path)
                    for dir_name in dirs if '+' in dir_name
                ]

                if analysis_function:
                    for domain_path in domain_folders:
                        domain_pair = os.path.basename(domain_path)
                        protein1_domain, protein2_domain = domain_pair.split('+')

                        pbar.set_postfix_str(f"{protein1}:{protein1_domain} - {protein2}:{protein2_domain}")

                        try:
                            analysis_results = analysis_function(domain_path)
                            csv_entry = {
                                'Protein1': protein1,
                                'Protein2': protein2,
                                'Protein1_Domain': protein1_domain,
                                'Protein2_Domain': protein2_domain,
                            }
                            if isinstance(analysis_results, dict):
                                analysis_results = [analysis_results]

                            for item in analysis_results:
                                csv_entry_copy = csv_entry.copy()
                                csv_entry_copy.update({
                                    k: v for k, v in item.items()
                                    if k not in ['json_file', 'model_file', 'structure_model', 'residue_pairs', 'abs_res_lookup_dict', 'pae_data', 'confident_pairs']
                                })
                                all_data.append(csv_entry_copy)
                        except Exception as e:
                            print(f"Error processing domain pair {domain_pair}: {e}")

                # Optionally generate the ipTM matrix if is_pdb is True
                if ipTM_graphic:
                    try:
                        iptm_matrix = create_iptm_matrix(protein_pair_path)
                        visualize_iptm_matrix(iptm_matrix, os.path.join(protein_pair_path, 'iptm_matrix.png'))
                    except Exception as e:
                        print(f"Error generating ipTM matrix for {protein1} and {protein2}: {e}")
                        traceback.print_exc()

            except Exception as e:
                print(f"Error processing protein pair {protein_pair_folder}: {e}")
                traceback.print_exc()

            pbar.update(1)

    if all_data:
        headers = set().union(*[data.keys() for data in all_data])
        # order the headers to put protein and domain names first
        headers = ['Protein1', 'Protein2', 'Protein1_Domain', 'Protein2_Domain'] + sorted([h for h in headers if h not in ['Protein1', 'Protein2', 'Protein1_Domain', 'Protein2_Domain']])
        df = pd.DataFrame(all_data)
        df = df[list(headers)]
        df.to_csv(os.path.join(base_folder, output_file), index=False)
        print(f"Results written to {output_file}")
    else:
        print("No data to write.")

####################################################################################################
# Example usage:
#from interface_analysis import score_interaction
#base_folder = '/Users/poppy/Dropbox/PCM'
#process_all_predictions(base_folder, analysis_function=score_interaction, output_file="PCM_results.csv")

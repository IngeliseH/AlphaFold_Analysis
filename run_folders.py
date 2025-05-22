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
                if '_output' in protein_pair_folder: # ignore suffix '_output' if present
                    protein_pair_folder = protein_pair_folder.split('_output')[0]

                # Determine the separator to split the folder name
                if '+' in protein_pair_folder:
                    protein1, protein2 = protein_pair_folder.split('+')
                else:
                    parts = protein_pair_folder.split('_')
                    parts = [part for part in parts if part not in ['fold', 'fl']]
                    # Determine the pattern: p1_p2, p1_dimer_p2, p1_p2_dimer, or p1_dimer_p2_dimer
                    if len(parts) == 2:  # p1_p2
                        protein1 = parts[0]
                        protein2 = parts[1]
                    elif len(parts) == 3:
                        if parts[1] == 'dimer':  # p1_dimer_p2
                            protein1 = f"{parts[0]}_dimer"
                            protein2 = parts[2]
                        else:  # p1_p2_dimer
                            protein1 = parts[0]
                            protein2 = f"{parts[1]}_dimer"
                    elif len(parts) == 4:  # p1_dimer_p2_dimer
                        protein1 = f"{parts[0]}_dimer"
                        protein2 = f"{parts[2]}_dimer"
                    else:
                        raise ValueError(f"Invalid folder name format: {protein_pair_folder}")

                # Look for domain pair folders
                domain_folders = [
                    os.path.join(root, dir_name) for root, dirs, _ in os.walk(protein_pair_path)
                    for dir_name in dirs if '+' in dir_name
                ]

                if analysis_function:
                    for domain_path in domain_folders:
                        domain_pair = os.path.basename(domain_path)
                        # if present, remove '_output' from the domain pair name
                        if '_output' in domain_pair:
                            domain_pair = domain_pair.split('_output')[0]
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
                                    # removed residue_pairs and confident_pairs from this list to make further analysis easier
                                    if k not in ['json_file', 'model_file', 'structure_model', 'abs_res_lookup_dict', 'pae_data']
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

def process_FL_predictions(
    base_folder,
    analysis_function,
    output_file="alphafold_predictions_results.csv"
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
                if '_output' in protein_pair_folder: # ignore suffix '_output' if present
                    protein_pair_folder = protein_pair_folder.split('_output')[0]

                # Determine the separator to split the folder name
                if '+' in protein_pair_folder:
                    protein1, protein2 = protein_pair_folder.split('+')
                else:
                    parts = protein_pair_folder.split('_')
                    parts = [part for part in parts if part not in ['fold', 'fl']]
                    # Determine the pattern: p1_p2, p1_dimer, p1_dimer_p2, p1_p2_dimer, or p1_dimer_p2_dimer
                    if len(parts) == 2:  # p1_p2 or p1_dimer
                        protein1 = parts[0]
                        if parts[1] == 'dimer':
                            protein2 = parts[0]
                        else:
                            protein2 = parts[1]
                    elif len(parts) == 3:
                        if parts[1] == 'dimer':  # p1_dimer_p2
                            protein1 = f"{parts[0]}_dimer"
                            protein2 = parts[2]
                        else:  # p1_p2_dimer
                            protein1 = parts[0]
                            protein2 = f"{parts[1]}_dimer"
                    elif len(parts) == 4:  # p1_dimer_p2_dimer
                        protein1 = f"{parts[0]}_dimer"
                        protein2 = f"{parts[2]}_dimer"
                    else:
                        raise ValueError(f"Invalid folder name format: {protein_pair_folder}")

                if analysis_function:
                    pbar.set_postfix_str(f"{protein1} - {protein2}")
                    try:
                        analysis_results = analysis_function(protein_pair_path)
                        csv_entry = {
                            'Protein1': protein1,
                            'Protein2': protein2,
                        }
                        if isinstance(analysis_results, dict):
                            analysis_results = [analysis_results]

                        for item in analysis_results:
                            csv_entry_copy = csv_entry.copy()
                            csv_entry_copy.update({
                                k: v for k, v in item.items()
                                # leaving residue_pairs but leaving secondary_pairs as often exceed excel cell limit, causing issues
                                if k not in ['json_file', 'model_file', 'structure_model', 'abs_res_lookup_dict', 'pae_data', 'secondary_pairs']
                            })
                            all_data.append(csv_entry_copy)
                    except Exception as e:
                        print(f"Error processing protein pair {protein_pair_folder}: {e}")
            except Exception as e:
                print(f"Error processing protein pair {protein_pair_folder}: {e}")
                traceback.print_exc()

            pbar.update(1)

    if all_data:
        headers = set().union(*[data.keys() for data in all_data])
        # order the headers to put protein and domain names first
        headers = ['Protein1', 'Protein2'] + sorted([h for h in headers if h not in ['Protein1', 'Protein2']])
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

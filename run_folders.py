"""
Functions to process multiple folders of AlphaFold predictions

Functions:
process_alphafold_prediction
process_all_predictions
"""
import os
import csv
from analysis_utility import find_rank_001_files, parse_structure_file, extract_pae
from phosphorylation_handling import correct_cif_pae
from iptm_visualisation import extract_iptm
from process_pae import find_min_pae
from repeatability_from_pdb import measure_repeatability
from pdockq_calc import compute_pdockq

def process_alphafold_prediction(folder_path, is_pdb=True, **kwargs):
    """
    Process files from an AlphaFold prediction folder. Handles both AF2 (pdb) and AF3 (cif) outputs.

    Parameters:
    - folder_path (str): Path to the folder containing AlphaFold predictions.
    - is_pdb (bool): Whether the structure file is in PDB format. Default is True.
    - **kwargs: Additional function options can be specified as keyword arguments.

    Keyword Arguments (kwargs):
    - repeatability_params (dict): Parameters for measuring repeatability of predictions, including 'distance_cutoff'
      and 'pae_cutoff', and 'all_atom' to specify the atomic level of detail considered.
        Example: {'distance_cutoff': 5, 'pae_cutoff': 10, 'all_atom': True}

    Returns:
    - iptm (str): The iptm score of the top prediction
    - min_pae (float): The minimum interface PAE value
    - pdockq (float): The calculated pDockQ score
    - ppv (float): The calculated PPV score
    - num_consistent (int): The number of consistent predictions
    - level_consistent (str): The level of consistency (average score of >50% consistent predictions)
    """
    iptm, min_pae, pdockq, ppv, num_consistent, level_consistent = None, None, None, None, None, None
    structure_file, json_file, log_file, _ ,_ = find_rank_001_files(folder_path)
    # Parse structure file
    if structure_file:
        structure_model = parse_structure_file(structure_file, is_pdb)
    # Extract and potentially correct PAE matrix
    if json_file:
        pae_matrix = extract_pae(json_file)
        if not is_pdb and structure_model:
            pae_matrix = correct_cif_pae(structure_model, pae_matrix)

    # Find iptm score
    if log_file:
        iptm = extract_iptm(log_file)

    # Find minimum interface PAE
    if structure_model:
        min_pae, _ = find_min_pae(structure_model, pae_matrix)

    # Calculate pdockq and ppv
    if structure_file and json_file:
        pdockq, ppv = compute_pdockq(structure_file, json_file)

    # Measure repeatability with either provided or default parameters
    # get the pae_cutoff from the kwargs, if not provided, use the default value of 15
    pae_cutoff = kwargs.get('repeatability_params', {}).get('pae_cutoff', 15)
    if min_pae < pae_cutoff:
        num_consistent, level_consistent = measure_repeatability(folder_path, **kwargs.get('repeatability_params', {}))
    else:
        num_consistent = 0
        level_consistent = 'N/A'
    
    return iptm, min_pae, pdockq, ppv, num_consistent, level_consistent

def process_all_predictions(base_folder, output_path="alphafold_predictions_results.csv", **kwargs):
    """
    Process all AlphaFold predictions in the given base folder and write results to a CSV file.
    Expects base folder to contain protein pair folders, each of which contains domain pair folders.

    Parameters:
        - base_folder (str): Path to the base folder containing all protein pair folders.
        - output_path (str): Path to the output CSV file. Default is "alphafold_predictions_results.csv".
        - **kwargs: Additional function options can be specified as keyword arguments.

    Keyword Arguments (kwargs):
        - repeatability_params (dict): Parameters for measuring repeatability of predictions, including 'distance_cutoff'
          and 'pae_cutoff', and 'all_atom' to specify the atomic level of detail considered.
          Example: {'distance_cutoff': 5, 'pae_cutoff': 10, 'all_atom': True}

    Returns:
        - None, but writes the results to a CSV file.
    """
    headers = ['Protein1', 'Protein2', 'Protein1_Domain', 'Protein2_Domain', 'ipTM', 'min_PAE', 'pDockQ', 'PPV', 'Num_Consistent', 'Level_Consistent']

    with open(output_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(headers)  # Write the header row

        # Walk through the base folder containing all protein pair folders
        for root, dirs, _ in os.walk(base_folder):
            for dir in dirs:
                if '+' in dir:  # This is a protein pair folder
                    protein1, protein2 = dir.split('+')
                    print(f"Processing {protein1} and {protein2}...")
                    domain_folder_path = os.path.join(root, dir, 'Results')
                    
                    # Process each domain pair folder within the 'Results' folder
                    if os.path.exists(domain_folder_path):
                        for domain_pair in os.listdir(domain_folder_path):
                            if domain_pair.startswith(protein1) and '+' in domain_pair:
                                domain_path = os.path.join(domain_folder_path, domain_pair)
                                # Extract domain information
                                protein1_domain, protein2_domain = domain_pair.split('+')
                                print(f"Processing {protein1_domain} and {protein2_domain}...")
                                # Process the domain pair folder
                                iptm, min_pae, pdockq, ppv, num_consistent, level_consistent = process_alphafold_prediction(domain_path, **kwargs)

                                # Write results to CSV
                                writer.writerow([protein1, protein2, protein1_domain, protein2_domain, iptm, min_pae, pdockq, ppv, num_consistent, level_consistent])

####################################################################################################
# Example usage
#base_folder = "../../../../../Dropbox/2022.10.20_Drosophila_Version_1"
#process_all_predictions(base_folder)

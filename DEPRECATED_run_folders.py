"""
Functions to process multiple folders of AlphaFold predictions

Functions:
process_alphafold_prediction
process_all_predictions
"""
import os
import csv
import subprocess
from pathlib import Path
from analysis_utility import find_rank_001_files, parse_structure_file, extract_pae, map_chains_and_residues
from phosphorylation_handling import correct_cif_pae
from iptm_visualisation import extract_iptm, create_iptm_matrix, visualize_iptm_matrix
from process_pae import find_min_pae, compute_average_interface_pae
from repeatability_from_pdb import measure_repeatability, find_confident_pairs
from pdockq_calc import compute_pdockq
from process_MSA import process_alignment, calculate_interface_conservation_score, modify_bfactors

def process_alphafold_prediction(folder_path, is_pdb=True, **kwargs):
    """
    Process files from an AlphaFold prediction folder. Handles both AF2 (pdb) and AF3 (cif) outputs.

    Parameters:
        - folder_path (str): Path to the folder containing AlphaFold predictions.
        - is_pdb (bool): Whether the structure file is in PDB format. Default is True.
        - **kwargs: Additional function options can be specified as keyword arguments.

    Returns:
        - iptm (str): The iptm score of the top prediction
        - min_pae (float): The minimum interface PAE value
        - avg_int_pae (float): The average PAE value of the interface residue pairs with PAE < 15
        - pdockq (float): The calculated pDockQ score
        - ppv (float): The calculated PPV score
        - rop (int): The number of consistent predictions
        - percent_rop (str): The level of consistency (average score of >10% consistent predictions)
        - interface_size (int): The number of interface residue pairs found
        - conservation_score (float): The overall conservation score of the interface
    """
    iptm, min_pae, avg_int_pae, pdockq, ppv, rop, percent_rop, interface_size, conservation_score = None, None, None, None, None, None, None, None, None
    residue_pairs = None
    structure_file, structure_model, json_file, log_file = None, None, None, None
    structure_file, json_file, log_file, _, _ = find_rank_001_files(folder_path)
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
        
        chain_residue_map = map_chains_and_residues(structure_model)
        # Create mappings between absolute residue indices and (chain_id, res_num)
        abs_res_lookup_dict = {(entry[0], entry[1]): entry[3] for entry in chain_residue_map}
        abs_res_reverse_lookup = {entry[3]: (entry[0], entry[1]) for entry in chain_residue_map}
        residue_pairs = find_confident_pairs(
            structure_model, json_file, distance_cutoff=5, pae_cutoff=15,
            abs_res_lookup_dict=abs_res_lookup_dict, is_pdb=is_pdb, all_atom=True
        )
        interface_size = len(residue_pairs)

        # Find average interface PAE if confident_pairs are found
        if residue_pairs:
            avg_int_pae = compute_average_interface_pae(pae_matrix, residue_pairs)

    # Measure repeatability
    pae_cutoff = kwargs.get('repeatability_params', {}).get('pae_cutoff', 15)
    if min_pae and min_pae < pae_cutoff:
        rop, percent_rop = measure_repeatability(folder_path, **kwargs.get('repeatability_params', {}))
    else:
        rop = 0
        percent_rop = 'N/A'

    # Process MSA and calculate conservation score
    # Convert .a3m to FASTA and calculate conservation score if confident interface residues are found
    if residue_pairs and rop >= 5: # blocking this from running as was failing for some and meaing they were left out of the overall file
        # Find the .a3m file in the folder
        msa_files = [f for f in os.listdir(folder_path) if f.endswith('.a3m')]
        if msa_files:
            a3m_file = os.path.join(folder_path, msa_files[0])  # Take the first .a3m file found
            fasta_file = os.path.join(folder_path, msa_files[0].replace('.a3m', '.fasta'))

            # Use reformat.pl to convert .a3m to FASTA
            reformat_command = f"perl reformat.pl {a3m_file} {fasta_file}"
            subprocess.run(reformat_command, shell=True, check=True)

            # Check if the FASTA file was created successfully
            if os.path.exists(fasta_file):
                # Get residue conservation scores
                residue_similarity_scores_list, start_positions = process_alignment(fasta_file, structure_file)

                # Compute conservation score using the new function
                conservation_score = calculate_interface_conservation_score(
                    residue_similarity_scores_list,
                    residue_pairs,
                    abs_res_reverse_lookup,
                    start_positions
                )

                # Optionally save the PDB file
                if rop >= 3:
                    protein_model = parse_structure_file(structure_file)
                    if isinstance(structure_file, str):
                        output_name = structure_file.split(".")[0] + "_conservation.pdb"
                    elif isinstance(structure_file, Path):
                        output_name = structure_file.parent / (structure_file.stem + "_conservation.pdb")
                    modify_bfactors(protein_model, residue_similarity_scores_list, output_name)

            else:
                print(f"FASTA conversion failed for {a3m_file}")
        else:
            print(f"No .a3m file found in {os.path.join(folder_path, 'msas')}")
    else:
        conservation_score = None

    return iptm, min_pae, avg_int_pae, pdockq, ppv, rop, percent_rop, interface_size, conservation_score

def process_all_predictions(base_folder, output_file="alphafold_predictions_results.csv", is_pdb=True, ipTM_graphic=True, **kwargs):
    """
    Process all AlphaFold predictions in the given base folder and write results to a CSV file.
    Expects base folder to contain protein pair folders, each of which contains domain pair folders.

    Parameters:
        - base_folder (str): Path to the base folder containing all protein pair folders.
        - output_file (str): Output CSV file name. Default is "alphafold_predictions_results.csv".
        - is_pdb (bool): Whether the structure file is in PDB format. Default is True.
        - ipTM_graphic (bool): Whether to generate an ipTM matrix graphic. Default is True.
        - **kwargs: Additional function options can be specified as keyword arguments.

    Keyword Arguments (kwargs):
        - repeatability_params (dict): Parameters for measuring repeatability of predictions, including 'distance_cutoff'
          and 'pae_cutoff', and 'all_atom' to specify the atomic level of detail considered.
          Example: {'distance_cutoff': 5, 'pae_cutoff': 10, 'all_atom': True}

    Returns:
        - None, but writes the results to a CSV file.
    """
    headers = ['Protein1', 'Protein2', 'Protein1_Domain', 'Protein2_Domain', 'ipTM', 'min_PAE', 'average_interface_PAE', 'pDockQ', 'ppv', 'ROP', 'percent_ROP', 'interface_size', 'BLOSUM_conservation_score']
    output_path = os.path.join(base_folder, output_file)
    with open(output_path, mode='w', newline='') as file:
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
                    print(f"Processing {protein1} and {protein2}...")

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

                    # Process each domain pair folder identified
                    for domain_path in domain_folders:
                        try:
                            domain_pair = os.path.basename(domain_path)
                            protein1_domain, protein2_domain = domain_pair.split('+')
                            # Process the domain pair folder
                            iptm, min_pae, avg_int_pae, pdockq, ppv, rop, percent_rop, interface_size, conservation_score = process_alphafold_prediction(domain_path, is_pdb, **kwargs)

                            # Write results to CSV
                            writer.writerow([protein1, protein2, protein1_domain, protein2_domain, iptm, min_pae, avg_int_pae, pdockq, ppv, rop, percent_rop, interface_size, conservation_score])
                        except Exception as e:
                            print(f"Error processing domain pair {domain_pair}: {e}")

                    # Optionally generate the ipTM matrix if is_pdb is True
                    if ipTM_graphic and is_pdb:
                        try:
                            iptm_matrix = create_iptm_matrix(protein_pair_path)
                            png_file_path = os.path.join(protein_pair_path, 'iptm_matrix.png')
                            visualize_iptm_matrix(iptm_matrix, png_file_path)
                        except Exception as e:
                            print(f"Error generating ipTM matrix for {protein1} and {protein2}: {e}")
                
                except Exception as e:
                    print(f"Error processing protein pair {protein1} and {protein2}: {e}")

####################################################################################################
# Example usage
#base_folder = "../../../../../Dropbox/2022.10.20_Drosophila_Version_1"
base_folder = '/Users/poppy/Dropbox/to_analyse'
process_all_predictions(base_folder, output_file="alphafold_predictions_results.csv", is_pdb=True, ipTM_graphic=False)

# For testing completion bias
# path = "data/CENPJ_PALB2/CENPJ_D5+PALB2_D5"
# print(process_alphafold_prediction(path))

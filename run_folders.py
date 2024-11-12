"""
Script to process AlphaFold predictions for protein-protein interactions.
For each domain pair prediction:
    - Finds files, pae data, interface residue pairs (and confident subset), structure model object, and
      abs_res_lookup_dict for each model
    - Selects the best model based on ROP scores and min_pae increase criteria (the best model is the one with
      the highest ROP score, unless it has significantly higher min_pae than the original top-ranked model)
    - Computes additional metrics for the best model (avg_interface_pae, pdockq, iptm, interface size,
      interface residues with chain-specific numbering, pae evenness)

TODO: - iptm matrix generation feature currently not working
      - Check compatibility with AF3 or multichain predictions

Functions:
    - collect_model_data
    - select_best_model
    - compute_additional_metrics
    - score_interaction
    - process_all_predictions
"""

import re
import os
import pandas as pd
import numpy as np
from tqdm import tqdm  # Progress bar
import traceback
from pathlib import Path
from analysis_utility import extract_pae, parse_structure_file, map_chains_and_residues
from repeatability_from_pdb import get_residue_pairs, find_confident_pairs, calculate_rop_scores, calculate_percent_rop
from process_pae import model_dictionary_min_pae, compute_average_interface_pae, compute_pae_evenness
from iptm_visualisation import extract_iptm, create_iptm_matrix, visualize_iptm_matrix
from pdockq_calc import compute_pdockq

def collect_model_data(folder_path, distance_cutoff=10.0, pae_cutoff=15.0, all_atom=True):
    """
    Collects essential model data (all_pairs and confident_pairs) required for ROP calculations.
    Returns a list of model_data dictionaries and a boolean is_pdb indicating file type.
    """
    # Determine file extension and set is_pdb flag
    is_pdb = any(Path(folder_path).glob('*.pdb'))
    file_extension = '*.pdb' if is_pdb else '*.cif'

    # Collect structure files
    structure_files = sorted(
        [x for x in Path(folder_path).glob(file_extension) if 'conservation' not in x.name],
        key=lambda x: int(re.search(r'rank_(\d+)', x.name).group(1))
    )

    model_data = []
    for model_file in structure_files:
        # Extract model rank from filename
        rank_match = re.search(r'rank_(\d+)', model_file.name)
        if rank_match:
            model_rank = rank_match.group(1)
        else:
            model_rank = None  # Handle as needed

        # Parse structure model
        structure_model = parse_structure_file(model_file, is_pdb)

        # Adjust JSON file name if it includes 'unrelaxed_' prefix
        json_file = model_file.with_name(model_file.name.replace("unrelaxed", "scores")).with_suffix('.json')
        # Extract PAE data
        pae_data = extract_pae(json_file)

        # Map chains and residues
        chain_residue_map = map_chains_and_residues(structure_model)
        abs_res_lookup_dict = {(chain_id, res_id): abs_res_id for chain_id, res_id, _, abs_res_id in chain_residue_map}

        # Get interprotein residue pairs within distance cutoff (all_pairs)
        all_pairs = set(get_residue_pairs(structure_model, distance_cutoff, abs_res_lookup_dict, all_atom))
        confident_pairs = find_confident_pairs(pae_data, pae_cutoff, all_pairs)

        # Store model data
        model_data.append({
            'model_file': model_file,
            'model_rank': model_rank,
            'json_file': json_file,
            'structure_model': structure_model,
            'pae_data': pae_data,
            'all_pairs': all_pairs,
            'confident_pairs': confident_pairs,
            'rop_score': 0, # Placeholder for ROP score
            'abs_res_lookup_dict': abs_res_lookup_dict
            # Other metrics will be computed and added later
        })

    return model_data

def select_best_model(model_data):
    """
    Selects the best model based on ROP scores and min_pae increase criteria.
    Returns the best model data dictionary.
    """
    # Identify the top-ranked model (the first in the list)
    top_model = model_data[0]
    top_model_rank = top_model['model_rank']

    models_with_confident_pairs = [m for m in model_data if m['confident_pairs']]

    if not models_with_confident_pairs:
        # If none of the models have any confident pairs, keep the initial top-ranked model
        #print("No models have confident pairs; keeping the initial top-ranked model.")
        top_model['min_pae'] = model_dictionary_min_pae(top_model)
        return top_model

    # Find model(s) with highest ROP score
    max_rop = max(model['rop_score'] for model in models_with_confident_pairs)
    models_with_max_rop = [model for model in models_with_confident_pairs if model['rop_score'] == max_rop]

    # If multiple models have equal highest ROP, pick the original highest-ranked one
    best_model = min(models_with_max_rop, key=lambda m: abs(int(m['model_rank']) - int(top_model_rank)))

    # Now compute min_pae for top_model and best_model
    top_model_min_pae = model_dictionary_min_pae(top_model)
    best_model_min_pae = model_dictionary_min_pae(best_model)

    # Check if picking this model leads to a significant rise in minPAE
    min_pae_increase = best_model_min_pae - top_model_min_pae

    # Compute the threshold for significant rise in minPAE
    min_pae_threshold = max(3, int(np.ceil(0.3 * top_model_min_pae)))

    if min_pae_increase > min_pae_threshold:
        # Significant rise in minPAE; keep the original top-ranked model
        #print(f"Keeping rank1 - top model minPAE = {top_model_min_pae}, best model minPAE = {best_model_min_pae}")
        best_model = top_model
    #else:
        #if best_model != top_model:
            #print(f"Model {best_model['model_rank']} selected, ROP {best_model['rop_score']}, vs rank1 {top_model['rop_score']}.")

    # Update min_pae in best_model data
    best_model['min_pae'] = best_model_min_pae

    return best_model

def compute_additional_metrics(model):
    """
    Computes avg_interface_pae, pdockq, and iptm for the given model.
    Updates the model data dictionary with these values.
    """
    confident_pairs = model['confident_pairs']
    pae_data = model['pae_data']
    model_file = model['model_file']
    json_file = model['json_file']

    # Compute min_pae if not already found (ie if not using best_model selection, as this finds pae)
    if 'min_pae' not in model:
        model['min_pae'] = model_dictionary_min_pae(model)
    # Compute avg_interface_pae
    model['avg_interface_pae'] = compute_average_interface_pae(pae_data, confident_pairs)
    # Compute pae_evenness
    model['pae_evenness'] = compute_pae_evenness(pae_data, confident_pairs)
    # Compute pdockq
    model['pdockq'], _ = compute_pdockq(model_file, json_file)
    # Extract iptm
    model['iptm'] = extract_iptm(json_file, model['model_rank'])
    # Find number of confident interface residue pairs
    model['interface_size'] = len(confident_pairs)
    
    # Find interface regions in terms of chain-specific residue IDs
    chain_ids = [chain.id for chain in model['structure_model'].get_chains()]
    # Map absolute residue IDs back to chain-specific numbering with chain IDs
    inv_abs_res_lookup_dict = {v: k for k, v in model['abs_res_lookup_dict'].items()}
    confident_pairs_chain_specific = [
        (f"{inv_abs_res_lookup_dict[res1][0]} {inv_abs_res_lookup_dict[res1][1]}",
         f"{inv_abs_res_lookup_dict[res2][0]} {inv_abs_res_lookup_dict[res2][1]}")
        for res1, res2 in confident_pairs if res1 in inv_abs_res_lookup_dict and res2 in inv_abs_res_lookup_dict
    ]
    model['confident_pairs_chain_specific'] = confident_pairs_chain_specific
    # Group residues by protein and format continuous ranges
    residues_by_protein = {'Protein1': set(), 'Protein2': set()}
    for res1, res2 in confident_pairs_chain_specific:
        chain1, res_num1 = res1.split()
        chain2, res_num2 = res2.split()
        if chain1 == chain_ids[0]:
            residues_by_protein['Protein1'].add(int(res_num1))
            residues_by_protein['Protein2'].add(int(res_num2))
        else:
            residues_by_protein['Protein1'].add(int(res_num2))
            residues_by_protein['Protein2'].add(int(res_num1))
    def format_ranges(numbers):
        sorted_nums = sorted(numbers)
        ranges = []
        start = end = sorted_nums[0]
        for n in sorted_nums[1:]:
            if n == end + 1:
                end = n
            else:
                ranges.append(f"{start}-{end}" if start != end else str(start))
                start = end = n
        ranges.append(f"{start}-{end}" if start != end else str(start))
        return ", ".join(ranges)

    model['confident_residues_by_protein'] = {
        'Protein1': format_ranges(residues_by_protein['Protein1']),
        'Protein2': format_ranges(residues_by_protein['Protein2'])
    }

def score_interaction(folder_path, distance_cutoff=5.0, pae_cutoff=15.0, all_atom=True):
    """
    Main function that selects and scores the best model.
    """
    # Collect model data (only essential data)
    model_data = collect_model_data(folder_path, distance_cutoff, pae_cutoff, all_atom)

    # Calculate ROP scores
    calculate_rop_scores(model_data)

    # Select the best model
    best_model = select_best_model(model_data)

    # Calculate percent_rop
    best_model['percent_rop'] = calculate_percent_rop(model_data, best_model)

    # Compute additional metrics for the best model
    compute_additional_metrics(best_model)

    return best_model

def process_all_predictions(base_folder, output_file="alphafold_predictions_results.csv", ipTM_graphic=True):
    """
    Process all AlphaFold predictions in the given base folder and write results to a CSV file.
    Expects base folder to contain protein pair folders, each of which contains domain pair folders.

    Parameters:
        - base_folder (str): Path to the base folder containing all protein pair folders.
        - output_file (str): Output CSV file name. Default is "alphafold_predictions_results.csv".
        - ipTM_graphic (bool): Whether to generate an ipTM matrix graphic. Default is True.

    Returns:
        - None, but writes the results to a CSV file.
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
                domain_folders = []
                for root, dirs, _ in os.walk(protein_pair_path):
                    for dir_name in dirs:
                        if '+' in dir_name:
                            domain_folders.append(os.path.join(root, dir_name))
                    break  # Only process the immediate subdirectories

                # Process each domain pair folder
                for domain_path in domain_folders:
                    domain_pair = os.path.basename(domain_path)
                    protein1_domain, protein2_domain = domain_pair.split('+')

                    current_task = f"{protein1}:{protein1_domain} - {protein2}:{protein2_domain}"
                    pbar.set_postfix_str(current_task)

                    try:
                        model = score_interaction(domain_path)

                        # model is expected to be a dictionary of attributes
                        # Add additional information
                        model_data = {
                            'Protein1': protein1,
                            'Protein2': protein2,
                            'Protein1_Domain': protein1_domain,
                            'Protein2_Domain': protein2_domain,
                        }

                        # Assuming model is a dictionary of attributes
                        if isinstance(model, dict):
                            # if present, exclude unnecessary keys
                            model_data.update({k: v for k, v in model.items() if k not in ['json_file', 'model_file', 'structure_model', 'confident_pairs', 'all_pairs', 'pae_data', 'abs_res_lookup_dict']})
                        else:
                            print(f"Model returned is not a dictionary for {current_task}. Skipping.")
                            continue

                        # Append to all_data
                        all_data.append(model_data)

                    except Exception as e:
                        print(f"Error processing domain pair {domain_pair}: {e}")
                        traceback.print_exc()

                # Optionally generate the ipTM matrix if is_pdb is True
                if ipTM_graphic:
                    try:
                        iptm_matrix = create_iptm_matrix(protein_pair_path)
                        png_file_path = os.path.join(protein_pair_path, 'iptm_matrix.png')
                        visualize_iptm_matrix(iptm_matrix, png_file_path)
                    except Exception as e:
                        print(f"Error generating ipTM matrix for {protein1} and {protein2}: {e}")
                        traceback.print_exc()

            except Exception as e:
                print(f"Error processing protein pair {protein_pair_folder}: {e}")
                traceback.print_exc()

            pbar.update(1)

    # After processing all data, write to CSV
    if all_data:
        # Dynamically determine the headers based on keys in all_data
        headers = set()
        for data in all_data:
            headers.update(data.keys())
        headers = sorted(headers)  # Sort headers for consistency

        # Write to CSV using pandas for simplicity
        df = pd.DataFrame(all_data)
        df = df[headers]  # Ensure columns are in the order of headers
        output_path = os.path.join(base_folder, output_file)
        df.to_csv(output_path, index=False)
        print(f"Results written to {output_path}")
    else:
        print("No data to write.")

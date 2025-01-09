"""
Script to process AlphaFold predictions as dictionaries, selecting the best model based on ROP scores and min_pae
For each prediction folder:
    - Finds files, pae data, interface residue pairs (and confident subset), structure model object, and
      abs_res_lookup_dict for each model
    - Selects the best model based on ROP scores and min_pae increase criteria (the best model is the one with
      the highest ROP score, unless it has significantly higher min_pae than the original top-ranked model)
    - Computes additional metrics for the best model (avg_interface_pae, pdockq, iptm, interface size,
      interface residues with chain-specific numbering, pae evenness)

TODO: - Check compatibility with AF3 or multichain predictions

Functions:
    - collect_model_data
    - select_best_model
    - compute_additional_metrics
    - score_interaction
"""

import re
import numpy as np
from pathlib import Path
from analysis_utility import extract_pae, parse_structure_file, map_chains_and_residues
from repeatability_from_pdb import get_residue_pairs, find_confident_pairs, calculate_rop_score, calculate_percent_rop
from process_pae import residue_pairs_min_pae, compute_average_interface_pae, compute_pae_evenness
from iptm_visualisation import extract_iptm
from pdockq_calc import compute_pdockq

def collect_model_data(folder_path, distance_cutoff=5.0, pae_cutoff=15.0, all_atom=True):
    """
    Collects essential model data (residue_pairs and confident_pairs) required for ROP calculations.
    Returns a list of model_data dictionaries and a boolean is_pdb indicating file type.

    Parameters:
        - folder_path (str): Path to the folder containing AlphaFold structure files
        - distance_cutoff (float): Distance cutoff for identifying interprotein residue pairs - 
          recommended to use 5 if all_atom is True, 10 if all_atom is False
        - pae_cutoff (float): PAE cutoff for identifying confident residue pairs

    """
    # Determine file extension and set is_pdb flag
    is_pdb = any(Path(folder_path).glob('*.pdb'))
    file_extension = '*.pdb' if is_pdb else '*.cif'

    # Collect structure files
    structure_files = sorted(
        [x for x in Path(folder_path).glob(file_extension) if 'conservation' not in x.name],
        key=lambda x: int(re.search(r'(rank|model)_(\d+)', x.name).group(2))
    )

    model_data = []
    for model_file in structure_files:
        # Extract model rank from filename
        rank_match = re.search(r'(rank|model)_(\d+)', model_file.name)
        if rank_match:
            model_rank = rank_match.group(2)
        else:
            model_rank = None  # Handle as needed

        # Parse structure model
        structure_model = parse_structure_file(model_file, is_pdb)

        if is_pdb:
            # Adjust JSON file name if it includes 'unrelaxed_' prefix
            json_file = model_file.with_name(model_file.name.replace("unrelaxed", "scores")).with_suffix('.json')
            # iptm data is in json file
            log_file = json_file
        else:
            # For CIF files, JSON file name is same as model file name with "model" replaced with "full_data"
            json_file = model_file.with_name(model_file.name.replace("model", "full_data")).with_suffix('.json')
            # log fie with iptm is named with summary_confidences instead of full_data
            log_file = model_file.with_name(model_file.name.replace("model", "summary_confidences")).with_suffix('.json')
        # Extract PAE data
        pae_data = extract_pae(json_file)

        # Map chains and residues
        chain_residue_map = map_chains_and_residues(structure_model)
        abs_res_lookup_dict = {(chain_id, res_id): abs_res_id for chain_id, res_id, _, abs_res_id in chain_residue_map}

        # Get interprotein residue pairs within distance cutoff (residue_pairs)
        residue_pairs = set(get_residue_pairs(structure_model, distance_cutoff, abs_res_lookup_dict, all_atom))
        confident_pairs = find_confident_pairs(residue_pairs, pae_data, pae_cutoff)
        secondary_pairs = get_residue_pairs(structure_model, distance_cutoff+1, abs_res_lookup_dict, all_atom)

        # Store model data
        model_data.append({
            'model_file': model_file,
            'model_rank': model_rank,
            'json_file': json_file,
            'log_file': log_file,
            'structure_model': structure_model,
            'pae_data': pae_data,
            'residue_pairs': residue_pairs,
            'confident_pairs': confident_pairs,
            'secondary_pairs': secondary_pairs,
            'rop': 0, # Placeholder for ROP score
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
    top_model['min_pae'] = residue_pairs_min_pae(top_model['residue_pairs'], top_model['pae_data'])

    models_with_confident_pairs = [m for m in model_data if m['confident_pairs']]

    if not models_with_confident_pairs:
        # If none of the models have any confident pairs, keep the initial top-ranked model
        return top_model

    # Find model(s) with highest ROP score
    max_rop = max(model['rop'] for model in models_with_confident_pairs)
    models_with_max_rop = [model for model in models_with_confident_pairs if model['rop'] == max_rop]

    # If multiple models have equal highest ROP, pick the original highest-ranked one
    best_model = min(models_with_max_rop, key=lambda m: abs(int(m['model_rank']) - int(top_model_rank)))

    # Now compute min_pae for best_model
    best_model['min_pae'] = residue_pairs_min_pae(best_model['confident_pairs'], best_model['pae_data'])

    # Check if picking this model leads to a significant rise in minPAE
    min_pae_increase = best_model['min_pae'] - top_model['min_pae']

    # Compute the threshold for significant rise in minPAE
    min_pae_threshold = max(3, int(np.ceil(0.3 * top_model['min_pae'])))

    if min_pae_increase > min_pae_threshold:
        # Significant rise in minPAE; keep the original top-ranked model
        best_model = top_model

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
    log_file = model['log_file']

    # Compute min_pae if not already found (ie if not using best_model selection, as this finds pae)
    if 'min_pae' not in model:
        model['min_pae'] = residue_pairs_min_pae(model['confident_pairs'], model['pae_data'])
    # Compute avg_interface_pae
    model['avg_interface_pae'] = compute_average_interface_pae(confident_pairs, pae_data)
    # Compute pae_evenness
    model['pae_evenness'] = compute_pae_evenness(confident_pairs, pae_data)
    # Compute pdockq if only 2 chains
    if len(model['structure_model']) == 2:
        model['pdockq'], _ = compute_pdockq(model_file, json_file)
    # Extract iptm
    model['iptm'] = extract_iptm(log_file, model['model_rank'])
    # Find number of confident interface residue pairs
    model['interface_size'] = len(confident_pairs)
    
    # Map absolute residue IDs back to chain-specific numbering with chain IDs
    inv_abs_res_lookup_dict = {v: k for k, v in model['abs_res_lookup_dict'].items()}
    confident_pairs_chain_specific = [
        (f"{inv_abs_res_lookup_dict[res1][0]} {inv_abs_res_lookup_dict[res1][1]}",
         f"{inv_abs_res_lookup_dict[res2][0]} {inv_abs_res_lookup_dict[res2][1]}")
        for res1, res2 in confident_pairs if res1 in inv_abs_res_lookup_dict and res2 in inv_abs_res_lookup_dict
    ]
    model['confident_pairs_chain_specific'] = confident_pairs_chain_specific

def score_interaction(folder_path, distance_cutoff=5.0, pae_cutoff=15.0, all_atom=True):
    """
    Main function that selects and scores the best model.
    """
    # Collect model data (only essential data)
    model_data = collect_model_data(folder_path, distance_cutoff, pae_cutoff, all_atom)

    # Calculate ROP scores
    for model in model_data:
        other_model_pairs = [m['confident_pairs'] for m in model_data if m != model]
        model['rop'] = calculate_rop_score(model['confident_pairs'], other_model_pairs)

    # Select the best model
    best_model = select_best_model(model_data)

    # Calculate percent_rop
    for model in model_data:
        other_model_pairs = [m['confident_pairs'] for m in model_data if m != model]
    best_model['percent_rop'] = calculate_percent_rop(best_model['confident_pairs'], other_model_pairs)

    # Compute additional metrics for the best model
    compute_additional_metrics(best_model)

    return best_model

####################################################################################################
# # Example usage:
# folder_path = 'data/Ana2_Sak/Ana2_D1+Sak_D1'
# distance_cutoff = 5.0
# pae_cutoff = 15.0
# all_atom = True

# # USING INDIVIDUAL FUNCTIONS
# # Collect model data
# model_data = collect_model_data(folder_path, distance_cutoff, pae_cutoff, all_atom)
# print(f"Collected data for {len(model_data)} models.")
# # Calculate ROP scores for each model
# for model in model_data:
#     other_model_pairs = [m['confident_pairs'] for m in model_data if m != model]
#     model['rop'] = calculate_rop_score(model['confident_pairs'], other_model_pairs)
# print("Calculated ROP scores for all models.")
# # Select the best model based on ROP scores and min_pae increase criteria
# best_model = select_best_model(model_data)
# print(f"Selected best model: {best_model['model_file']} with ROP score: {best_model['rop']}")
# # Compute additional metrics for the best model
# compute_additional_metrics(best_model)
# print("Computed additional metrics for the best model.")
# # Print out the best model's data
# print("Best model data:")
# for key, value in best_model.items():
#     print(f"{key}: {value}")

# # USING OVERALL SCORING FUNCTION
# # Score interaction for the best model
# best_model = score_interaction(folder_path, distance_cutoff, pae_cutoff, all_atom)
# print(f"Best model selected: {best_model['model_file']}")
# print(f"ROP score: {best_model['rop']}")
# print(f"Percent ROP: {best_model['percent_rop']}")
# print(f"Min PAE: {best_model['min_pae']}")
# print(f"Avg Interface PAE: {best_model['avg_interface_pae']}")
# print(f"PAE Evenness: {best_model['pae_evenness']}")
# print(f"pDockQ: {best_model['pdockq']}")
# print(f"ipTM: {best_model['iptm']}")
# print(f"Interface Size: {best_model['interface_size']}")
# print(f"Confident Pairs (Chain Specific): {best_model['confident_pairs_chain_specific']}")

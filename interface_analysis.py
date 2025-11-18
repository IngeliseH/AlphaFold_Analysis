"""
for each domain pair:
    - make model object for each model
    - find residue pairs, confident residue pairs and rop for each
    - find best model (based on rop and pae)
    - find pdockq, iptm for best model
    - find all interfaces in best model and size of each
    - for each interface > 1, avg_pae and confidence class (based on this avg_pae)
    - for each high or medium confidence interface, find min_pae and rop, and residues involved in readable format
        - find rop by set of residue pairs in interface to sets of confident residue pairs in each other model
        - set higher threshold - 70%
    - return list of dictionaries of data for each high or medium confidence interface
"""
import networkx as nx
import numpy as np
from scipy.spatial import cKDTree
from iptm_visualisation import extract_iptm
from process_pae import compute_average_interface_pae, compute_pae_evenness, residue_pairs_min_pae
from pdockq_calc import compute_pdockq
from analyse_models import collect_model_data, select_best_model
from repeatability_from_pdb import calculate_rop_score, calculate_rop_score, calculate_percent_rops
from analysis_utility import readable_ranges
#from interface_character import charge_score, hydrophobicity_score

def get_residue_coordinate(abs_res_id, inv_abs_res_lookup_dict, structure_model, use_side_chain=True):
    """
    Get the representative coordinate of a residue.

    Parameters:
        - abs_res_id (int): The absolute residue ID.
        - inv_abs_res_lookup_dict (dict): Mapping from abs_res_id to (chain_id, res_num).
        - structure_model (Bio.PDB.Model.Model): The protein model object.
        - use_side_chain (bool): If True, use centroid of side-chain atoms; if False, use the CA atom.

    Returns:
        - coord (numpy.ndarray): The representative coordinate of the residue.
    """
    # Retrieve chain ID and residue number
    chain_id, res_num = inv_abs_res_lookup_dict.get(abs_res_id, (None, None))
    if chain_id is None or res_num is None:
        raise ValueError(f"Residue with abs_res_id {abs_res_id} not found in inv_abs_res_lookup_dict.")

    # Access the chain
    chain = structure_model[chain_id]

    # Construct the residue ID tuple (hetfield, resseq, icode)
    residue_id = (' ', res_num, ' ')  # Adjust as necessary for your data

    # Retrieve the residue object
    residue = chain[residue_id]
    if residue is None:
        raise ValueError(f"Residue {residue_id} not found in chain {chain_id}.")

    # Compute the representative coordinate
    if use_side_chain:
        side_chain_atoms = [atom for atom in residue.get_atoms()
                            if atom.get_name() not in ('N', 'CA', 'C', 'O')]
        if side_chain_atoms:
            coords = np.array([atom.get_coord() for atom in side_chain_atoms])
            centroid = np.mean(coords, axis=0)
            return centroid
        elif 'CA' in residue:
            return residue['CA'].get_coord()
        else:
            raise ValueError(f"No side-chain or CA atoms found for residue {abs_res_id}.")
    else:
        if 'CA' in residue:
            return residue['CA'].get_coord()
        else:
            raise ValueError(f"CA atom not found for residue {abs_res_id}.")

def build_residue_pair_graph_with_spatial_proximity_and_pae(
    residue_pairs,
    inv_abs_res_lookup_dict,
    structure_model,
    pae_data,
    pae_cutoff=15.0,
    distance_cutoff=5.0
):
    """
    Build a graph where nodes are residue pairs, and edges connect pairs
    that are spatially close and have low within-chain PAE.

    Parameters:
        - residue_pairs (list of tuples): List of tuples of absolute residue IDs.
        - inv_abs_res_lookup_dict (dict): Mapping from abs_res_id to (chain_id, res_num).
        - structure_model (Bio.PDB.Model.Model): The protein model object.
        - pae_data (dict): PAE data with keys as (res1_id, res2_id) tuples.
        - pae_cutoff (float): Maximum PAE between residues within the same chain.
        - distance_cutoff (float): Maximum spatial distance between residue pairs
          to be considered part of the same interface.

    Returns:
        - G (networkx.Graph): The graph with residue pairs as nodes.
    """
    G = nx.Graph()

    # Add node for each interprotein residue pair
    # Node is associated with midpoint coordinate of the pair
    rep_coords = []
    for idx, (res1_id, res2_id) in enumerate(residue_pairs):
        coord1 = get_residue_coordinate(res1_id, inv_abs_res_lookup_dict, structure_model)
        coord2 = get_residue_coordinate(res2_id, inv_abs_res_lookup_dict, structure_model)
        rep_coord = (coord1 + coord2) / 2
        G.add_node(idx, pair=(res1_id, res2_id), coord=rep_coord)
        rep_coords.append(rep_coord)

    # Build edges based on spatial proximity and PAE
    rep_coords_array = np.array(rep_coords)
    tree = cKDTree(rep_coords_array)
    pairs = tree.query_pairs(r=distance_cutoff)
    for i, j in pairs:
        pair_i = G.nodes[i]['pair']
        pair_j = G.nodes[j]['pair']

        # Check within-chain PAE for residues in the same chain
        # Exclude input residue pairs (from different proteins) from PAE filtering
        residues_i = pair_i
        residues_j = pair_j

        # Collect all unique residues involved
        all_residues = set(residues_i + residues_j)

        # Group residues by chain
        chain_residues = {}
        for res_id in all_residues:
            chain_id, res_num = inv_abs_res_lookup_dict[res_id]
            chain_residues.setdefault(chain_id, []).append(res_id)

        # Check PAE between residues within the same chain
        low_pae_within = True
        for chain_id, res_ids in chain_residues.items():
            # If only one residue in the chain, skip PAE check
            if len(res_ids) <= 1:
                continue
            # Compute pairwise PAE between residues within the chain
            for res_id1 in res_ids:
                for res_id2 in res_ids:
                    if res_id1 >= res_id2:
                        continue  # Avoid duplicates
                    # Exclude input residue pairs (from different proteins)
                    if (res_id1, res_id2) in residue_pairs or (res_id2, res_id1) in residue_pairs:
                        continue
                    pae = min(pae_data[res_id1][res_id2], pae_data[res_id2][res_id1])
                    if pae > pae_cutoff:
                        low_pae_within = False
                        break
                if not low_pae_within:
                    break
            if not low_pae_within:
                break

        if low_pae_within:
            G.add_edge(i, j)

    return G

def get_interfaces_from_graph(G):
    """
    Identify interfaces as connected components in the residue pair graph.

    Parameters:
        - G (networkx.Graph): Graph of residue pairs.

    Returns:
        - interfaces (list): List of interfaces, each as a dictionary with residue pairs.
    """
    interfaces = []
    for component in nx.connected_components(G):
        interface_pairs = [G.nodes[idx]['pair'] for idx in component]
        interface = {
            'residue_pairs': interface_pairs,
            'size': len(interface_pairs)
            # Additional attributes can be added here if needed
        }
        interfaces.append(interface)
    return interfaces

def interface_confidence(avg_pae, min_pae):
    """
    Determine the confidence level of an interface based on the average PAE.

    Parameters:
        - avg_pae (float): Average PAE of the interface.
        - min_pae (float): Minimum PAE of the interface.
    """
    if not avg_pae or not min_pae:
        return None
    if avg_pae < 7.5 and min_pae < 3:
        return 'high'
    if avg_pae < 15 and min_pae < 5:
        return 'medium'
    return 'low'

def find_and_score_interfaces(folder_path, pair_distance=6.0, interface_separation_distance = 5.0, interprotein_pae=15.0, intraprotein_pae=15.0,  all_atom=True):
    """
    Main function that selects and scores the best model.

    Parameters:
        - folder_path (str): Path to the folder containing the model files.
        - pair_distance (float): Maximum distance between residues from different proteins to be considered as part of an interface.
        - interface_separation_distance (float): Maximum distance between pair midpoints to be considered as part of the same interface.
        - interprotein_pae (float): Maximum PAE between residues from different proteins for the interaction to be considered as confident
        - intraprotein_pae (float): Maximum PAE between residues within the same protein for them to be considered as part of the same interface

    Returns:
        - interfaces (list): List of dictionaries containing data for each interface.
    """
    # Find folder name from folder path
    folder = folder_path.split('/')[-1]
    # Find chain groupings from folder name
    if '_output' in folder: # ignore suffix '_output' if present
        folder = folder.split('_output')[0]
    # Split into the two protein parts
    parts = folder.split('+')
    if len(parts) < 2: # naming system used in FL screen in AF3 is fold_protein1_fl_protein2_fl
        parts = folder.split('_')
        parts = [p for p in parts if p.lower() not in ['fold', 'fl']]

    # Check each protein part for 'dimer'
    p1_is_dimer = '_dimer_' in parts[0] or (parts[1] == 'dimer' and not len(parts) == 2)
    p2_is_dimer = '_dimer_' in parts[1]
    if len(parts) >= 3:
        if parts[2] == 'dimer':
            p2_is_dimer = True
    if "dimer" in folder and not (p1_is_dimer or p2_is_dimer):
        if "fold" not in folder: # ok in case of af3 predictions of single diimerised fragment (naming format = fold_p1_fl_dimer)
            raise ValueError(f"Folder name '{folder}' indicates a dimer, but no dimer detected in parts: {parts}. Please check code and folder naming conventions to make sure dimers are correctly interpreted.")
    # if dimer is in name twice
    elif folder.count('dimer') > 1 and not (p1_is_dimer and p2_is_dimer):
        raise ValueError(f"Folder name '{folder}' indicates two dimers, but only one dimer detected in parts: {parts}. Please check code and folder naming conventions to make sure dimers are correctly interpreted.")
    elif folder.count('dimer') == 1 and (p1_is_dimer and p2_is_dimer):
        raise ValueError(f"Folder name '{folder}' indicates a single dimer, but both parts detected as dimers: {parts}. Please check code and folder naming conventions to make sure dimers are correctly interpreted.")
    elif folder.count('dimer') > 2:
        raise ValueError(f"Folder name '{folder}' indicates more than two dimers - not yet supported.")
    elif folder.count('dimer') == 0 and (p1_is_dimer or p2_is_dimer):
        raise ValueError(f"Folder name '{folder}' indicates no dimers, but one or both parts detected as dimers: {parts}. Please check code and folder naming conventions to make sure dimers are correctly interpreted.")

    # Determine chain groupings based on dimer presence
    chain_groupings = None
    if not p1_is_dimer and not p2_is_dimer: # p1_Fx+p2_Fy
        chain_groupings = [('A'), ('B')]
    elif p1_is_dimer and not p2_is_dimer: # p1_dimer_Fx+p2_Fy
        chain_groupings = [('A', 'B'), ('C')]
    elif not p1_is_dimer and p2_is_dimer: # p1_Fx+p2_dimer_Fy
        chain_groupings = [('A'), ('B', 'C')]
    elif p1_is_dimer and p2_is_dimer: # p1_dimer_Fx+p2_dimer_Fy
        chain_groupings = [('A', 'B'), ('C', 'D')]
    else:
        raise ValueError(f"Unrecognized folder format: {folder}. Expected format: p1_Fx+p2_Fy, p1_dimer_Fx+p2_Fy, p1_Fx+p2_dimer_Fy, or p1_dimer_Fx+p2_dimer_Fy.")

    # Collect model data (only essential data)
    model_data = collect_model_data(folder_path, pair_distance, interprotein_pae, all_atom, chain_groupings=chain_groupings)

    # Calculate ROP scores
    for model in model_data:
        other_model_pairs = [m['secondary_pairs'] for m in model_data if m != model]
        model['rop'] = calculate_rop_score(model['confident_pairs'], other_model_pairs)
        model['avg_pct_rop'] = np.mean(calculate_percent_rops(model['confident_pairs'], other_model_pairs))
        # new - %rop
        #for other_model in other_model_pairs:
        #    shared_pairs = model['confident_pairs'] & set(other_model)
        #    percentage = len(shared_pairs) / len(model['confident_pairs'])
        #    print(percentage)

    # Select the best model
    best_model = select_best_model(model_data)
    other_model_pairs = [set(m['secondary_pairs']) for m in model_data if m != best_model]

    # find num chains in model
    chains = [i for i in best_model['structure_model'].get_chains()]
    if len(chains) == 2:
        # Compute pdockq
        best_model['pdockq'], _ = compute_pdockq(best_model['model_file'], best_model['json_file'])  
    # Extract iptm
    best_model['iptm'] = extract_iptm(best_model['log_file'], best_model['model_rank'])

    if not best_model['residue_pairs']:
        return best_model

    # Identify interfaces
    inv_abs_res_lookup_dict = {v: k for k, v in best_model['abs_res_lookup_dict'].items()}
    G = build_residue_pair_graph_with_spatial_proximity_and_pae(
        best_model['residue_pairs'],
        inv_abs_res_lookup_dict,
        best_model['structure_model'],
        pae_data=best_model['pae_data'],  # Include PAE data
        pae_cutoff=intraprotein_pae,
        distance_cutoff=interface_separation_distance
    )
    interfaces = get_interfaces_from_graph(G)

    interfaces = [interface for interface in interfaces if interface['size'] > 1]
    for interface in interfaces:
        #print(f"new_interface: {interface['residue_pairs']}")
        # update interface with data from best_model - EXCEPT 'residue_pairs' as this will override interface residue pairs if left in
        interface.update({k: v for k, v in best_model.items() if k not in ['residue_pairs']})
        interface['avg_pae'] = compute_average_interface_pae(interface['residue_pairs'], best_model['pae_data'])
        interface['confidence_class'] = interface_confidence(interface['avg_pae'], interface['min_pae'])
        interface['min_pae'] = residue_pairs_min_pae(interface['residue_pairs'], best_model['pae_data'])
        interface['location'] = readable_ranges(interface['residue_pairs'], best_model['structure_model'], best_model['abs_res_lookup_dict'])
        interface['rop'] = calculate_rop_score(interface['residue_pairs'], other_model_pairs, threshold=0.7)
        if interface['confidence_class'] in ['high', 'medium']:
            interface['evenness'] = compute_pae_evenness(interface['residue_pairs'], best_model['pae_data'])
            interface['location'] = readable_ranges(interface['residue_pairs'], best_model['structure_model'], best_model['abs_res_lookup_dict'])
            #interface['charge_match'] = charge_score(interface['residue_pairs'], best_model['structure_model'], inv_abs_res_lookup_dict)
            #interface['hydrophobicity_match'] = hydrophobicity_score(interface['residue_pairs'], best_model['structure_model'], inv_abs_res_lookup_dict)

    return interfaces

##################################################################################################
# # Example usage
# # Process all predictions in a folder
# from run_folders import process_all_predictions
# folder_path = '/Users/poppy/Dropbox/PCM'
# process_all_predictions(folder_path, find_and_score_interfaces, output_file='PCM_interface_analysis.csv', ipTM_graphic=False)

# # single interface
# folder_path = '/Users/poppy/Dropbox/centriole_screen/Plk4_Ana2/Plk4_F1+Ana2_F2'
# folder_path = '/Users/poppy/Dropbox/PCM/Spd-2_GCP5/Spd-2_F4+GCP5_F2'
# interface_data = find_and_score_interfaces(folder_path)
# # Display the interfaces
# for idx, interface in enumerate(interface_data):
#     print(f"Interface {idx + 1}:")
#     print(f"  Number of Residue Pairs: {interface['size']}")
#     print(f"  Average PAE: {interface['avg_pae']:.2f}")
#     print(f"  Residue Pairs: {interface['residue_pairs']}")
#     print(f"  Confidence Class: {interface['confidence_class']}")
#     print(f"  ROP: {interface['rop']}")
#     print(f"  min_pae: {interface['min_pae']:.2f}")
#     print(f"  Location: {interface['location']}")
#     #if 'location' in interface:
#         #print(f"  Location: {interface['location']}")
#         #print(f"  ROP: {interface['rop']}")
#     print()

# folder_path = '/Users/poppy/Dropbox/Msps/Msps_TACC/Msps_F3+TACC_F4'
# folder_path = '/Users/poppy/Dropbox/centriole_screen/Plk4_Ana2/Plk4_F1+Ana2_F2'
# folder_path = '/Users/poppy/Dropbox/centriole_screen/Ana2_Sas-4/Ana2_F1+Sas-4_F3'
# folder_path = '/Users/poppy/Dropbox/centriole_screen/Sas-6_Sas-6/Sas-6_F1+Sas-6_F1'
# folder_path = '/Users/poppy/Dropbox/PCM/TACC_TACC/TACC_F4+TACC_F4'
# folder_path = '/Users/poppy/Dropbox/PCM/TACC_AurA/TACC_F3+AurA_F1'
# folder_path = '/Users/poppy/Dropbox/FZY/FZY_MAD1/FZY_F1+MAD1_F2'
# folder_path = '/Users/poppy/Dropbox/BUB1/BUB1_Grip71/BUB1_F1+Grip71_F2'
# folder_path = '/Users/poppy/Dropbox/PCM/Cnn_Cnn/Cnn_F2+Cnn_F2'
#folder_path = '/Users/poppy/Dropbox/atub_btub_Ana1/Ana1_BUB1/Ana1_F3+BUB1_F2'
# interface_data = find_and_score_interfaces(folder_path, pair_distance = 6, interface_separation_distance = 5
#                                            #intraprotein_pae=7
#                                            )
# # Display the interfaces
# for idx, interface in enumerate(interface_data):
#     print(f"Interface {idx + 1}:")
#     #print(f"  Number of Residue Pairs: {interface['size']}")
#     print(f"  Average PAE: {interface['avg_pae']:.2f}")
#     #print(f"  Residue Pairs: {interface['residue_pairs']}")
#     #print(f"  Confidence Class: {interface['confidence_class']}")
#     if 'location' in interface:
#         print(f"  Location: {interface['location']}")
#         print(f"  ROP: {interface['rop']:.2f}")
#         print(f"  Model Rank: {interface['model_rank']}")
#         #print(f"  avg_pct_rop: {interface['avg_pct_rop']:.2f}")
#         print(f"  min_pae: {interface['min_pae']:.2f}")
#         # print interface keys
#         #print(interface.keys())
#     print()

# from process_pae import visualize_pae_matrix, determine_chain_lengths
# from analysis_utility import extract_pae, find_rank_001_files, parse_structure_file
# structure_file, json_file, log_file, PAE_png, fasta_file = find_rank_001_files(folder_path)
# structure_model = parse_structure_file(structure_file)
# pae_matrix = extract_pae(json_file)
# chain_lengths = determine_chain_lengths(structure_model)
# visualize_pae_matrix(pae_matrix, chain_lengths)
# visualize_pae_matrix(pae_matrix[161:, 161:])
#visualize_pae_matrix(pae_matrix[397:443, 397:443])
#visualize_pae_matrix(pae_matrix[122+486:194+486, 122+486:194+486])

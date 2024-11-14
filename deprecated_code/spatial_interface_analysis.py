"""
Old version of interface analysis code.
Categorises interfaces based on spatial proximity of residue pairs
Can accidentally over classify interfaces if very close together, but (unlike previous sequence
based version,) correctly distinguishes interfaces involving eg opposite sides of a helix, by
using midpoints of residue pairs to cluster and ignoring proximity in sequence
"""
import networkx as nx
import numpy as np
from scipy.spatial import cKDTree
from iptm_visualisation import extract_iptm
from process_pae import compute_average_interface_pae, compute_pae_evenness, residue_pairs_min_pae
from pdockq_calc import compute_pdockq
from analyse_models import collect_model_data, select_best_model, calculate_rop_score
from analysis_utility import readable_ranges

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

def build_residue_pair_graph_with_spatial_proximity(residue_pairs, inv_abs_res_lookup_dict, structure_model, distance_cutoff=5.0):
    """
    Build a graph where nodes are residue pairs, and edges connect pairs
    that are spatially close based on their representative coordinates.

    Parameters:
        - residue_pairs (list of tuples): List of tuples of absolute residue IDs.
        - inv_abs_res_lookup_dict (dict): Mapping from abs_res_id to (chain_id, res_num).
        - structure_model (Bio.PDB.Model.Model): The protein model object.
        - distance_cutoff (float): Maximum spatial distance between residue pairs
          to be considered part of the same interface.

    Returns:
        - G (networkx.Graph): The graph with residue pairs as nodes.
    """
    G = nx.Graph()

    rep_coords = []
    for idx, (res1_id, res2_id) in enumerate(residue_pairs):
        coord1 = get_residue_coordinate(res1_id, inv_abs_res_lookup_dict, structure_model)
        coord2 = get_residue_coordinate(res2_id, inv_abs_res_lookup_dict, structure_model)
        rep_coord = (coord1 + coord2) / 2
        G.add_node(idx, pair=(res1_id, res2_id), coord=rep_coord)
        rep_coords.append(rep_coord)

    # Build edges based on spatial proximity
    rep_coords_array = np.array(rep_coords)
    tree = cKDTree(rep_coords_array)
    pairs = tree.query_pairs(r=distance_cutoff)
    for i, j in pairs:
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

def interface_confidence(pae):
    """
    Determine the confidence level of an interface based on the average PAE.

    Parameters:
        - pae (np.array): The avg pae value for the interface
    """
    if not pae:
        return None
    if pae > 15:
        return 'low'
    if pae > 7:
        return 'medium'
    return 'high'

def find_and_score_interfaces(folder_path, distance_cutoff=5.0, pae_cutoff=15.0, all_atom=True):
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
    other_model_pairs = [set(m['confident_pairs']) for m in model_data if m != best_model]

    # Compute pdockq
    best_model['pdockq'], _ = compute_pdockq(best_model['model_file'], best_model['json_file'])
    # Extract iptm
    best_model['iptm'] = extract_iptm(best_model['json_file'], best_model['model_rank'])

    if not best_model['residue_pairs']:
        return best_model

    # Identify interfaces
    inv_abs_res_lookup_dict = {v: k for k, v in best_model['abs_res_lookup_dict'].items()}
    G = build_residue_pair_graph_with_spatial_proximity(
        best_model['residue_pairs'],
        inv_abs_res_lookup_dict,
        best_model['structure_model'],
        distance_cutoff=5.0  # Adjust the cutoff as needed
    )
    interfaces = get_interfaces_from_graph(G)

    interfaces = [interface for interface in interfaces if interface['size'] > 1]
    for interface in interfaces:
        # update interface with data from best_model - EXCEPT 'residue_pairs' as this will override interface residue pairs if left in
        interface.update({k: v for k, v in best_model.items() if k in ['pdockq', 'iptm', 'model_rank']})
        interface['avg_pae'] = compute_average_interface_pae(interface['residue_pairs'], best_model['pae_data'])
        interface['confidence_class'] = interface_confidence(interface['avg_pae'])
        if interface['confidence_class'] in ['high', 'medium']:
            interface['min_pae'] = residue_pairs_min_pae(interface['residue_pairs'], best_model['pae_data'])
            interface['rop'] = calculate_rop_score(interface['residue_pairs'], other_model_pairs, threshold=0.7)
            interface['evenness'] = compute_pae_evenness(interface['residue_pairs'], best_model['pae_data'])
            interface['location'] = readable_ranges(interface['residue_pairs'], best_model['structure_model'], best_model['abs_res_lookup_dict'])


    return interfaces

##############################################################################################
# Example Usage
#folder_path = '/Users/poppy/Dropbox/centriole_screen/Plk4_Ana2/Plk4_F1+Ana2_F2'
#interface_data = find_and_score_interfaces(folder_path) # Returns a list of interfaces, each as a dictionary
# Display the interfaces
# for idx, interface in enumerate(interface_data):
#     print(f"Interface {idx + 1}:")
#     print(f"  Number of Residue Pairs: {interface['size']}")
#     print(f"  Average PAE: {interface['avg_pae']:.2f}")
#     print(f"  Residue Pairs: {interface['residue_pairs']}")
#     print(f"  Confidence Class: {interface['confidence_class']}")
#     if 'location' in interface:
#         print(f"  Location: {interface['location']}")
#         print(f"  ROP: {interface['rop']:.2f}")
#     print()
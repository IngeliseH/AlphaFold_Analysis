"""
Old version of interface analysis code.
Categorises interfaces based on sequence proximity of residue pairs
Eg - (eg pairs (1,2), (6,4), (7,11) exist between proteinA proteinB, proteinA 6 and 7 are
continuous, proteinB 2 and 4 are nearby enough to count as continuous, so 1, 6-7 of ProteinA is
involved in a single interface with 2-4, 11 of ProteinB

Not ideal as frequently combines interfaces - eg if two different sides of a helix form separate
interfaces, they will be combined into one by sequence proximity
"""
import networkx as nx
from iptm_visualisation import extract_iptm
from process_pae import compute_average_interface_pae, residue_pairs_min_pae, compute_pae_evenness
from pdockq_calc import compute_pdockq
from analyse_models import collect_model_data, select_best_model, calculate_rop_score
from analysis_utility import readable_ranges


def build_residue_pair_graph_with_sequence_proximity(residue_pairs, sequence_distance=5):
    """
    Build a graph where nodes are residue pairs, and edges connect pairs
    that share a residue or have residues within sequence_distance on the same protein.
    
    Parameters:
        - residue_pairs (list of tuples): List of tuples where the first item is
          a residue from protein1 and the second is from protein2.
        - sequence_distance (int): Maximum sequence distance within which residues
          are considered to be part of the same interface.

    Returns:
        - G (networkx.Graph): The graph with residue pairs as nodes.
    """
    if not isinstance(residue_pairs, list):
        residue_pairs = list(residue_pairs)

    G = nx.Graph()
    
    # Add nodes to the graph
    for idx, pair in enumerate(residue_pairs):
        G.add_node(idx, pair=pair)
    
    # Build edges based on shared residues or sequence proximity
    for i in range(len(residue_pairs)):
        res1_i, res2_i = residue_pairs[i]  # protein1 residue, protein2 residue

        for j in range(i + 1, len(residue_pairs)):
            res1_j, res2_j = residue_pairs[j]  # next pair's protein1 and protein2 residues

            # Check if residue pairs share a residue
            if res1_i == res1_j or res2_i == res2_j:
                G.add_edge(i, j)
                continue

            # Check sequence proximity within the same protein
            if abs(res1_i - res1_j) <= sequence_distance:  # Within protein1
                G.add_edge(i, j)
                continue
            if abs(res2_i - res2_j) <= sequence_distance:  # Within protein2
                G.add_edge(i, j)
                continue

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
            'size': len(interface_pairs),
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
    G = build_residue_pair_graph_with_sequence_proximity(
        best_model['residue_pairs'],
        sequence_distance=1
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

####################################################################################################
# Example usage
# folder_path = '/Users/poppy/Dropbox/centriole_screen/Plk4_Ana2/Plk4_F1+Ana2_F2'
# interfaces = find_and_score_interfaces(folder_path)

# # Display the interfaces
# for idx, interface in enumerate(interfaces):
#     print(f"Interface {idx + 1}:")
#     print(f"  Number of Residue Pairs: {interface['size']}")
#     print(f"  Average PAE: {interface['avg_pae']:.2f}")
#     if 'location' in interface:
#         print(f"  Location: {interface['location']}")
#     print()

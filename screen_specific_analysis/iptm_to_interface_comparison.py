import pandas as pd

# read in csvs in data folder - iptm and interface networks - no headings but should be give names weight, domain1, domain2
iptm_network = pd.read_csv('data/iptm_network.csv', header=None)
interface_network = pd.read_csv('data/interface_network.csv', header=None)
iptm_network.columns = ['weight', 'domain1', 'domain2']
interface_network.columns = ['weight', 'domain1', 'domain2']

# Define a function to analyze the agreement and exclusive rates with given weight cutoff
def compare_networks(network1, network2, cutoff1=None, cutoff2=None):
    # Apply weight cutoff if specified
    if cutoff1 is not None:
        network1 = network1[network1['weight'] >= cutoff1]
    if cutoff2 is not None:
        network2 = network2[network2['weight'] >= cutoff2]
    
    # Create sets of edges
    edges1 = set(zip(network1['domain1'], network1['domain2']))
    edges2 = set(zip(network2['domain1'], network2['domain2']))
    
    # Calculate metrics
    agreement = len(edges1.intersection(edges2))
    only_in_network1 = len(edges1 - edges2)
    only_in_network2 = len(edges2 - edges1)
    total_edges = len(edges1.union(edges2))
    agreement_percentage = (agreement / total_edges) * 100 if total_edges > 0 else 0
    only_in_network1_percentage = (only_in_network1 / total_edges) * 100 if total_edges > 0 else 0
    only_in_network2_percentage = (only_in_network2 / total_edges) * 100 if total_edges > 0 else 0
    
    return {
        'agreement%': agreement_percentage,
        'n1%': only_in_network1_percentage,
        'n2%': only_in_network2_percentage,
        'total': total_edges
    }

# Define various cutoff combinations to test
cutoffs = [(None, None), (None, 2), (None, 3), (None, 4), (2, None), (2, 2), (2, 3), (2, 4), (3, None), (3, 2), (3, 3), (3, 4)]

# Run comparisons and collect results
results = []
for cutoff1, cutoff2 in cutoffs:
    result = compare_networks(iptm_network, interface_network, cutoff1, cutoff2)
    result['n1_cutoff'] = cutoff1
    result['n2_cutoff'] = cutoff2
    results.append(result)

# Convert results to DataFrame for better display
results_df = pd.DataFrame(results)
print(results_df)

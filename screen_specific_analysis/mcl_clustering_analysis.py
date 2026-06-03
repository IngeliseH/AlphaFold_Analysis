# assess clustering of proteins
# trying to assess whether PPIs are seen more within group than between group
# proxy for this is determining whether known groups are clustered together
import pandas as pd
import numpy as np
from sklearn.metrics import pair_confusion_matrix
import matplotlib.pyplot as plt
import markov_clustering as mc
import scipy.sparse
import networkx as nx

interface_data = pd.read_csv('/Users/Poppy/Dropbox/t7_interface_analysis_2026.04.11.csv')
interface_data = pd.read_csv('/Users/Poppy/Dropbox/t7_interface_analysis_with_coil.csv')

protein_data = pd.read_csv('/Users/poppy/Dropbox/all_fragments.csv')
# map protein names to true groups
protein_groups = protein_data.set_index('name')['expected_cluster']
# identify any unlabelled proteins
missing_labels = protein_groups[protein_groups.isna()]
if not missing_labels.empty:
    print("Proteins with missing true group labels:")
    print(missing_labels)
# make any group containing 'kinetochore' into 'kinetochore'
#protein_groups = protein_groups.apply(lambda x: 'kinetochore' if 'kinetochore' in str(x).lower() else x)
# combine 'centriole' and 'PCM' groups into 'centrosome'
#protein_groups = protein_groups.apply(lambda x: 'centrosome' if str(x).lower() in ['centriole', 'pcm'] else x)

# make a matrix of presence/absence of interfaces between proteins
proteins = protein_groups.index
interface_matrix = pd.DataFrame(0, index=proteins, columns=proteins)
interface_matrix = interface_matrix.astype(float)
min_pae_cutoff = 11

filtered_data = interface_data[
    (interface_data['size'] >= 15) &
    (interface_data['rop'] >= 3) &
    (interface_data['avg_pae'] <= 7.5) &
    (interface_data['min_pae'] <= min_pae_cutoff) &
    (interface_data['max_promiscuity'] <= 20) &
    (interface_data['coiledness_score'] <= 0.25) &
    (interface_data['pdockq'] >= 0.23)
]

for _, row in filtered_data.iterrows():
    p1 = row['Protein1']
    p2 = row['Protein2']
    if p1 in proteins and p2 in proteins:
        original_val = interface_matrix.loc[p1, p2]
        #normalise min_pae - lower min_pae means stronger interaction, so invert this to make higher values mean stronger interaction, and normalise between 0-1 with cutoff at min_pae_cutoff
        new_val = max(0, min_pae_cutoff - row['min_pae']) / min_pae_cutoff
        interface_matrix.loc[p1, p2] = max(original_val, new_val)
        interface_matrix.loc[p2, p1] = max(original_val, new_val)  # symmetric

"""
MCL_INFLATION_SWEEP = [1.1, 1.2, 1.5, 2, 2.5, 3, 4]
print("MCL inflation parameter sweep:")
for infl in MCL_INFLATION_SWEEP:
    result = mc.run_mcl(scipy.sparse.csr_matrix(interface_matrix.values), inflation=infl)
    sweep_clusters = mc.get_clusters(result)
    print(f"Inflation={infl:.1f}: n_clusters={len(sweep_clusters)}, largest_cluster={max((len(c) for c in sweep_clusters), default=0)}")

MCL_INFLATION = 1.5  # higher => smaller/tighter MCL clusters
"""

# align to valid true labels (avoid NaN and index mismatch in metrics)
true_labels = protein_groups.loc[interface_matrix.index]
valid_mask = true_labels.notna()
if (~valid_mask).any():
    print(f"Dropping {(~valid_mask).sum()} proteins with missing true labels for MCL evaluation.")
true_labels = true_labels.loc[valid_mask]
mcl_matrix = interface_matrix.loc[valid_mask, valid_mask]

sparse_matrix = scipy.sparse.csr_matrix(mcl_matrix.values)
#result = mc.run_mcl(sparse_matrix, inflation=MCL_INFLATION)
result = mc.run_mcl(sparse_matrix)
clusters = mc.get_clusters(result)

# assign predicted MCL cluster per protein BEFORE plotting
mcl_pred_by_protein = pd.Series(-1, index=mcl_matrix.index, name="mcl_pred_cluster")
for i, cluster in enumerate(clusters):
    cluster_idx = list(cluster)
    cluster_names = mcl_matrix.index.take(cluster_idx).tolist()
    print(f"Cluster {i}: {cluster_names}")
    mcl_pred_by_protein.loc[cluster_names] = i

# pairwise metrics
tn, fp, fn, tp = pair_confusion_matrix(true_labels, mcl_pred_by_protein).ravel()
rand_index = (tp + tn) / (tp + tn + fp + fn)
jaccard_index = tp / (tp + fp + fn) if (tp + fp + fn) else float("nan")
print(f"\nMCL clustering:")
print(f"Pair counts: a=SS={tp}, b=SD={fp}, c=DS={fn}, d=DD={tn}")
print(f"Rand index: {rand_index:.4f}")
print(f"Jaccard index: {jaccard_index:.4f}")

# network plot coloured by MCL clusters
G = nx.from_pandas_adjacency(mcl_matrix)
node_colors = []
i = 1
cluster_id_colour_map = {}
for node in G.nodes():
    cluster_id = mcl_pred_by_protein.loc[node]
    if cluster_id in cluster_id_colour_map:
        node_colors.append(cluster_id_colour_map[cluster_id])
        continue
    if cluster_id < 0:
        node_colors.append("#BDBDBD")  # grey for unclustered
    else:
        cluster_size = (mcl_pred_by_protein == cluster_id).sum()
        if cluster_size <= 1:
            node_colors.append("#BDBDBD")  # grey for clusters of size 1
        else:
            node_colors.append(f"C{i % 10}")  # color by cluster ID, cycling through matplotlib default colors
            cluster_id_colour_map[cluster_id] = f"C{i % 20}"
            i += 1
unique_clusters = set(mcl_pred_by_protein[mcl_pred_by_protein >= 0].unique())
plt.figure(figsize=(8, 8))
pos = nx.spring_layout(G, k=0.5, seed=0)
nx.draw(G, pos, with_labels=True, node_color=node_colors, node_size=100, font_size=8, edge_color='lightgray')
plt.title(f"Protein Interaction Network Colored by MCL Clusters")
plt.show()

# get pair confusion matrix and rand and jaccard indices for MCL clustering vs true groups
true_labels = protein_groups.loc[interface_matrix.index]
mcl_pred_by_protein = pd.Series(-1, index=interface_matrix.index, name="mcl_pred_cluster")
for i, cluster in enumerate(clusters):
    cluster_idx = list(cluster)  # tuple -> list for safe positional indexing
    cluster_names = interface_matrix.index.take(cluster_idx).tolist()
    mcl_pred_by_protein.loc[cluster_names] = i
tn, fp, fn, tp = pair_confusion_matrix(true_labels, mcl_pred_by_protein).ravel()
rand_index = (tp + tn) / (tp + tn + fp + fn)
jaccard_index = tp / (tp + fp + fn) if (tp + fp + fn) else float("nan")
print(f"\nMCL clustering:")
print(f"Pair counts: a=SS={tp}, b=SD={fp}, c=DS={fn}, d=DD={tn}")
print(f"Rand index: {rand_index:.4f}")
print(f"Jaccard index: {jaccard_index:.4f}")

# plot matrix of clustering for each group
# matrix split by category
categories = true_labels.unique()
# find proportion of pairs in same predicted cluster within each true group
category_within_group_scores = {}
for category in categories:
    category_proteins = true_labels[true_labels == category].index
    if len(category_proteins) < 2:
        print(f"Skipping category '{category}' with less than 2 proteins for within-group clustering score.")
        continue
    same_cluster_count = sum(
        1 for i in range(len(category_proteins)) for j in range(i + 1, len(category_proteins))
        if mcl_pred_by_protein.loc[category_proteins[i]] == mcl_pred_by_protein.loc[category_proteins[j]]
    )
    total_pairs = len(category_proteins) * (len(category_proteins) - 1) / 2
    category_within_group_scores[category] = same_cluster_count / total_pairs if total_pairs > 0 else float("nan")
# sort categories by proportion of pairs in same predicted cluster
categories = sorted(categories, key=lambda cat: category_within_group_scores.get(cat, 0), reverse=True)
category_matrix = pd.DataFrame(0.0, index=categories, columns=categories)

# for every pair of proteins between these two categories, count how many are classified in the same category
# get proportion by dividing by total number of pairs between these two categories
from itertools import combinations
for cat1 in categories:
    for cat2 in categories:
        if cat1 == cat2:
            pairs = {comb for comb in combinations(true_labels[true_labels == cat1].index, 2)}
        else:
            pairs = {(p1, p2) for p1 in true_labels[true_labels == cat1].index for p2 in true_labels[true_labels == cat2].index}
        if pairs:
            same_cluster_count = sum(1 for p1, p2 in pairs if mcl_pred_by_protein[p1] == mcl_pred_by_protein[p2])
            category_matrix.loc[cat1, cat2] = same_cluster_count / len(pairs)

# keep only lower triangle
for i in range(len(category_matrix)):
    for j in range(i + 1, len(category_matrix)):
        category_matrix.iat[i, j] = np.nan

# plot
plt.figure(figsize=(8, 6))
im = plt.imshow(category_matrix, cmap='viridis', vmin=0, vmax=1)
plt.colorbar(im, label='Proportion of pairs in same predicted cluster')
plt.xticks(range(len(categories)), categories, rotation=90)
plt.yticks(range(len(categories)), categories)
plt.title('Clustering agreement between true groups')
plt.tight_layout()
plt.show()
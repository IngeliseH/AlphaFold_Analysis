import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# =============================================
# 1. Load Data and compartment Lists
# =============================================
interactions = pd.read_csv("/Users/poppy/Dropbox/all_interface_analysis_2025.02.26.csv")

compartment_proteins = {
    'centriole':['alpha-tubulin', 'Ana1', 'Ana2', 'Ana3', 'Asl', 'beta-tubulin', 'Centrin',
                 'Centrobin', 'Cep135', 'Cep97', 'CP110', 'Ctp', 'PLK4', 'Poc1',
                 'Rcd4', 'Sas-4', 'Sas-6'],
    'PCM':['Asp', 'AurA', 'CHC', 'Cnn', 'Msps', 'Mud', 'PLK1', 'PLP', 'Spd-2', 'TACC'],
    'gTURC':['gamma-tubulin', 'GCP2', 'GCP3', 'GCP4', 'GCP5', 'Gcp6', 'Grip71', 'Mzt1'],
    'cdks':['Cdk1', 'Cdk2', 'CycA', 'CycB1', 'CycB3', 'CycE'],
    'kinetochore':['AurB', 'BUB1', 'BUBR1', 'CENPC', 'CID', 'FZY', 'MAD1', 'MAD2', 'MPS1',
                   'ROD', 'ZW10', 'ZWILCH']
}

# =============================================
# 2. Map Proteins to compartments & Remove Duplicates
# =============================================
protein_to_compartment = {}
for compartment, proteins in compartment_proteins.items():
    for protein in proteins:
        protein_to_compartment[protein] = compartment

interactions["compartment1"] = interactions["Protein1"].map(protein_to_compartment)
interactions["compartment2"] = interactions["Protein2"].map(protein_to_compartment)
interactions = interactions.dropna(subset=["compartment1", "compartment2"])

# Sort compartments to ensure upper triangular matrix
interactions[["compartment1", "compartment2"]] = np.sort(interactions[["compartment1", "compartment2"]], axis=1)

# Remove self-interactions and condense pairs
interactions = interactions[
    (interactions["Protein1"] != interactions["Protein2"])
].reset_index(drop=True)

# =============================================
# 3. Filter Interactions for All Methods
# =============================================
def remove_duplicates(df):
    df[["Protein1", "Protein2"]] = np.sort(df[["Protein1", "Protein2"]], axis=1)
    return df.drop_duplicates(subset=["Protein1", "Protein2"])

# Accepted method
iptm_threshold = 0.55
accepted_interactions = interactions[interactions["iptm"] > iptm_threshold]
accepted_interactions = remove_duplicates(accepted_interactions)

# Higher threshold
higher_threshold = 0.7
stronger_interactions = interactions[interactions["iptm"] > higher_threshold]
stronger_interactions = remove_duplicates(stronger_interactions)

# Custom method
criteria = {
    'rop': lambda x: x >= 2,
    'avg_pae': lambda x: x <= 10,
    'size': lambda x: x >= 5,
    'p1d_promiscuity': lambda x: x <= 30,
    'p2d_promiscuity': lambda x: x <= 30,
    'iptm': lambda x: x >= 0.35,
    'min_pae': lambda x: x <= 5
}

custom_interactions = interactions[
    interactions["rop"].apply(criteria["rop"]) &
    interactions["avg_pae"].apply(criteria["avg_pae"]) &
    interactions["size"].apply(criteria["size"]) &
    interactions["p1d_promiscuity"].apply(criteria["p1d_promiscuity"]) &
    interactions["p2d_promiscuity"].apply(criteria["p2d_promiscuity"]) &
    interactions["iptm"].apply(criteria["iptm"]) &
    interactions["min_pae"].apply(criteria["min_pae"])
]
custom_interactions = remove_duplicates(custom_interactions)

# =============================================
# 4. Generate Matrices with Size Normalization
# =============================================
def get_normalized_matrix(df, compartment_proteins):
    compartments = sorted(compartment_proteins.keys())
    size_dict = {org: len(proteins) for org, proteins in compartment_proteins.items()}
    
    # Get raw counts with proper dtype
    count_matrix = pd.crosstab(
        index=df["compartment2"],
        columns=df["compartment1"]
    ).reindex(index=compartments, columns=compartments, fill_value=0).astype(float)
    print(count_matrix)
    
    # Initialize possible_pairs with float dtype
    possible_pairs = pd.DataFrame(index=compartments, columns=compartments, dtype=float)
    
    for i in compartments:
        for j in compartments:
            n = size_dict[i]
            m = size_dict[j]
            if i == j:
                possible_pairs.loc[i,j] = n * (n-1) / 2  # Use float division
            else:
                possible_pairs.loc[i,j] = n * m
    
    # Handle division safely
    norm_matrix = count_matrix / possible_pairs.replace(0, np.nan)  # Prevent division by zero
    norm_matrix = norm_matrix.fillna(0)
    
    return count_matrix.astype(int), norm_matrix

accepted_counts, accepted_norm = get_normalized_matrix(accepted_interactions, compartment_proteins)
stronger_counts, stronger_norm = get_normalized_matrix(stronger_interactions, compartment_proteins)
custom_counts, custom_norm = get_normalized_matrix(custom_interactions, compartment_proteins)

# =============================================
# 5. Create Combined Annotation Labels
# =============================================
def create_annotations(count_mat, norm_mat):
    annot = pd.DataFrame(index=count_mat.index, columns=count_mat.columns)
    for i in count_mat.index:
        for j in count_mat.columns:
            count = count_mat.loc[i,j]
            norm = f"{norm_mat.loc[i,j]:.2f}" if norm_mat.loc[i,j] > 0 else "0"
            annot.loc[i,j] = f"{count}\n({norm})"
    return annot.values

accepted_annot = create_annotations(accepted_counts, accepted_norm)
stronger_annot = create_annotations(stronger_counts, stronger_norm)
custom_annot = create_annotations(custom_counts, custom_norm)

# =============================================
# 6. Plot Heatmaps with Dual Annotations
# =============================================
mask = np.triu(np.ones_like(custom_norm, dtype=bool), k=1)  # Mask UPPER triangle EXCLUDING diagonal

# heat map for lower triangle with diagonal
fig, axes = plt.subplots(1, 3, figsize=(24, 7))

# Find global maximum for consistent coloring
max_norm = max([accepted_norm.max().max(), custom_norm.max().max(), stronger_norm.max().max()])
# Accepted Method
sns.heatmap(
    accepted_norm,
    annot=accepted_annot,
    fmt="",
    cmap="YlGnBu",
    ax=axes[0],
    vmin=0,
    vmax=max_norm,
    cbar_kws={"label": "Size-Normalized Interaction Rate"},
    linewidths=0.5,
    annot_kws={"size": 9},
    mask=mask
)

axes[0].set_title("ipTM > 0.55\nRaw Counts + Normalized Rates")

# Custom Method
sns.heatmap(
    custom_norm,
    annot=custom_annot,
    fmt="",
    cmap="YlGnBu",
    ax=axes[1],
    vmin=0,
    vmax=max_norm,
    cbar_kws={"label": "Size-Normalized Interaction Rate"},
    linewidths=0.5,
    annot_kws={"size": 9},
    mask=mask
)
axes[1].set_title("Custom Method\nRaw Counts + Normalized Rates")

# Stronger Method
sns.heatmap(
    stronger_norm,
    annot=stronger_annot,
    fmt="",
    cmap="YlGnBu",
    ax=axes[2],
    vmin=0,
    vmax=max_norm,
    cbar_kws={"label": "Size-Normalized Interaction Rate"},
    linewidths=0.5,
    annot_kws={"size": 9},
    mask=mask
)
axes[2].set_title("ipTM > 0.7\nRaw Counts + Normalized Rates")

plt.tight_layout()
plt.show()

# =============================================
# 7. Create Difference Heatmaps
# =============================================
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Calculate differences
diff_accepted_custom = custom_counts - accepted_counts
diff_stronger_custom = custom_counts - stronger_counts

# Create colormap
#cmap = sns.diverging_palette(240, 10, as_cmap=True)  # Blue to Red with yellow midpoint
cmap = sns.diverging_palette(240, 10, as_cmap=True)

# Plot Accepted vs Custom
sns.heatmap(
    diff_accepted_custom,
    annot=True,
    fmt="d",
    cmap=cmap,
    ax=axes[0],
    center=0,
    cbar_kws={"label": "Change in Interaction Count\n(Custom - Accepted)"},
    linewidths=0.5,
    mask=mask,
    annot_kws={"size": 10}
)
axes[0].set_title("Change from Accepted to Custom Method\n(Blue = Less in Custom, Red = More in Custom)")

# Plot Stronger vs Custom
sns.heatmap(
    diff_stronger_custom,
    annot=True,
    fmt="d",
    cmap=cmap,
    ax=axes[1],
    center=0,
    cbar_kws={"label": "Change in Interaction Count\n(Custom - Stronger)"},
    linewidths=0.5,
    mask=mask,
    annot_kws={"size": 10}
)
axes[1].set_title("Change from Stronger to Custom Method\n(Blue = Less in Custom, Red = More in Custom)")

plt.tight_layout()
plt.show()

# =============================================
# 8. Make difference values proportional to accepted counts
# =============================================

fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Calculate differences
diff_accepted_custom = ((custom_counts-accepted_counts) / accepted_counts) * 100
diff_stronger_custom = ((custom_counts-stronger_counts) / stronger_counts) * 100

# Create colormap with 100 as diverging point
cmap = sns.diverging_palette(240, 10, as_cmap=True)

# Plot Accepted vs Custom
sns.heatmap(
    diff_accepted_custom,
    annot=True,
    fmt=".2f",
    cmap=cmap,
    ax=axes[0],
    center=1,
    cbar_kws={"label": "Percentage Change in Interaction Count\n(Change / Accepted * 100)"},
    linewidths=0.5,
    mask=mask,
    annot_kws={"size": 10}
)
axes[0].set_title("Percent change in interaction number from Accepted to Custom Method\n(Blue = Less in Custom, Red = More in Custom)")

# Plot Stronger vs Custom
sns.heatmap(
    diff_stronger_custom,
    annot=True,
    fmt=".2f",
    cmap=cmap,
    ax=axes[1],
    center=1,
    cbar_kws={"label": "Percentage Change in Interaction Count\n(Change / Stronger * 100)"},
    linewidths=0.5,
    mask=mask,
    annot_kws={"size": 10}
)
axes[1].set_title("Percent change in interaction number from Stronger to Custom Method\n(Blue = Less in Custom, Red = More in Custom)")

plt.tight_layout()
plt.show()

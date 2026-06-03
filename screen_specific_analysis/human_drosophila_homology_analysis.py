# prepping df listing all the protein pairs with interactions that appear in drosophila and human screens
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

#import and format human and drosophila interaction data
drosophila_df = pd.read_csv("/Users/poppy/Dropbox/t7_interface_analysis_2026.04.29.csv")
human_df = pd.read_csv("/Users/poppy/Dropbox/human_interface_analysis_with_coil_2026.04.30.csv")
human_df = human_df.rename(columns={"Source": "Protein1", "Target": "Protein2"})

# make list of all unique proteins in both datasets
drosophila_proteins = set(drosophila_df["Protein1"]).union(set(drosophila_df["Protein2"]))
human_proteins = set(human_df["Protein1"]).union(set(human_df["Protein2"]))

# map homologous proteins
homolog_dict = {'Ana1': ['CEP295'], 'CP110': ['CCP110'], 'Centrin': ['CETN2', 'CETN3'],
                'Centrobin': ['CNTROB'], 'Cep97': ['CEP97'], 'Cep135': ['CEP135'],
                'Fzr': ['FZR'], 'Mzt1': ['MZT1'], 'PLK1': ['PLK1'], 'PLK4': ['PLK4'],
                'Poc1': ['POC1A', 'POC1B'], 'Rcd4': ['PPP1R35'], 'Ana3': ['RTTN'], 'Sas-6': ['SASS6'],
                'Sas-4': ['CENPJ'], 'Ana2': ['STIL'], 'AurA': ['AURKA'], 'Cnn': ['CDK5RAP2'],
                'Asl': ['CEP152'], 'Spd-2': ['CEP192'], 'Msps': ['CKAP5'], 'Grip71': ['NEDD1'],
                'Plp': ['PCNT'], 'GCP2': ['TUBGCP2'], 'GCP3': ['TUBGCP3'], 'GCP4': ['TUBGCP4'],
                'GCP5': ['TUBGCP5'], 'Gcp6': ['TUBGCP6'], 'Mud': ['NUMA1'], 'gamma-tubulin': ['TUBG1']}
#print("Remaining drosophila proteins without homologs:")
#print([prot for prot in drosophila_proteins if prot not in homolog_dict])
#print("Remaining human proteins without homologs:")
#print([prot for prot in human_proteins if prot not in [item for sublist in homolog_dict.values() for item in sublist]])

# filter dataframes
criteria = {
    "min_pae": 11, #must be <=
    "avg_pae": 14, #must be <=
    "rop": 1, #must be >=
    "size": 15, #must be >=
    #"max_promiscuity": 20, #must be <=
    #"pdockq": 0.23, #must be >=
}
dm_proteins_with_homologs = set(homolog_dict.keys())
filtered_drosophila_df = drosophila_df[
    (drosophila_df["min_pae"] <= criteria["min_pae"])
    & (drosophila_df["avg_pae"] <= criteria["avg_pae"])
    & (drosophila_df["rop"] >= criteria["rop"])
    & (drosophila_df["size"] >= criteria["size"])
    & (drosophila_df["Protein1"].isin(dm_proteins_with_homologs))
    & (drosophila_df["Protein2"].isin(dm_proteins_with_homologs))
    #& (drosophila_df["max_promiscuity"] <= criteria["max_promiscuity"])
    #& (drosophila_df["pdockq"] >= criteria["pdockq"])
]
hs_proteins_with_homologs = {hs for hs_list in homolog_dict.values() for hs in hs_list}
filtered_human_df = human_df[
    (human_df["min_pae"] <= criteria["min_pae"])
    & (human_df["avg_pae"] <= criteria["avg_pae"])
    & (human_df["rop"] >= criteria["rop"])
    & (human_df["size"] >= criteria["size"])
    & (human_df["Protein1"].isin(hs_proteins_with_homologs))
    & (human_df["Protein2"].isin(hs_proteins_with_homologs))
    #& (human_df["max_promiscuity"] <= criteria["max_promiscuity"])
    #& (human_df["pdockq"] >= criteria["pdockq"])
]
print(f"Number of interactions in filtered drosophila dataset: {len(filtered_drosophila_df)}")
print(f"Number of interactions in filtered human dataset: {len(filtered_human_df)}")

# make a new df listing number of human and drosophila interactions seen for each protein pair
interaction_data = []
# check every drosophila protein pairs
for prot1 in dm_proteins_with_homologs:
    for prot2 in dm_proteins_with_homologs:
        if prot1 < prot2:
            dm_interactions = 0
            dm_interactions += len(filtered_drosophila_df[
                ((filtered_drosophila_df["Protein1"] == prot1) & (filtered_drosophila_df["Protein2"] == prot2))
                | ((filtered_drosophila_df["Protein1"] == prot2) & (filtered_drosophila_df["Protein2"] == prot1))
            ])
            # if no interactions, check with protein names in other order
            if dm_interactions == 0:
                dm_interactions += len(filtered_drosophila_df[
                    ((filtered_drosophila_df["Protein1"] == prot2) & (filtered_drosophila_df["Protein2"] == prot1))
                ])

            # find human homologs for each protein
            hs_prot1_candidates = homolog_dict.get(prot1, [])
            hs_prot2_candidates = homolog_dict.get(prot2, [])
            print(f"Checking interactions for {prot1} and {prot2}, human homologs are {hs_prot1} and {hs_prot2}")
            hs_total_interactions = 0
            for hs_prot1 in hs_prot1_candidates:
                for hs_prot2 in hs_prot2_candidates:
                    hs_interactions = 0
                    if hs_prot1 is not None and hs_prot2 is not None:
                        print(f"Checking human interactions for {hs_prot1} and {hs_prot2}")
                        hs_interactions += len(filtered_human_df[
                            ((filtered_human_df["Protein1"] == hs_prot1) & (filtered_human_df["Protein2"] == hs_prot2))
                            | ((filtered_human_df["Protein1"] == hs_prot2) & (filtered_human_df["Protein2"] == hs_prot1))
                        ])
                        if hs_interactions == 0:
                            hs_interactions += len(filtered_human_df[
                                ((filtered_human_df["Protein1"] == hs_prot2) & (filtered_human_df["Protein2"] == hs_prot1))
                            ])
                        hs_total_interactions += hs_interactions
            interaction_data.append({
                "Protein1": prot1,
                "Protein2": prot2,
                "hs_interactions": hs_total_interactions,
                "dm_interactions": dm_interactions,
            })
interaction_df = pd.DataFrame(interaction_data)
interaction_df["Interaction_Type"] = interaction_df.apply(lambda row: "Both" if row["hs_interactions"] > 0 and row["dm_interactions"] > 0 else ("Human" if row["hs_interactions"] > 0 else ("Drosophila" if row["dm_interactions"] > 0 else "None")), axis=1)
# drop all rows where neither human or drosophila protein pairs have interactions
interaction_df = interaction_df[interaction_df["Interaction_Type"] != "None"]


# make network plot
G = nx.Graph()
for index, row in interaction_df.iterrows():
    prot1 = row["Protein1"]
    prot2 = row["Protein2"]
    interaction_type = row["Interaction_Type"]
    if interaction_type == "Both":
        G.add_edge(prot1, prot2, color="green")
    elif interaction_type == "Human":
        G.add_edge(prot1, prot2, color="orange")
    elif interaction_type == "Drosophila":
        G.add_edge(prot1, prot2, color="blue")
edge_colors = [G[u][v]["color"] for u, v in G.edges()]
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True, edge_color=edge_colors)
plt.show()

#################################################
# Add interface contributor residues to the dataframes for later use
def find_interface_contributors(interface_data):
    """Get interface contributor residues from each side for each interface."""
    interface_data = interface_data.copy()
    interface_data['residue_pairs'] = interface_data['residue_pairs'].apply(lambda x: literal_eval(x) if isinstance(x, str) else x)
    interface_data['chain_a_contributors'] = interface_data['residue_pairs'].apply(lambda pairs: set(pair[0] for pair in pairs))
    interface_data['chain_b_contributors'] = interface_data['residue_pairs'].apply(lambda pairs: set(pair[1] for pair in pairs))
    return interface_data

filtered_drosophila_df = find_interface_contributors(filtered_drosophila_df)
filtered_human_df = find_interface_contributors(filtered_human_df)

################################################
# initial method, using sequence homology
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from ast import literal_eval

def map_fragments_to_seqs(protein_data_csv, proteins_to_include=None):
    """Reads a CSV file containing protein fragment information and creates a mapping of fragment names to their sequences."""
    protein_data = pd.read_csv(protein_data_csv)
    seq_dict = {}
    for _, row in protein_data.iterrows():
        if 'dimer' in row['name'].lower():
            # skip dimers as sequence layout and pair numbering is different
            continue
        if proteins_to_include and row['name'] not in proteins_to_include:
            continue
        name = row['name']
        fragment_sequences = literal_eval(row['fragment_sequences'])
        for i, seq in enumerate(fragment_sequences):
            seq_dict[f'{name}_F{i+1}'] = seq
    return seq_dict

def get_alignment_mapping(seq1, seq2):
    """Aligns two sequences and maps 1-based indices.

    Uses local alignment plus a substitution matrix so conservative residue
    changes can still map even when the full sequences are diverged.
    """
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM45")

    alignments = aligner.align(seq1, seq2)
    best_align = alignments[0]
    aligned_indices_1, aligned_indices_2 = best_align.aligned

    index_map = {}
    for (start1, end1), (start2, end2) in zip(aligned_indices_1, aligned_indices_2):
        for idx1, idx2 in zip(range(start1, end1), range(start2, end2)):
            index_map[idx1 + 1] = idx2 + 1

    return index_map

def compare_interfaces(pairs_1, pairs_2, map_A, map_B):
    """Maps interaction pairs from Species 1 to Species 2 and finds overlap."""
    mapped_pairs_1 = set()
    for res_A, res_B in pairs_1:
        if res_A in map_A and res_B in map_B:
            mapped_pairs_1.add((map_A[res_A], map_B[res_B]))
    
    set_2 = set(pairs_2)
    intersection = mapped_pairs_1.intersection(set_2)
    union = mapped_pairs_1.union(set_2)
    
    return {
        "jaccard_index": len(intersection) / len(union) if union else 0.0,
        "shared_contacts": len(intersection),
        "precision": len(intersection) / len(pairs_2) if pairs_2 else 0.0,
        "recall": len(intersection) / len(pairs_1) if pairs_1 else 0.0
    }


drosophila_seq_dict = map_fragments_to_seqs("/Users/poppy/Dropbox/all_fragments.csv", proteins_to_include=dm_proteins_with_homologs)
human_seq_dict = map_fragments_to_seqs("/Users/poppy/Dropbox/human_screen/alphafragment_human_protein_fragments.csv", proteins_to_include=hs_proteins_with_homologs)

"""
for index, row in interaction_df.iterrows():
    if row["Interaction_Type"] == "Both":
        prot1 = row["Protein1"]
        prot2 = row["Protein2"]
        # for testing, only compare one pair known to be homologous
        #if prot1 != "Ana2" or prot2 != "Sas-4":
        #    continue

        # find all interactions for that drosophila protein pair
        dm_interactions = filtered_drosophila_df[
            ((filtered_drosophila_df["Protein1"] == prot1) & (filtered_drosophila_df["Protein2"] == prot2))
            | ((filtered_drosophila_df["Protein1"] == prot2) & (filtered_drosophila_df["Protein2"] == prot1))
        ]

        # find human interactions for same protein pair
        hs_protA_candidates = homolog_dict.get(prot1, [])
        hs_protB_candidates = homolog_dict.get(prot2, [])
        hs_interactions_list = []
        for hs_protA in hs_protA_candidates:
            for hs_protB in hs_protB_candidates:
                hs_interactions_list.append(filtered_human_df[
                    ((filtered_human_df["Protein1"] == hs_protA) & (filtered_human_df["Protein2"] == hs_protB))
                    | ((filtered_human_df["Protein1"] == hs_protB) & (filtered_human_df["Protein2"] == hs_protA))
                ])
        hs_interactions = pd.concat(hs_interactions_list, ignore_index=True) if hs_interactions_list else pd.DataFrame()
        
        # for each drosophila/human pair, align fragment sequences
        for dm_index, dm_row in dm_interactions.iterrows():
            for hs_index, hs_row in hs_interactions.iterrows():
                protA_dm = dm_row["Protein1"]
                protB_dm = dm_row["Protein2"]
                protA_hs = hs_row["Protein1"]
                protB_hs = hs_row["Protein2"]
                fragA_dm = dm_row["Protein1_Domain"]
                fragB_dm = dm_row["Protein2_Domain"]
                fragA_hs = hs_row["Protein1_Domain"]
                fragB_hs = hs_row["Protein2_Domain"]
                seqA_dm = drosophila_seq_dict.get(fragA_dm)
                seqB_dm = drosophila_seq_dict.get(fragB_dm)
                seqA_hs = human_seq_dict.get(fragA_hs)
                seqB_hs = human_seq_dict.get(fragB_hs)
                if None in [seqA_dm, seqB_dm, seqA_hs, seqB_hs]:
                    print(f"Missing sequence for one of the fragments: {fragA_dm}, {fragB_dm}, {fragA_hs}, {fragB_hs}")
                    continue

                order_flipped = False
                if (protA_hs not in homolog_dict.get(protA_dm, [])) or (protB_hs not in homolog_dict.get(protB_dm, [])):
                    order_flipped = True
                    protA_hs, protB_hs = protB_hs, protA_hs
                    fragA_hs, fragB_hs = fragB_hs, fragA_hs
                    seqA_hs, seqB_hs = seqB_hs, seqA_hs

                # align sequences and get homologous residue maps
                map_A = get_alignment_mapping(seqA_dm, seqA_hs)
                map_B = get_alignment_mapping(seqB_dm, seqB_hs)
                if not map_A or not map_B:
                    print(f"No alignment possible for {fragA_dm} and {fragA_hs} or {fragB_dm} and {fragB_hs}. Skipping comparison.")
                    continue

                # for each potentially homologous protein pair, check if interacting residue pairs match
                pairs_dm = literal_eval(dm_row["residue_pairs"])
                pairs_hs = literal_eval(hs_row["residue_pairs"])
                if order_flipped:
                    pairs_hs = [(resB, resA) for resA, resB in pairs_hs]
                results = compare_interfaces(pairs_hs, pairs_dm, map_A, map_B)
                print(f"Comparison results for interaction between {fragA_dm} and {fragB_dm}:")
                print(f"Jaccard Similarity Score: {results['jaccard_index']:.3f}")
                print(f"Shared Conserved Contacts: {results['shared_contacts']}")
"""

#########################################################
# sequence homology at interaction site
def residue_ids_to_sequence(fragment_sequence, residue_ids):
    """Gets the sequence of residues at the interface"""
    site_residues = []
    for residue_id in sorted(set(residue_ids)):
        sequence_index = residue_id - 1
        if 0 <= sequence_index < len(fragment_sequence):
            site_residues.append(fragment_sequence[sequence_index])
    return "".join(site_residues)

def alignment_homology_score(sequence_dm, sequence_hs, substitution_matrix_name="BLOSUM45"):
    """Return a substitution-matrix alignment score for two sequences."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load(substitution_matrix_name)

    best_alignment = aligner.align(sequence_dm, sequence_hs)[0]
    best_alignment_score = float(best_alignment.score)
    normalised_score = best_alignment_score / max(len(sequence_dm), len(sequence_hs))
    
    return normalised_score


drosophila_seq_dict = map_fragments_to_seqs("/Users/poppy/Dropbox/all_fragments.csv", proteins_to_include=dm_proteins_with_homologs)
human_seq_dict = map_fragments_to_seqs("/Users/poppy/Dropbox/human_screen/alphafragment_human_protein_fragments.csv", proteins_to_include=hs_proteins_with_homologs)

# make df for output
"""sequence_results_db = pd.DataFrame(columns=["dm_prot1", "dm_prot2", "dm_fragA", "dm_fragB", "hs_prot1", "hs_prot2", "hs_fragA", "hs_fragB", "chain_a_alignment_score", "chain_b_alignment_score", "alignment_score"])
for index, row in interaction_df.iterrows():
    if row["Interaction_Type"] == "Both":
        prot1 = row["Protein1"]
        prot2 = row["Protein2"]

        # test with specific protein pair
        #if prot1 != "AurA" or prot2 != "Spd-2":
        #    continue

        # find all interactions in the human and drosophila datasets for that protein pair meeting criteria
        dm_interactions = filtered_drosophila_df[
            ((filtered_drosophila_df["Protein1"] == prot1) & (filtered_drosophila_df["Protein2"] == prot2))
            | ((filtered_drosophila_df["Protein1"] == prot2) & (filtered_drosophila_df["Protein2"] == prot1))
        ]

        hs_prot1_candidates = homolog_dict.get(prot1, [])
        hs_prot2_candidates = homolog_dict.get(prot2, [])
        hs_interactions_list = []
        for hs_prot1 in hs_prot1_candidates:
            for hs_prot2 in hs_prot2_candidates:
                hs_interactions_list.append(filtered_human_df[
                    ((filtered_human_df["Protein1"] == hs_prot1) & (filtered_human_df["Protein2"] == hs_prot2))
                    | ((filtered_human_df["Protein1"] == hs_prot2) & (filtered_human_df["Protein2"] == hs_prot1))
                ])
        hs_interactions = pd.concat(hs_interactions_list, ignore_index=True) if hs_interactions_list else pd.DataFrame()

        # iterate through each homologous pair
        for dm_index, dm_row in dm_interactions.iterrows():
            for hs_index, hs_row in hs_interactions.iterrows():
                order_flipped = False
                protA_dm = dm_row["Protein1"]
                protB_dm = dm_row["Protein2"]
                protA_hs = hs_row["Protein1"]
                protB_hs = hs_row["Protein2"]

                # if human proteins are in opposite order to drosophila, flip human order
                if (protA_hs not in homolog_dict.get(protA_dm, [])) or (protB_hs not in homolog_dict.get(protB_dm, [])):
                    order_flipped = True
                fragA_dm = dm_row["Protein1_Domain"]
                fragB_dm = dm_row["Protein2_Domain"]
                fragA_hs = hs_row["Protein1_Domain"]
                fragB_hs = hs_row["Protein2_Domain"]
                seqA_dm = drosophila_seq_dict.get(fragA_dm)
                seqB_dm = drosophila_seq_dict.get(fragB_dm)
                seqA_hs = human_seq_dict.get(fragA_hs)
                seqB_hs = human_seq_dict.get(fragB_hs)
                if None in [seqA_dm, seqB_dm, seqA_hs, seqB_hs]:
                    print(f"Missing sequence for one of the fragments: {fragA_dm}, {fragB_dm}, {fragA_hs}, {fragB_hs}")
                    continue

                # get interface sequences for each chain in each species
                # adjust b chain numbering to chain relative
                int_seq_dm_a = residue_ids_to_sequence(seqA_dm, dm_row["chain_a_contributors"])
                int_seq_dm_b = residue_ids_to_sequence(seqB_dm, [res - len(seqA_dm) for res in dm_row["chain_b_contributors"]])
                int_seq_hs_a = residue_ids_to_sequence(seqA_hs, hs_row["chain_a_contributors"])
                int_seq_hs_b = residue_ids_to_sequence(seqB_hs, [res - len(seqA_hs) for res in hs_row["chain_b_contributors"]])

                if order_flipped:
                    protA_hs, protB_hs = protB_hs, protA_hs
                    fragA_hs, fragB_hs = fragB_hs, fragA_hs
                    seqA_hs, seqB_hs = seqB_hs, seqA_hs
                    int_seq_hs_a, int_seq_hs_b = int_seq_hs_b, int_seq_hs_a

                score_a = alignment_homology_score(int_seq_dm_a, int_seq_hs_a)
                score_b = alignment_homology_score(int_seq_dm_b, int_seq_hs_b)
                score = (score_a + score_b)/2

                sequence_results_db = pd.concat([sequence_results_db, pd.DataFrame([{
                    "dm_prot1": protA_dm,
                    "dm_prot2": protB_dm,
                    "dm_fragA": fragA_dm,
                    "dm_fragB": fragB_dm,
                    "hs_prot1": protA_hs,
                    "hs_prot2": protB_hs,
                    "hs_fragA": fragA_hs,
                    "hs_fragB": fragB_hs,
                    "chain_a_alignment_score": score_a,
                    "chain_b_alignment_score": score_b,
                    "alignment_score": score,
                }])], ignore_index=True)

# plot histogram of alignment scores
plt.hist(sequence_results_db["alignment_score"].dropna(), bins=20, edgecolor='black')
plt.title("Distribution of Interaction Site Alignment Scores")
plt.xlabel("Average Chain Alignment Score")
plt.ylabel("Frequency")
plt.show()

# make network plot of all > 2.5
high_score_interactions = sequence_results_db[sequence_results_db["alignment_score"] > 2.5]
G = nx.Graph()
for index, row in high_score_interactions.iterrows():
    protA = row["dm_prot1"]
    protB = row["dm_prot2"]
    G.add_edge(protA, protB, weight=row["alignment_score"])
edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True, node_color="lightblue", edge_color='gray', node_size=200, font_size=10)
plt.title("Network of Interactions with High Sequence Homology at Interface")
plt.show()

# combine with structure_results_db, and plot scater plot of sequence homology score vs structural homology score
# condense each dataframe to one row per fragment_pair, taking max score for sequence homology and low score for structural
condensed_sequence_db = sequence_results_db.groupby(["dm_fragA", "dm_fragB", "hs_fragA", "hs_fragB", "dm_prot1", "dm_prot2", "hs_prot1", "hs_prot2"]).agg({
    "alignment_score": "max",
}).reset_index()
condensed_structure_db = results_db.groupby(["dm_fragA", "dm_fragB", "hs_fragA", "hs_fragB", "dm_prot1", "dm_prot2", "hs_prot1", "hs_prot2"]).agg({
    "rmsd": "min",
}).reset_index()

combined_db = pd.merge(condensed_sequence_db, condensed_structure_db, on=["dm_fragA", "dm_fragB", "hs_fragA", "hs_fragB", "dm_prot1", "dm_prot2", "hs_prot1", "hs_prot2"], how="inner")
plt.scatter(combined_db["alignment_score"], combined_db["rmsd"])
plt.title("Sequence Homology Score vs Structural Homology")
plt.xlabel("Sequence Homology Alignment Score")
plt.ylabel("Structural Homology RMSD")
plt.show()

# plot network of all with high sequence homology but low structural homology
high_seq_low_struct = combined_db[(combined_db["alignment_score"] > 2.5) & (combined_db["rmsd"] < 4)]
G = nx.Graph()
for index, row in high_seq_low_struct.iterrows():
    protA = row["dm_prot1"]
    protB = row["dm_prot2"]
    G.add_edge(protA, protB, weight=row["alignment_score"])
edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True, node_color="lightblue", edge_color='gray', node_size=200, font_size=10)
plt.title("Network of Interactions with High Sequence Homology but Low Structural Homology")
plt.show()"""

################################################
# Using structural homology
import numpy as np
from Bio.PDB import PDBParser
import numpy as np
from sklearn.neighbors import NearestNeighbors
import open3d as o3d
import plotly.graph_objects as go

def get_residue_coords_for_ids(pdb_path, residue_ids):
    """Return coordinates for the requested residues.

    `residue_ids` are interpreted as absolute residue positions across the full structure
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)
    model = next(structure.get_models())

    residue_id_set = {int(residue_id) for residue_id in residue_ids}
    coords = []

    global_residue_index = 0
    for model_chain in model:
        for residue in model_chain:
            if residue.id[0] != " ":
                print(f"Skipping non-standard residue {residue} in chain {model_chain.id}")
                continue
            if global_residue_index in residue_id_set:
                if 'CA' in residue:
                    coords.append(residue['CA'].coord)
            global_residue_index += 1

    if not coords:
        return np.empty((0, 3))
    if len(coords) < len(residue_id_set):
        print(f"Warning: Only found coordinates for {len(coords)} out of {len(residue_id_set)} requested residues in {pdb_path}")
    return np.asarray(coords)

def principal_axes_initial_transforms(coords_moving, coords_reference):
    """Computes a set of initial transformations to align two sets of coordinates based on their principal axes."""
    # center both sets of coordinates
    moving_center = np.mean(coords_moving, axis=0)
    reference_center = np.mean(coords_reference, axis=0)
    moving_centered = coords_moving - moving_center
    reference_centered = coords_reference - reference_center

    # compute covariance matrices and their eigen-decomposition to get principal axes
    moving_cov = np.cov(moving_centered, rowvar=False)
    reference_cov = np.cov(reference_centered, rowvar=False)
    moving_eigenvalues, moving_axes = np.linalg.eigh(moving_cov)
    reference_eigenvalues, reference_axes = np.linalg.eigh(reference_cov)

    # sort eigenvalues and corresponding axes in descending order
    moving_order = np.argsort(moving_eigenvalues)[::-1]
    reference_order = np.argsort(reference_eigenvalues)[::-1]
    moving_axes = moving_axes[:, moving_order]
    reference_axes = reference_axes[:, reference_order]

    # ensure right-handed coordinate systems for both sets of axes
    if np.linalg.det(moving_axes) < 0:
        moving_axes[:, -1] *= -1
    if np.linalg.det(reference_axes) < 0:
        reference_axes[:, -1] *= -1

    sign_variants = (
        np.diag([1, 1, 1]),
        np.diag([1, -1, -1]),
        np.diag([-1, 1, -1]),
        np.diag([-1, -1, 1]),
    )

    initial_transforms = []
    for sign_matrix in sign_variants:
        rotation = reference_axes @ sign_matrix @ moving_axes.T
        # skip any transformations that would produce a reflection (negative determinant)
        if np.linalg.det(rotation) < 0:
            continue

        initial_transform = np.eye(4)
        initial_transform[:3, :3] = rotation
        initial_transform[:3, 3] = reference_center - rotation @ moving_center
        initial_transforms.append(initial_transform)

    if not initial_transforms:
        fallback_transform = np.eye(4)
        fallback_transform[:3, 3] = reference_center - moving_center
        initial_transforms.append(fallback_transform)

    return initial_transforms

def align_point_clouds(src, dst, max_iterations=50, tolerance=1e-5):
    """Aligns src (moving coord set) to dst(reference coord set) without needing 1-to-1 mapping."""
    source_pcd = o3d.geometry.PointCloud()
    target_pcd = o3d.geometry.PointCloud()
    source_pcd.points = o3d.utility.Vector3dVector(src)
    target_pcd.points = o3d.utility.Vector3dVector(dst)

    # set max correspondence distance based on the size of the point clouds, with a minimum threshold to allow some flexibility
    all_points = np.vstack([src, dst])
    bbox_diag = float(np.linalg.norm(np.ptp(all_points, axis=0)))
    max_correspondence_distance = max(10.0, 0.25 * bbox_diag) if bbox_diag > 0 else 10.0

    # use point-to-point ICP estimation
    estimation = o3d.pipelines.registration.TransformationEstimationPointToPoint()
    criteria = o3d.pipelines.registration.ICPConvergenceCriteria(
        max_iteration=max_iterations,
        relative_fitness=tolerance,
        relative_rmse=tolerance,
    )

    # run refinement from multiple initial rotations
    best_result = None
    best_key = (-np.inf, np.inf)
    for initial_transform in principal_axes_initial_transforms(src, dst):
        result = o3d.pipelines.registration.registration_icp(
            source_pcd,
            target_pcd,
            max_correspondence_distance,
            initial_transform,
            estimation,
            criteria,
        )
        # assess alignment success based on fitness and inlier RMSE
        result_key = (result.fitness, -result.inlier_rmse)
        if result_key > best_key:
            best_key = result_key
            best_result = result

    if best_result is None:
        return np.copy(src)

    src_homogeneous = np.c_[src, np.ones(len(src))]
    aligned_src = (best_result.transformation @ src_homogeneous.T).T[:, :3]
    return aligned_src

def calculate_point_cloud_similarity(coords_1, coords_2_aligned):
    """
    Calculates the bidirectional nearest-neighbor distance between two aligned protein interface point clouds.
    """
    # Distance from Cloud 1 to nearest in Cloud 2
    neigh1 = NearestNeighbors(n_neighbors=1).fit(coords_2_aligned)
    distances1, _ = neigh1.kneighbors(coords_1)
    
    # Distance from Cloud 2 to nearest in Cloud 1
    neigh2 = NearestNeighbors(n_neighbors=1).fit(coords_1)
    distances2, _ = neigh2.kneighbors(coords_2_aligned)
    
    # Combine them using a root-mean-square approach to mimic RMSD
    combined_squared_distances = np.concatenate([distances1**2, distances2**2])
    return np.sqrt(np.mean(combined_squared_distances))

def plot_interfaces(coords_dm_a, coords_dm_b, coords_hs_a, coords_hs_b,
    name_dm_a="Drosophila - Protein A",
    name_dm_b="Drosophila - Protein B",
    name_hs_a="Human - Protein A",
    name_hs_b="Human - Protein B"
):
    # adjust coordinates so that all axes start at 0
    def adjust_coords(all_coords, coords):
        min_vals = all_coords.min(axis=0)
        return coords - min_vals
    all_coords = np.vstack([coords_dm_a, coords_dm_b, coords_hs_a, coords_hs_b])
    coords_dm_a = adjust_coords(all_coords, coords_dm_a)
    coords_dm_b = adjust_coords(all_coords, coords_dm_b)
    coords_hs_a = adjust_coords(all_coords, coords_hs_a)
    coords_hs_b = adjust_coords(all_coords, coords_hs_b)

    # plot each interface on 3D plot
    fig = go.Figure()
    def add_interface_trace(coords, color, name):
        if len(coords) > 0:
            fig.add_trace(go.Scatter3d(
                x=coords[:, 0],
                y=coords[:, 1],
                z=coords[:, 2],
                mode='markers',
                marker=dict(size=5, color=color),
                name=name
            ))
    add_interface_trace(coords_dm_a, 'blue', name_dm_a)
    add_interface_trace(coords_dm_b, 'orange', name_dm_b)
    add_interface_trace(coords_hs_a, 'green', name_hs_a)
    add_interface_trace(coords_hs_b, 'red', name_hs_b)
    
    fig.update_layout(
        scene=dict(
            aspectmode='data',
        ),
    )
    fig.show()

"""
structure_results_db = pd.DataFrame(columns=["dm_prot1", "dm_prot2", "dm_fragA", "dm_fragB", "hs_prot1", "hs_prot2", "hs_fragA", "hs_fragB", "rmsd", "res_coords_1_a", "res_coords_1_b", "res_coords_2_a", "res_coords_2_b", "human_int_coiledness"])
for index, row in interaction_df.iterrows():
    if row["Interaction_Type"] == "Both":
        prot1 = row["Protein1"]
        prot2 = row["Protein2"]

        # test with single protein pair
        #if prot1 != "Ana2" or prot2 != "Sas-4":
        #    continue

        # get all drosophila interactions for current protein pair
        dm_interactions = filtered_drosophila_df[
            ((filtered_drosophila_df["Protein1"] == prot1) & (filtered_drosophila_df["Protein2"] == prot2))
            | ((filtered_drosophila_df["Protein1"] == prot2) & (filtered_drosophila_df["Protein2"] == prot1))
        ]

        # get all human interactions for homologous protein pairs
        hs_prot1_candidates = homolog_dict.get(prot1, [])
        hs_prot2_candidates = homolog_dict.get(prot2, [])
        hs_interactions_list = []
        for hs_prot1 in hs_prot1_candidates:
            for hs_prot2 in hs_prot2_candidates:
                hs_interactions_list.append(filtered_human_df[
                    ((filtered_human_df["Protein1"] == hs_prot1) & (filtered_human_df["Protein2"] == hs_prot2))
                    | ((filtered_human_df["Protein1"] == hs_prot2) & (filtered_human_df["Protein2"] == hs_prot1))
                ])
        hs_interactions = pd.concat(hs_interactions_list, ignore_index=True) if hs_interactions_list else pd.DataFrame()

        for dm_index, dm_row in dm_interactions.iterrows():
            for hs_index, hs_row in hs_interactions.iterrows():
                try:
                    protA_dm = dm_row["Protein1"]
                    protB_dm = dm_row["Protein2"]
                    protA_hs = hs_row["Protein1"]
                    protB_hs = hs_row["Protein2"]
                    fragA_dm = dm_row["Protein1_Domain"]
                    fragB_dm = dm_row["Protein2_Domain"]
                    fragA_hs = hs_row["Protein1_Domain"]
                    fragB_hs = hs_row["Protein2_Domain"]
                    chain_a_contributors_dm = dm_row["chain_a_contributors"]
                    chain_b_contributors_dm = dm_row["chain_b_contributors"]
                    chain_a_contributors_hs = hs_row["chain_a_contributors"]
                    chain_b_contributors_hs = hs_row["chain_b_contributors"]

                    json_dm = dm_row["log_file"]
                    json_hs = hs_row["log_file"]
                    pdb_dm = json_dm.replace("scores", "unrelaxed").replace(".json", ".pdb")
                    pdb_hs = json_hs.replace("scores", "unrelaxed").replace(".json", ".pdb")

                    order_flipped = False
                    if (protA_hs not in homolog_dict.get(protA_dm, [])) or (protB_hs not in homolog_dict.get(protB_dm, [])):
                        order_flipped = True
                        protA_hs, protB_hs = protB_hs, protA_hs
                        fragA_hs, fragB_hs = fragB_hs, fragA_hs
                        chain_a_contributors_hs, chain_b_contributors_hs = chain_b_contributors_hs, chain_a_contributors_hs

                    coords_dm_a = get_residue_coords_for_ids(pdb_dm, chain_a_contributors_dm)
                    coords_dm_b = get_residue_coords_for_ids(pdb_dm, chain_b_contributors_dm)
                    coords_hs_a = get_residue_coords_for_ids(pdb_hs, chain_a_contributors_hs)
                    coords_hs_b = get_residue_coords_for_ids(pdb_hs, chain_b_contributors_hs)

                    if len(coords_dm_a) < 5 or len(coords_dm_b) < 5 or len(coords_hs_a) < 5 or len(coords_hs_b) < 5:
                        print(f"Warning: Insufficient coordinates for one of the chains in {pdb_dm} or {pdb_hs}. Skipping this pair.")
                        continue
                    
                    coords_dm_combined = np.vstack([coords_dm_a, coords_dm_b])
                    coords_hs_combined = np.vstack([coords_hs_a, coords_hs_b])

                    coords_hs_combined = align_point_clouds(coords_hs_combined, coords_dm_combined, max_iterations=50)

                    chain_a_count_hs = len(coords_hs_a)
                    chain_b_count_hs = len(coords_hs_b)
                    coords_hs_a_aligned = coords_hs_combined[:chain_a_count_hs]
                    coords_hs_b_aligned = coords_hs_combined[chain_a_count_hs:chain_a_count_hs + chain_b_count_hs]

                    similarity_a = calculate_point_cloud_similarity(coords_dm_a, coords_hs_a_aligned)
                    similarity_b = calculate_point_cloud_similarity(coords_dm_b, coords_hs_b_aligned)
                    similarity_score = (similarity_a + similarity_b) / 2.0

                    structure_results_db = pd.concat([structure_results_db, pd.DataFrame([{
                        "dm_prot1": protA_dm,
                        "dm_prot2": protB_dm,
                        "dm_fragA": fragA_dm,
                        "dm_fragB": fragB_dm,
                        "hs_prot1": protA_hs,
                        "hs_prot2": protB_hs,
                        "hs_fragA": fragA_hs,
                        "hs_fragB": fragB_hs,
                        "rmsd": similarity_score,
                        "res_coords_1_a": coords_dm_a,
                        "res_coords_1_b": coords_dm_b,
                        "res_coords_2_a": coords_hs_a_aligned,
                        "res_coords_2_b": coords_hs_b_aligned,
                        "human_int_coiledness": hs_row["coiledness_score"], # include so that coiled coil interactions, which will likely show as homologous, can be filtered out if desired
                    }])], ignore_index=True)

                except Exception as e:
                    print(f"Error processing {fragA_dm}+{fragB_dm} vs {fragA_hs}+{fragB_hs}: {str(e)}")
                    continue

# Overview of results
plt.figure(figsize=(10, 6))
plt.hist(structure_results_db["rmsd"], bins=100, color='green', alpha=0.7, edgecolor='black')
plt.title('Comparison of ICP vs PCA Interface Similarity Scores')
plt.xlabel('Weighted Interface Similarity (Å)')
plt.ylabel('Frequency')
plt.legend()
plt.grid(axis='y')
plt.xlim(0, 10)
plt.show()

# Comparisons below 4.0 Å
below_4 = structure_results_db["rmsd"].lt(4.0).sum()
print(f"\nComparisons with ICP similarity < 4.0 Å: {below_4}")

# plot interfaces in 3D to assess method accuracy
for index, row in structure_results_db.iterrows():
    #if row["rmsd"] < 4.0 and row["human_int_coiledness"] <= 0.25:
    if row["dm_prot1"] == "Spd-2" and row["dm_prot2"] == "AurA" and row["rmsd"] < 4.0:
        print(f"Plotting aligned interfaces for {row['dm_fragA']}+{row['dm_fragB']}/{row['hs_fragA']}+{row['hs_fragB']}")
        print(f"with similarity score {row['rmsd']:.2f} Å")
        plot_interfaces(
            row["res_coords_1_a"],
            row["res_coords_1_b"],
            row["res_coords_2_a"],
            row["res_coords_2_b"],
            name_dm_a=row['dm_prot1'],
            name_dm_b=row['dm_prot2'],
            name_hs_a=row['hs_prot1'],
            name_hs_b=row['hs_prot2'],
        )

# make network plot 
G_icp = nx.Graph()
for index, row in structure_results_db.iterrows():
    if row["rmsd"] <= 4.0 and row["human_int_coiledness"] <= 0.25:
        dm_node1 = f"{row['dm_prot1']}"
        dm_node2 = f"{row['dm_prot2']}"
        G_icp.add_edge(dm_node1, dm_node2, weight=1/row["rmsd"])
plt.figure(figsize=(10, 8))
pos = nx.spring_layout(G_icp, seed=51)
nx.draw(G_icp, pos, with_labels=True, node_color='lightblue', edge_color='gray', node_size=200, font_size=10)
plt.title('Network of Protein Fragments with Similar Interaction Interfaces (ICP Similarity <= 4.0 Å, Coiledness <= 0.25)')
plt.show()
"""

##################################################
# Combine sequence and structural homology scoring
results_db = pd.DataFrame(columns=["dm_prot1", "dm_prot2", "dm_fragA", "dm_fragB", "hs_prot1", "hs_prot2", "hs_fragA", "hs_fragB", "rmsd", "res_coords_1_a", "res_coords_1_b", "res_coords_2_a", "res_coords_2_b", "alignment_score", "human_int_coiledness"])
for index, row in interaction_df.iterrows():
    if row["Interaction_Type"] == "Both":
        prot1 = row["Protein1"]
        prot2 = row["Protein2"]

        # test with single protein pair
        #if prot1 != "Ana2" or prot2 != "Sas-4":
        #    continue

        # get all drosophila interactions for current protein pair
        dm_interactions = filtered_drosophila_df[
            ((filtered_drosophila_df["Protein1"] == prot1) & (filtered_drosophila_df["Protein2"] == prot2))
            | ((filtered_drosophila_df["Protein1"] == prot2) & (filtered_drosophila_df["Protein2"] == prot1))
        ]

        # get all human interactions for homologous protein pairs
        hs_prot1_candidates = homolog_dict.get(prot1, [])
        hs_prot2_candidates = homolog_dict.get(prot2, [])
        hs_interactions_list = []
        for hs_prot1 in hs_prot1_candidates:
            for hs_prot2 in hs_prot2_candidates:
                hs_interactions_list.append(filtered_human_df[
                    ((filtered_human_df["Protein1"] == hs_prot1) & (filtered_human_df["Protein2"] == hs_prot2))
                    | ((filtered_human_df["Protein1"] == hs_prot2) & (filtered_human_df["Protein2"] == hs_prot1))
                ])
        hs_interactions = pd.concat(hs_interactions_list, ignore_index=True) if hs_interactions_list else pd.DataFrame()

        for dm_index, dm_row in dm_interactions.iterrows():
            for hs_index, hs_row in hs_interactions.iterrows():
                try:
                    protA_dm = dm_row["Protein1"]
                    protB_dm = dm_row["Protein2"]
                    protA_hs = hs_row["Protein1"]
                    protB_hs = hs_row["Protein2"]

                    fragA_dm = dm_row["Protein1_Domain"]
                    fragB_dm = dm_row["Protein2_Domain"]
                    fragA_hs = hs_row["Protein1_Domain"]
                    fragB_hs = hs_row["Protein2_Domain"]

                    chain_a_contributors_dm = dm_row["chain_a_contributors"]
                    chain_b_contributors_dm = dm_row["chain_b_contributors"]
                    chain_a_contributors_hs = hs_row["chain_a_contributors"]
                    chain_b_contributors_hs = hs_row["chain_b_contributors"]

                    seqA_dm = drosophila_seq_dict.get(fragA_dm)
                    seqB_dm = drosophila_seq_dict.get(fragB_dm)
                    seqA_hs = human_seq_dict.get(fragA_hs)
                    seqB_hs = human_seq_dict.get(fragB_hs)
                    if None in [seqA_dm, seqB_dm, seqA_hs, seqB_hs]:
                        print(f"Missing sequence for one of the fragments: {fragA_dm}, {fragB_dm}, {fragA_hs}, {fragB_hs}")
                        continue

                    # get interface sequences for each chain in each species
                    # adjust b chain numbering to chain relative
                    int_seq_dm_a = residue_ids_to_sequence(seqA_dm, chain_a_contributors_dm)
                    int_seq_dm_b = residue_ids_to_sequence(seqB_dm, [res - len(seqA_dm) for res in chain_b_contributors_dm])
                    int_seq_hs_a = residue_ids_to_sequence(seqA_hs, chain_a_contributors_hs)
                    int_seq_hs_b = residue_ids_to_sequence(seqB_hs, [res - len(seqA_hs) for res in chain_b_contributors_hs])

                    json_dm = dm_row["log_file"]
                    json_hs = hs_row["log_file"]
                    pdb_dm = json_dm.replace("scores", "unrelaxed").replace(".json", ".pdb")
                    pdb_hs = json_hs.replace("scores", "unrelaxed").replace(".json", ".pdb")

                    order_flipped = False
                    if (protA_hs not in homolog_dict.get(protA_dm, [])) or (protB_hs not in homolog_dict.get(protB_dm, [])):
                        order_flipped = True
                        protA_hs, protB_hs = protB_hs, protA_hs
                        fragA_hs, fragB_hs = fragB_hs, fragA_hs
                        chain_a_contributors_hs, chain_b_contributors_hs = chain_b_contributors_hs, chain_a_contributors_hs
                        seqA_hs, seqB_hs = seqB_hs, seqA_hs
                        int_seq_hs_a, int_seq_hs_b = int_seq_hs_b, int_seq_hs_a
                    
                    score_a = alignment_homology_score(int_seq_dm_a, int_seq_hs_a)
                    score_b = alignment_homology_score(int_seq_dm_b, int_seq_hs_b)
                    score = (score_a + score_b)/2

                    coords_dm_a = get_residue_coords_for_ids(pdb_dm, chain_a_contributors_dm)
                    coords_dm_b = get_residue_coords_for_ids(pdb_dm, chain_b_contributors_dm)
                    coords_hs_a = get_residue_coords_for_ids(pdb_hs, chain_a_contributors_hs)
                    coords_hs_b = get_residue_coords_for_ids(pdb_hs, chain_b_contributors_hs)
                    
                    if len(coords_dm_a) < 5 or len(coords_dm_b) < 5 or len(coords_hs_a) < 5 or len(coords_hs_b) < 5:
                        print(f"Warning: Insufficient coordinates for one of the chains in {pdb_dm} or {pdb_hs}. Skipping this pair.")
                        continue
                    
                    coords_dm_combined = np.vstack([coords_dm_a, coords_dm_b])
                    coords_hs_combined = np.vstack([coords_hs_a, coords_hs_b])

                    coords_hs_combined = align_point_clouds(coords_hs_combined, coords_dm_combined, max_iterations=50)

                    chain_a_count_hs = len(coords_hs_a)
                    chain_b_count_hs = len(coords_hs_b)
                    coords_hs_a_aligned = coords_hs_combined[:chain_a_count_hs]
                    coords_hs_b_aligned = coords_hs_combined[chain_a_count_hs:chain_a_count_hs + chain_b_count_hs]

                    similarity_a = calculate_point_cloud_similarity(coords_dm_a, coords_hs_a_aligned)
                    similarity_b = calculate_point_cloud_similarity(coords_dm_b, coords_hs_b_aligned)
                    similarity_score = (similarity_a + similarity_b) / 2.0

                    results_db = pd.concat([results_db, pd.DataFrame([{
                        "dm_prot1": protA_dm,
                        "dm_prot2": protB_dm,
                        "dm_fragA": fragA_dm,
                        "dm_fragB": fragB_dm,
                        "hs_prot1": protA_hs,
                        "hs_prot2": protB_hs,
                        "hs_fragA": fragA_hs,
                        "hs_fragB": fragB_hs,
                        "rmsd": similarity_score,
                        "res_coords_1_a": coords_dm_a,
                        "res_coords_1_b": coords_dm_b,
                        "res_coords_2_a": coords_hs_a_aligned,
                        "res_coords_2_b": coords_hs_b_aligned,
                        "alignment_score": score,
                        "human_int_coiledness": hs_row["coiledness_score"], # include so that coiled coil interactions, which will likely show as homologous, can be filtered out if desired
                    }])], ignore_index=True)

                except Exception as e:
                    print(f"Error processing {fragA_dm}+{fragB_dm} vs {fragA_hs}+{fragB_hs}: {str(e)}")
                    continue
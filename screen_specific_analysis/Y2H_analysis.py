import pandas as pd
y2h_fragments_csv = "/Users/poppy/Dropbox/Y2H_analysis/Y2H_fragments.csv"
y2h_data_csv = "/Users/poppy/Dropbox/Y2H_analysis/curated_Y2H_fragment_data.csv"
y2h_fragments_df = pd.read_csv(y2h_fragments_csv)
y2h_data_df = pd.read_csv(y2h_data_csv)

# y2h data is a matrix format, where the same fragments are in the top row and leftmost column, and the left/bottom half of the matrix has a score where an interaction has been identified between two fragments (and is blank if no interaction was identified), while the top/right half of the matrix is empty. We want to convert this into a long format, where each row is a pair of fragments and either the interaction score or 0 if no interaction was identified between them.
rows = []

# First column contains row fragment IDs; remaining columns are column fragment IDs
row_fragments = y2h_data_df.iloc[:, 0]
col_fragments = y2h_data_df.columns[1:]

for i, fragment2 in enumerate(row_fragments):
    for j, fragment1 in enumerate(col_fragments):
        # Keep only left/bottom triangle, including the central diagonal
        if j > i:
            continue
        interaction_score = y2h_data_df.iloc[i, j + 1]
        #trim out any whitespace and any '()' (with anything in between) from the interaction score, and convert to int if possible
        if pd.isna(interaction_score):
            interaction_score = 0
        else:
            # use regex to strip '(*)' from the interaction score, and any whitespace
            import re
            interaction_score = re.sub(r"\(.*\)", "", str(interaction_score)).strip()
            # strip any '.0' from the interaction score
            interaction_score = re.sub(r"\.0$", "", interaction_score)
            try:
                interaction_score = int(interaction_score)
            except ValueError:
                print(f"Warning: Could not convert interaction score '{interaction_score}' to int for fragments {fragment1} and {fragment2}. Setting interaction score to 0.")
        rows.append(
            {
                "Fragment1": fragment1,
                "Fragment2": fragment2,
                "Interaction_Score": interaction_score,
            }
        )

print(rows)
y2h_data_long_df = pd.DataFrame(rows)
print(y2h_data_long_df.head())

#now, create a dictionary mapping each fragment to its protein and indices
fragment_mapping = {}
for _, row in y2h_fragments_df.iterrows():
    fragment_id = row["Fragment"]
    protein = row["Protein"]
    start = row["Start"]
    end = row["End"]
    #only add if all values are present
    if pd.notna(fragment_id) and pd.notna(protein) and pd.notna(start) and pd.notna(end):
        fragment_mapping[fragment_id] = {"protein": protein, "start": start, "end": end}
proteins = set()
for fragment_id, info in fragment_mapping.items():
    proteins.add(info["protein"])
print(proteins)

# for each row, add columns for Protein1, Protein2, Protein1_start, Protein1_end, Protein2_start, Protein2_end based on the Y2H fragments df
def get_protein_and_indices(fragment_id):
    if fragment_id in fragment_mapping:
        return (
            fragment_mapping[fragment_id]["protein"],
            fragment_mapping[fragment_id]["start"],
            fragment_mapping[fragment_id]["end"],
        )
    else:
        return None, None, None
y2h_data_long_df[["Protein1", "Fragment1_start", "Fragment1_end"]] = y2h_data_long_df["Fragment1"].apply(
    lambda x: pd.Series(get_protein_and_indices(x))
)
y2h_data_long_df[["Protein2", "Fragment2_start", "Fragment2_end"]] = y2h_data_long_df["Fragment2"].apply(
    lambda x: pd.Series(get_protein_and_indices(x))
)
print(y2h_data_long_df.head())


#read in screen data
"""
af_screen_csv = "/Users/poppy/Desktop/all_interface_analysis_2025.11.26_absolute.csv"
af_screen_df = pd.read_csv(af_screen_csv)
#identify all interaction betwwen the subset of proteins in the Y2H screen, and meeting set criteria
criteria = {
    "min_pae": 11, #must be <=
    "avg_pae": 14, #must be <=
    "rop": 1, #must be >=
    "size": 15, #must be >=
}
filtered_af_screen_df = af_screen_df[
    (af_screen_df["min_pae"] <= criteria["min_pae"])
    & (af_screen_df["avg_pae"] <= criteria["avg_pae"])
    & (af_screen_df["rop"] >= criteria["rop"])
    & (af_screen_df["size"] >= criteria["size"])
    & (af_screen_df["Protein1"].isin(proteins))
    & (af_screen_df["Protein2"].isin(proteins))
]
"""
af_screen_csv = '/Users/Poppy/Dropbox/t7_interface_analysis_coil_absolute_no_dimers_2026.04.29.csv'
af_screen_df = pd.read_csv(af_screen_csv)
#identify all interaction betwwen the subset of proteins in the Y2H screen, and meeting set criteria
criteria = {
    "min_pae": 11, #must be <=
    "avg_pae": 14, #must be <=
    "rop": 1, #must be >=
    "size": 15, #must be >=
    "max_promiscuity": 20, #must be <=
    "pdockq": 0.23, #must be >=
}

filtered_af_screen_df = af_screen_df[
    (af_screen_df["min_pae"] <= criteria["min_pae"])
    & (af_screen_df["avg_pae"] <= criteria["avg_pae"])
    & (af_screen_df["rop"] >= criteria["rop"])
    & (af_screen_df["size"] >= criteria["size"])
    & (af_screen_df["Protein1"].isin(proteins))
    & (af_screen_df["Protein2"].isin(proteins))
    & (af_screen_df["max_promiscuity"] <= criteria["max_promiscuity"])
    & (af_screen_df["pdockq"] >= criteria["pdockq"])
    #& (af_screen_df["chain_a_avg_coil_prob"] <= 0.5)
    #& (af_screen_df["chain_b_avg_coil_prob"] <= 0.5)
]
                                                         
# only keep Protein1, Protein2 and absolute_location columns
filtered_af_screen_df = filtered_af_screen_df[["Protein1", "Protein2", "absolute_location"]]

#absolute location is in format {'Chain A': 'a, c, d-f, p', 'Chain B': 'x, y-z'}.
# Chain A will map to Protein1 and Chain B will map to Protein2.
# We want to convert this into a format where we have Protein1, Protein2, Protein1_start, Protein1_end, Protein2_start, Protein2_end.
# We can do this by parsing the absolute_location column, splitting the list of indices for each chain by any ', ' or '-' and taking the min and max of the resulting list of indices for each chain.
def parse_absolute_location(absolute_location):
    chain_indices = []
    parts = absolute_location.split(", ")
    for part in parts:
        if "-" in part:
            start, end = part.split("-")
            chain_indices.extend(range(int(start), int(end) + 1))
        else:
            chain_indices.append(int(part))
    return min(chain_indices), max(chain_indices)

def get_interaction_indices(absolute_location):
    if pd.isna(absolute_location):
        return None, None, None, None
    try:
        location_dict = eval(absolute_location)
    except Exception as e:
        print(f"Error parsing absolute_location {absolute_location}: {e}")
        return None, None, None, None
    if "Chain A" in location_dict:
        chain_a_location = location_dict["Chain A"]
        chain_a_indices = parse_absolute_location(chain_a_location)
    else:
        chain_a_indices = (None, None)
    if "Chain B" in location_dict:
        chain_b_location = location_dict["Chain B"]
        chain_b_indices = parse_absolute_location(chain_b_location)
    else:
        chain_b_indices = (None, None)

    return chain_a_indices[0], chain_a_indices[1], chain_b_indices[0], chain_b_indices[1]

# apply this function to the absolute_location column to create new columns Protein1_start, Protein1_end, Protein2_start, Protein2_end
filtered_af_screen_df[["Protein1_start", "Protein1_end", "Protein2_start", "Protein2_end"]] = filtered_af_screen_df["absolute_location"].apply(
    lambda x: pd.Series(get_interaction_indices(x))
)

# i want to add a column to the y2h data long df, af_interactions
# this column will be incremented by 1 for each interaction identified in the filtered_af_screen_df between the samepair of proteins as in the y2h data long df
# for each row in the af dataset, find any rows in the y2h dataset where protein1 and protein2 are the same (or protein1=protein2 and protein2=protein1)
# then check if the interaction indices in the af dataset overlap with the fragment indices in the y2h dataset for both proteins. If there is an overlap, increment the af_interactions column for that row in the y2h data long df by 1
y2h_data_long_df["af_interactions"] = 0
for _, af_row in filtered_af_screen_df.iterrows():
    af_protein1 = af_row["Protein1"]
    af_protein2 = af_row["Protein2"]
    af_protein1_start = af_row["Protein1_start"]
    af_protein1_end = af_row["Protein1_end"]
    af_protein2_start = af_row["Protein2_start"]
    af_protein2_end = af_row["Protein2_end"]

    for idx, y2h_row in y2h_data_long_df.iterrows():
        y2h_protein1 = y2h_row["Protein1"]
        y2h_protein2 = y2h_row["Protein2"]
        y2h_fragment1_start = y2h_row["Fragment1_start"]
        y2h_fragment1_end = y2h_row["Fragment1_end"]
        y2h_fragment2_start = y2h_row["Fragment2_start"]
        y2h_fragment2_end = y2h_row["Fragment2_end"]

        # check if proteins match (in either order)
        if (af_protein1 == y2h_protein1 and af_protein2 == y2h_protein2):
            # check for overlap in interaction indices and fragment indices for both proteins
            if (
                max(af_protein1_start, y2h_fragment1_start) <= min(af_protein1_end, y2h_fragment1_end)
                and max(af_protein2_start, y2h_fragment2_start) <= min(af_protein2_end, y2h_fragment2_end)
            ):
                # increment the af_interactions column for this row in the y2h data long df
                y2h_data_long_df.at[idx, "af_interactions"] += 1
        elif (af_protein1 == y2h_protein2 and af_protein2 == y2h_protein1):
            # check for overlap in interaction indices and fragment indices for both proteins (swapping protein1 and protein2)
            if (
                max(af_protein1_start, y2h_fragment2_start) <= min(af_protein1_end, y2h_fragment2_end)
                and max(af_protein2_start, y2h_fragment1_start) <= min(af_protein2_end, y2h_fragment1_end)
            ):
                # increment the af_interactions column for this row in the y2h data long df
                y2h_data_long_df.at[idx, "af_interactions"] += 1

print(y2h_data_long_df.head())

# add another column, status, which is "None" if af_interactions and y2h_interactions are 0, "Y2H" if af_interactions is 0, "AF" if y2h_interactions is 0, and "Both" if both af_interactions and y2h_interactions are >0
def determine_status(row):
    af = int(row["af_interactions"])
    y2h = int(row["Interaction_Score"])
    if af == 0 and y2h == 0:
        return "None"
    elif af == 0 and y2h > 0:
        return "Y2H"
    elif af > 0 and y2h == 0:
        return "AF"
    else:
        return "Both"
y2h_data_long_df["status"] = y2h_data_long_df.apply(determine_status, axis=1)
print(y2h_data_long_df.head())

# save this dataframe to a csv on the desktop
y2h_data_long_df.to_csv("/Users/poppy/Desktop/Y2H_AF_comparison.csv", index=False)

# plot a bar chart showing the count of each status category
import matplotlib.pyplot as plt
status_counts = y2h_data_long_df["status"].value_counts()
plt.bar(status_counts.index, status_counts.values)
plt.xlabel("Status")
plt.ylabel("Count")
plt.title("Comparison of Y2H and AF interactions")
plt.show()

print(status_counts)

y2h_data_long_df_filtered = y2h_data_long_df[y2h_data_long_df["status"] != "None"]
status_counts_filtered = y2h_data_long_df_filtered["status"].value_counts()
plt.bar(status_counts_filtered.index, status_counts_filtered.values)
plt.xlabel("Status")
plt.ylabel("Count")
plt.title("Comparison of Y2H and AF interactions")
plt.show()

# stratify bars by y2h interaction score, showing a bar for each status category with bands for each interaction score
y2h_data_long_df_filtered["Interaction_Score"] = y2h_data_long_df_filtered["Interaction_Score"].astype(int)
interaction_score_bins = [-1, 0, 1, 2, 3]
y2h_data_long_df_filtered["Interaction_Score_Bin"] = pd.cut(y2h_data_long_df_filtered["Interaction_Score"], bins=interaction_score_bins, labels=["0", "1", "2", "3"])
interaction_score_status_counts = y2h_data_long_df_filtered.groupby(["status", "Interaction_Score_Bin"]).size().unstack(fill_value=0)
interaction_score_status_counts.plot(kind="bar", stacked=True)
plt.xlabel("Status")
plt.ylabel("Count")
plt.title("Comparison of Y2H and AF interactions stratified by Y2H interaction score")
plt.legend(title="Y2H Interaction Score")
plt.show()



#collapse to one row per protein pair, where status = "None" if all rows for that protein pair have status "None", "Y2H" if any row for that protein pair has status "Y2H" and no rows have status "AF" or "Both", "AF" if any row for that protein pair has status "AF" and no rows have status "Y2H" or "Both", "Both" if any row for that protein pair has status "Both", or mixed if there are rows with both "Y2H" and "AF" and no rows with "Both"
def determine_overall_status(group):
    statuses = set(group["status"])
    if statuses == {"None"}:
        return "None"
    elif "Both" in statuses:
        return "Both"
    elif "Y2H" in statuses and "AF" in statuses:
        return "Mixed"
    elif "Y2H" in statuses:
        return "Y2H"
    elif "AF" in statuses:
        return "AF"
protein_pair_status_df = y2h_data_long_df.groupby(["Protein1", "Protein2"]).apply(determine_overall_status).reset_index()
protein_pair_status_df.columns = ["Protein1", "Protein2", "Overall_Status"]
print(protein_pair_status_df.head())
#rename to work in gephi
protein_pair_status_df = protein_pair_status_df.rename(columns={"Protein1": "Source", "Protein2": "Target", "Overall_Status": "Status"})
protein_pair_status_df["Weight"] = 1
print(protein_pair_status_df.head())

# make network plot
import networkx as nx
import matplotlib.patches as mpatches
G = nx.from_pandas_edgelist(protein_pair_status_df[protein_pair_status_df["Status"] != "None"], "Source", "Target", ["Status", "Weight"])
color_map = {"Y2H": "blue", "AF": "orange", "Both": "green", "Mixed": "red"}
edge_colors = [color_map[G.edges[edge]["Status"]] for edge in G.edges()]
pos = nx.spring_layout(G)
plt.figure(figsize=(10, 10))
nx.draw(G, pos, with_labels=True, edge_color=edge_colors, node_color="lightgrey", node_size=500)
plt.title("Comparison of Y2H and AF interactions")
blue_patch = mpatches.Patch(color='blue', label='Y2H')
purple_patch = mpatches.Patch(color='orange', label='AF')
green_patch = mpatches.Patch(color='green', label='Both')
red_patch = mpatches.Patch(color='red', label='Mixed')
plt.legend(handles=[blue_patch, purple_patch, green_patch, red_patch])
plt.show()

import datetime
date = datetime.datetime.now().strftime("%Y-%m-%d")
protein_pair_status_df.to_csv(f"/Users/poppy/Desktop/Y2H_AF_comparison_protein_pairs_{date}.csv", index=False)

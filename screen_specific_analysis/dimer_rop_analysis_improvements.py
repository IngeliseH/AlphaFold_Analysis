
import pandas as pd
df_old = pd.read_csv('/Users/poppy/Dropbox/all_dimer_interface_analysis_2025.06.05.csv')
df_new = pd.read_csv('/Users/poppy/Dropbox/all_dimer_interface_analysis_2025.11.24.csv')
# compare num of rows (as some reruns added in between)
print(f"Old df rows: {len(df_old)}, New df rows: {len(df_new)}") # 9968 vs 9995 - 27 more rows in new df
# compare rop distribution
import matplotlib.pyplot as plt
plt.hist(df_old['rop'], bins=20, alpha=0.5, label='Old DF')
plt.hist(df_new['rop'], bins=20, alpha=0.5, label='New DF')
# count num of rop 0, 1, 2, 3, 4 in old and in new
for i in range(5):
    old_count = len(df_old[df_old['rop'] == i])
    new_count = len(df_new[df_new['rop'] == i])
    print(f"Rop {i}: Old count = {old_count}, New count = {new_count}.")
# Rop 0: Old count = 8950, New count = 8679. Difference = -271
# Rop 1: Old count = 739, New count = 828. Difference = 89
# Rop 2: Old count = 184, New count = 266. Difference = 82
# Rop 3: Old count = 75, New count = 147. Difference = 72
# Rop 4: Old count = 20, New count = 75. Difference = 55
# compare rop of rows directly (using Protein1_Domain and Protein2_Domain columns as id) and make df of rop changes
df_changes = pd.DataFrame(columns=['Protein1_Domain', 'Protein2_Domain', 'Location', 'Old_rop', 'New_rop', 'rop_change'])
for _, row in df_new.iterrows():
    p1d = row['Protein1_Domain']
    p2d = row['Protein2_Domain']
    location = row['location']
    new_rop = row['rop']
    old_row = df_old[(df_old['Protein1_Domain'] == p1d) & (df_old['Protein2_Domain'] == p2d) & (df_old['location'] == location)]
    if not old_row.empty:
        old_rop = old_row.iloc[0]['rop']
        new_rop = row['rop']
        rop_change = new_rop - old_rop
        df_changes = pd.concat([df_changes, pd.DataFrame({'Protein1_Domain': [p1d], 'Protein2_Domain': [p2d], 'Location': [location], 'Old_rop': [old_rop], 'New_rop': [new_rop], 'rop_change': [rop_change]})], ignore_index=True)
print(f"Number of rows with changed rop: {len(df_changes['rop_change'][df_changes['rop_change'] != 0])}")
print(df_changes)
# plot rop changes (all values will be integers between -4 and +4)
plt.figure()
plt.hist(df_changes['rop_change'], bins=range(-4, 6), align='left', rwidth=0.8)
# print count of each rop change value
for i in range(-4, 5):
    count = len(df_changes[df_changes['rop_change'] == i])
    print(f"Rop change {i}: Count = {count}")
# Rop change -4: Count = 0
# Rop change -3: Count = 0
# Rop change -2: Count = 0
# Rop change -1: Count = 10
# Rop change 0: Count = 7638
# Rop change 1: Count = 400
# Rop change 2: Count = 68
# Rop change 3: Count = 13
# Rop change 4: Count = 0
# surprising that there are any decreases at all - need to investigate these cases individually
for _, row in df_changes[df_changes['rop_change'] < 0].iterrows():
    print(row['Protein1_Domain'], row['Protein2_Domain'], row['Old_rop'], row['New_rop'])
# Asl_dimer_F2 FZY_F2 1 0 # same model number
# Asl_dimer_F2 FZY_F2 1 0
# Ana1_F5 Asl_dimer_F2 2 1
# Ana1_F5 Asl_dimer_F2 2 1
# PLK4_dimer_F1 PLP_F1 1 0
# Msps_F2 TACC_dimer_F2 1 0
# MPS1_F2 TACC_dimer_F2 1 0
# MPS1_F2 TACC_dimer_F2 1 0
# GCP5_F2 Sas-6_dimer_F1 1 0
# GCP5_F2 Sas-6_dimer_F1 1 0

# compare dimeric to monomeric predictions - look for changes in pae and rop
# make plot of pae changes and plot of rop changes - for cahnges need to pick one interface per domain pair- and may not be the same for old and new so not direct comparision
# could instead just plot all pae and rop values for monomeric vs dimeric predictions - but not as informative?
# pick best rop for each pair / best min pae for each pair / best avg pae for each pair
# for each, plot best dimeric vs best monomeric value on x and y axis
# start by making a df of best values for each pair for monomeric and dimeric predictions
monomeric_df = pd.read_csv('/Users/poppy/Dropbox/all_interface_analysis_2025.11.24.csv')
dimeric_df = pd.read_csv('/Users/poppy/Dropbox/all_dimer_interface_analysis_2025.11.24.csv')
# make summarised version of monomeric df with best scores for each domain pair
monomeric_summary = pd.DataFrame(columns=['Protein1_Domain', 'Protein2_Domain', 'Best_monomeric_rop', 'Best_monomeric_min_pae', 'Best_monomeric_avg_pae'])
for _, row in monomeric_df.iterrows():
    p1d = row['Protein1_Domain']
    p2d = row['Protein2_Domain']
    rop = row['rop']
    min_pae = row['min_pae']
    avg_pae = row['avg_pae']
    existing_row = monomeric_summary[((monomeric_summary['Protein1_Domain'] == p1d) & (monomeric_summary['Protein2_Domain'] == p2d)) | ((monomeric_summary['Protein1_Domain'] == p2d) & (monomeric_summary['Protein2_Domain'] == p1d))]
    if existing_row.empty:
        monomeric_summary = pd.concat([monomeric_summary, pd.DataFrame({'Protein1_Domain': [p1d], 'Protein2_Domain': [p2d], 'Best_monomeric_rop': [rop], 'Best_monomeric_min_pae': [min_pae], 'Best_monomeric_avg_pae': [avg_pae]})], ignore_index=True)
    else:
        idx = existing_row.index[0]
        if rop > monomeric_summary.at[idx, 'Best_monomeric_rop']:
            monomeric_summary.at[idx, 'Best_monomeric_rop'] = rop
        if min_pae < monomeric_summary.at[idx, 'Best_monomeric_min_pae']:
            monomeric_summary.at[idx, 'Best_monomeric_min_pae'] = min_pae
        if avg_pae < monomeric_summary.at[idx, 'Best_monomeric_avg_pae']:
            monomeric_summary.at[idx, 'Best_monomeric_avg_pae'] = avg_pae
# now for each dimeric row, find corresponding monomeric row and get best scores
dimeric_summary = pd.DataFrame(columns=['Protein1_Domain', 'Protein2_Domain', 'Best_dimeric_rop', 'Best_dimeric_min_pae', 'Best_dimeric_avg_pae'])
for _, row in dimeric_df.iterrows():
    if row['Protein1_Domain'].endswith('_dimer') & row['Protein2_Domain'].endswith('_dimer'):
        continue
    p1d = row['Protein1_Domain'].replace('_dimer', '')
    p2d = row['Protein2_Domain'].replace('_dimer', '')
    rop = row['rop']
    min_pae = row['min_pae']
    avg_pae = row['avg_pae']
    existing_row = dimeric_summary[((dimeric_summary['Protein1_Domain'] == p1d) & (dimeric_summary['Protein2_Domain'] == p2d)) | ((dimeric_summary['Protein1_Domain'] == p2d) & (dimeric_summary['Protein2_Domain'] == p1d))]
    if existing_row.empty:
        dimeric_summary = pd.concat([dimeric_summary, pd.DataFrame({'Protein1_Domain': [p1d], 'Protein2_Domain': [p2d], 'Best_dimeric_rop': [rop], 'Best_dimeric_min_pae': [min_pae], 'Best_dimeric_avg_pae': [avg_pae]})], ignore_index=True)
    else:
        idx = existing_row.index[0]
        if rop > dimeric_summary.at[idx, 'Best_dimeric_rop']:
            dimeric_summary.at[idx, 'Best_dimeric_rop'] = rop
        if min_pae < dimeric_summary.at[idx, 'Best_dimeric_min_pae']:
            dimeric_summary.at[idx, 'Best_dimeric_min_pae'] = min_pae
        if avg_pae < dimeric_summary.at[idx, 'Best_dimeric_avg_pae']:
            dimeric_summary.at[idx, 'Best_dimeric_avg_pae'] = avg_pae
# now merge the two summaries on Protein1_Domain and Protein2_Domain, skipping those that don't have both monomeric and dimeric entries
merged_summary = pd.merge(monomeric_summary, dimeric_summary, on=['Protein1_Domain', 'Protein2_Domain'])
print(merged_summary)
# plot best dimeric vs best monomeric rop
for metric in ['min_pae', 'avg_pae']:
    plt.figure()
    plt.scatter(merged_summary[f'Best_monomeric_{metric}'], merged_summary[f'Best_dimeric_{metric}'], alpha=0.2)
    plt.xlabel(f'Best Monomeric {metric.capitalize()}')
    plt.ylabel(f'Best Dimeric {metric.capitalize()}')
    plt.title(f'Best Dimeric vs Monomeric {metric.capitalize()}')
# rop is discrete so needs different plotting
plt.figure()
plt.hist2d(merged_summary['Best_monomeric_rop'], merged_summary['Best_dimeric_rop'], bins=[range(0,5), range(0,5)], cmap='RdBu_r')
plt.xlabel('Best Monomeric Rop')
plt.ylabel('Best Dimeric Rop')
plt.title('Best Dimeric vs Monomeric Rop')
plt.colorbar(label='Counts')
plt.show()

# print rows where dimeric is better than monomeric
print("Rows where dimeric min_pae is better than monomeric:")
for _, row in merged_summary.iterrows():
    if (row['Best_dimeric_min_pae'] < 5) and (row['Best_monomeric_min_pae'] > 5):
        print(row['Protein1_Domain'], row['Protein2_Domain'], row['Best_dimeric_min_pae'], row['Best_monomeric_min_pae'])
print("Rows where dimeric rop is better than monomeric:")
for _, row in merged_summary.iterrows():
    if row['Best_dimeric_rop'] - row['Best_monomeric_rop'] >= 2:
        print(row['Protein1_Domain'], row['Protein2_Domain'], row['Best_dimeric_rop'], row['Best_monomeric_rop'])

# compare incidence of predictions with controls in dimeric vs equivalent monomeric predictions
# list control proteins
controls = ['FUS', 'cpcB', 'bamC', 'PIK3CA', 'ERD10', 'MFP1']
dimers = ['Asl', 'PLK4', 'Sas-6', 'TACC', 'Cnn']
# make df of only rows involving controls, and these proteins as either dimers or monomers
dimer_monomer_control_df = pd.DataFrame()
for control in controls:
    for protein in dimers:
        dimeric_rows = dimeric_df[((dimeric_df['Protein1'] == f"{protein}_dimer") & (dimeric_df['Protein2'] == control)) | ((dimeric_df['Protein1'] == control) & (dimeric_df['Protein2'] == f"{protein}_dimer"))]
        monomeric_rows = monomeric_df[((monomeric_df['Protein1'] == protein) & (monomeric_df['Protein2'] == control)) | ((monomeric_df['Protein1'] == control) & (monomeric_df['Protein2'] == protein))]
        # add column to indicate dimeric or monomeric
        dimeric_rows = dimeric_rows.copy()
        dimeric_rows['Prediction_Type'] = 'Dimeric'
        monomeric_rows = monomeric_rows.copy()
        monomeric_rows['Prediction_Type'] = 'Monomeric'
        dimer_monomer_control_df = pd.concat([dimer_monomer_control_df, dimeric_rows, monomeric_rows], ignore_index=True)
print(f"Total rows involving controls: {len(dimer_monomer_control_df)}")
print(dimer_monomer_control_df)
# plot bar chart of rop distribution for dimeric vs monomeric predictions involving controls
import numpy as np
rop_values = [0, 1, 2, 3, 4]
dimeric_counts = [len(dimer_monomer_control_df[(dimer_monomer_control_df['Prediction_Type'] == 'Dimeric') & (dimer_monomer_control_df['rop'] == rop)]) for rop in rop_values]
monomeric_counts = [len(dimer_monomer_control_df[(dimer_monomer_control_df['Prediction_Type'] == 'Monomeric') & (dimer_monomer_control_df['rop'] == rop)]) for rop in rop_values]
x = np.arange(len(rop_values))
width = 0.35
plt.figure()
plt.bar(x - width/2, dimeric_counts, width, label='Dimeric')
plt.bar(x + width/2, monomeric_counts, width, label='Monomeric')
plt.xlabel('ROP Value')
plt.ylabel('Counts')
plt.title('ROP Distribution for Dimeric vs Monomeric Predictions Involving Controls')
plt.xticks(x, rop_values)
plt.legend()
plt.show()

# repeat for min_pae (non discrete so use histograms)
plt.figure()
plt.hist(dimer_monomer_control_df[dimer_monomer_control_df['Prediction_Type'] == 'Dimeric']['min_pae'], bins=30, alpha=0.5, label='Dimeric')
plt.hist(dimer_monomer_control_df[dimer_monomer_control_df['Prediction_Type'] == 'Monomeric']['min_pae'], bins=30, alpha=0.5, label='Monomeric')
plt.xlabel('Min PAE')
plt.ylabel('Counts')
plt.title('Min PAE Distribution for Dimeric vs Monomeric Predictions Involving Controls')
plt.legend()
plt.show()

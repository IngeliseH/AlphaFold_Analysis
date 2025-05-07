##################################################################################################################################
# COMPARING DIMER AND MONOMER PREDICTIONS
# most dimer predictions seem to be same as monomer predictions (only actively involving one chain of each protein)
# most interested in those that are not - so looking for domain pairs where a stronger prediction is made when a dimer is included
##################################################################################################################################
import pandas as pd
# read in dimer and monomer dataframes
monomer_csv = "/Users/poppy/Dropbox/all_interface_analysis_2025.04.29.csv"
centriole_dimer_csv = "/Users/poppy/Dropbox/Centriole_dimers/Centriole_dimers_interface_analysis.csv"
pcm_dimer_csv = "/Users/poppy/Dropbox/PCM_dimers/PCM_dimers_interface_analysis.csv"
other_dimer_csv = "/Users/poppy/Dropbox/more_dimers/more_dimers_interface_analysis.csv"
# read in dimer dataframes
monomers = pd.read_csv(monomer_csv)
centriole_dimers = pd.read_csv(centriole_dimer_csv)
pcm_dimers = pd.read_csv(pcm_dimer_csv)
other_dimers = pd.read_csv(other_dimer_csv)
# combine dimer dataframes
dimers = pd.concat([centriole_dimers, pcm_dimers, other_dimers], ignore_index=True)

# condense dataframes to best stats for each domain pair (Protein1_Domain and Protein2_Domain)
def condense_df(df, col1, col2):
    # define whether to keep max or min values for each column
    min = ['min_pae', 'avg_pae']
    max = ['rop', 'iptm', 'pdockq', 'evenness', 'avg_pct_rop']
    # create a new dataframe to store the condensed data, with col1, col2, and columns from min and max lists if present in df
    condensed_df = pd.DataFrame(columns=[col1, col2] + [col for col in df.columns if col in min + max])
    # iterate through each unique pair of col1 and col2
    for protein1, protein2 in zip(df[col1], df[col2]):
        # check if the pair is already in condensed_df
        if not ((condensed_df[col1] == protein1) & (condensed_df[col2] == protein2)).any():
            # if not,  find all instances of the pair in df
            df_pair = df[(df[col1] == protein1) & (df[col2] == protein2)]
            # make a new row for the condensed dataframe
            new_row = {col1: protein1, col2: protein2}
            # for each column in df_pair, check if it is min or max (ignore if not in min or max)
            for col in df_pair.columns:
                if col in min:
                    # if it is min, take the min value from df_pair
                    new_row[col] = df_pair[col].min()
                elif col in max:
                    # if it is max, take the max value from df_pair
                    new_row[col] = df_pair[col].max()
            # append the new row to the condensed dataframe
            condensed_df = pd.concat([condensed_df, pd.DataFrame([new_row])], ignore_index=True)
    # return the condensed dataframe
    return condensed_df

# condense monomer and dimer dataframes
monomers = condense_df(monomers, 'Protein1_Domain', 'Protein2_Domain')
dimers = condense_df(dimers, 'Protein1_Domain', 'Protein2_Domain')
                    
# rename all columns of dimer dataframe to dimer_name
dimers = dimers.rename(columns=lambda x: x + '_dimer')
# add columns 'Protein1_Domain', 'Protein2_Domain' to dimer dataframe for domain names with '_dimer' excluded, to match monomers
dimers['Protein1_Domain'] = dimers['Protein1_Domain_dimer'].str.replace('_dimer', '')
dimers['Protein2_Domain'] = dimers['Protein2_Domain_dimer'].str.replace('_dimer', '')
# print first few rows of the dimer dataframe to check the columns
print(dimers.head())
# merge dataframes by domain pair names
merged_df = pd.merge(monomers, dimers, on=['Protein1_Domain', 'Protein2_Domain'], how='inner')
# identify all rows where min_pae_dimer is less than min_pae and less than 5
dimer_improvement = merged_df[(merged_df['min_pae_dimer'] < merged_df['min_pae']) & (merged_df['min_pae_dimer'] < 5)]
# print protein1 and protein2 domain names, min_pae_dimer, and min_pae
print(dimer_improvement[['Protein1_Domain_dimer', 'Protein2_Domain_dimer', 'min_pae_dimer', 'min_pae']])


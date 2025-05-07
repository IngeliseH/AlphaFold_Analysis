import requests
import csv
import pandas as pd

################################################
# Sourcing interactions from BioGRID
################################################

# Get interactions from BioGRID
accessKey = '17b87d81ccfddcd5de129aee48d49106'
# Define list of proteins to search for
geneList = ['alphaTub67C', 'ana1', 'ana2', 'ana3', 'asl', 'betaTub56D',
            'Cnb', 'Cep135', 'Cep97', 'Cp110', 'ctp', 'SAK',
            'Rcd4', 'Sas-4', 'Sas-6', 'asp', 'aur', 'Chc', 'cnn',
            'msps', 'mud', 'polo', 'cp309', 'spd-2', 'tacc', 'gammaTub37C',
            'Grip84', 'Grip91', 'Grip75', 'Grip128', 'Grip163', 'Grip71',
            'CG42787', 'cdc2', 'cdc2c', 'CycA', 'CycB', 'CycB3', 'CycE', 'ial',
            'Bub1', 'BubR1', 'Cenp-C', 'cid', 'fzy', 'Mad1', 'mad2', 'ald',
            'rod', 'Zw10', 'Zwilch']
searchNames=True # search by gene name
geneListJoined = '|'.join(geneList)
taxId = '7227' # Drosophila melanogaster

request_url = "https://webservice.thebiogrid.org/interactions/?taxid=" + taxId + "&searchNames=true&geneList=" + geneListJoined + "&accesskey=" + accessKey
r = requests.get(request_url)
interactions = r.text

# Format interactions into dataframe
lines = interactions.split('\n')
interactions_df = pd.DataFrame(columns=['Protein1AF', 'Protein2AF', 'Gene1bg', 'Gene2BG', 'Evidence_count', 'Evidence'])
galletta_interactions_df = pd.DataFrame(columns=['Protein1AF', 'Protein2AF', 'Gene1BG', 'Gene2BG', 'Evidence_count'])

# Dict for conversion of biogrid naming format to naming used in my AF screen
bg_af = {'ald': 'MPS1', 'alphaTub67C': 'alpha-tubulin', 'ana1': 'Ana1', 'ana2': 'Ana2',
            'ana3': 'Ana3', 'asl': 'Asl', 'asp': 'Asp', 'aur': 'AurA', 'betaTub56D': 'beta-tubulin',
            'Bub1': 'BUB1', 'BubR1': 'BUBR1', 'cdc2': 'Cdk1', 'cdc2c': 'Cdk2',
            'Cenp-C': 'CENPC', 'CG42787': 'Mzt1', 'Chc': 'CHC', 'cid': 'CID',
            'Cnb': 'Centrobin', 'cnn': 'Cnn', 'Cp110': 'CP110', 'cp309': 'PLP',
            'ctp': 'Ctp', 'CycB': 'CycB1', 'fzy': 'FZY', 'gammaTub37C': 'gamma-tubulin',
            'Grip128': 'GCP5', 'Grip163': 'Gcp6', 'Grip75': 'GCP4', 'Grip84': 'GCP2',
            'Grip91': 'GCP3', 'ial': 'AurB', 'Mad1': 'MAD1', 'mad2': 'MAD2',
            'msps': 'Msps', 'mud': 'Mud', 'polo': 'PLK1', 'rod': 'ROD', 'SAK': 'PLK4',
            'spd-2': 'Spd-2', 'tacc': 'TACC', 'Zw10': 'ZW10', 'Zwilch': 'ZWILCH',
            'Cep135': 'Cep135', 'Cep97': 'Cep97', 'CycA': 'CycA', 'CycB3': 'CycB3',
            'CycE': 'CycE', 'Grip71': 'Grip71', 'Rcd4': 'Rcd4', 'Sas-4': 'Sas-4',
            'Sas-6': 'Sas-6'}

# iterate through lines and add to dataframe
for line in lines:
    fields = line.split('\t')
    if len(fields) > 1:
        gene1bg = fields[7]
        gene2bg = fields[8]
        # only include lines where both proteins are of interest
        if gene1bg in geneList and gene2bg in geneList:
            # convert gene names to AF naming format
            gene1af = bg_af[gene1bg]
            gene2af = bg_af[gene2bg]
            # sort gene names alphabetically by AF naming format
            if gene1af.lower() > gene2af.lower():
                gene1af, gene2af =  gene2af, gene1af
                gene1bg, gene2bg = gene2bg, gene1bg
            # add data to dataframe
            if fields[12] == 'physical': # only including evidence of physical interactions
                # join fields 11 (method) + 13 (ref) into one string
                evidence = fields[11] + ' ' + fields[13]
                interactions_df = pd.concat([interactions_df, pd.DataFrame([[gene1af, gene2af, gene1bg, gene2bg, 1, evidence]], columns=['Protein1AF', 'Protein2AF', 'Gene1BG', 'Gene2BG', 'Evidence_count', 'Evidence'])])
                if fields[13] == 'Galletta BJ (2016)':
                    galletta_interactions_df = pd.concat([galletta_interactions_df, pd.DataFrame([[gene1af, gene2af, gene1bg, gene2bg, 1]], columns=['Protein1AF', 'Protein2AF', 'Gene1BG', 'Gene2BG', 'Evidence_count'])])

# Merge duplicate interactions, summing evidence count
interactions_df = interactions_df.groupby(['Protein1AF', 'Protein2AF', 'Gene1BG', 'Gene2BG'], as_index=False).agg({'Evidence_count': 'sum', 'Evidence': lambda x: ', '.join(x)})
galletta_interactions_df = galletta_interactions_df.groupby(['Protein1AF', 'Protein2AF', 'Gene1BG', 'Gene2BG'], as_index=False).agg({'Evidence_count': 'sum'})

# Write all interactions to csv
bg_data_filename = 'BioGRID_interactions.csv'
with open(f'/Users/poppy/Dropbox/{bg_data_filename}', 'w') as f:
    writer = csv.writer(f)
    # write all columns from pd dataframe
    writer.writerow(['Protein1AF', 'Protein2AF', 'Gene1BG', 'Gene2BG', 'Evidence_count', 'Evidence'])
    for index, row in interactions_df.iterrows():
        writer.writerow([row['Protein1AF'], row['Protein2AF'], row['Gene1BG'], row['Gene2BG'], row['Evidence_count'], row['Evidence']])

# Write Galletta interactions to csv
Galletta_data_filename = 'Galletta2016_interactions.csv'
with open(f'/Users/poppy/Dropbox/{Galletta_data_filename}', 'w') as f:
    writer = csv.writer(f)
    # write all columns from pd dataframe
    writer.writerow(['Protein1AF', 'Protein2AF', 'Gene1BG', 'Gene2BG', 'Evidence_count'])
    for index, row in galletta_interactions_df.iterrows():
        writer.writerow([row['Protein1AF'], row['Protein2AF'], row['Gene1BG'], row['Gene2BG'], row['Evidence_count']])


############################################
# Comparing BioGRID to AlphaFold
############################################

# Read in AlphaFold data
af = pd.read_csv('/Users/poppy/Dropbox/gephi_input_strict_filtered_2025.03.19.csv')
af2 = pd.read_csv('/Users/poppy/Dropbox/gephi_input_2025.04.04.csv')

# Function to combine two dataframes to gephi format - one with BioGRID data and one with AlphaFold data
def compare_interactions(bg_df, af_df, save_name=None):
    """
    Compare BioGRID and AlphaFold interaction dataframes.

    Parameters:
        bg_df (pd.DataFrame): BioGRID dataframe with columns ['Protein1AF', 'Protein2AF', 'Evidence_count'].
        af_df (pd.DataFrame): AlphaFold dataframe with columns ['Source', 'Target', 'Weight'].
        save_name (str): Optional name to save the comparison dataframe.

    Returns:
        pd.DataFrame: Comparison dataframe with columns ['Source', 'Target', 'Weight', 'AF_ROP', 'seen_in', 'AF', 'BG'].
    """
    comparison_df = pd.DataFrame(columns=['Source', 'Target', 'Weight', 'AF_ROP', 'seen_in', 'AF', 'BG'])

    # Iterate through BioGRID interactions and add to comparison_df
    for index, row in bg_df.iterrows():
        source = row['Protein1AF']
        target = row['Protein2AF']
        weight = row['Evidence_count']
        af_rop = 0
        af = False
        bg = True
        seen_in = 'bg'

        # Check if interaction is in AF data - if so, combine info from both
        af_match = af_df[(af_df['Source'] == source) & (af_df['Target'] == target)]
        if not af_match.empty:
            af_rop = af_match['Weight'].values[0]
            af = True
            seen_in = 'both'
        comparison_df = pd.concat([comparison_df, pd.DataFrame([[source, target, weight, af_rop, seen_in, af, bg]], 
                                                               columns=['Source', 'Target', 'Weight', 'AF_ROP', 'seen_in', 'AF', 'BG'])])

    # Iterate through AF interactions and add to comparison_df if not already present
    for index, row in af_df.iterrows():
        source = row['Source']
        target = row['Target']
        if comparison_df[(comparison_df['Source'] == source) & (comparison_df['Target'] == target)].empty:
            weight = 0.01
            af_rop = row['Weight']
            af = True
            bg = False
            seen_in = 'af'
            comparison_df = pd.concat([comparison_df, pd.DataFrame([[source, target, weight, af_rop, seen_in, af, bg]], 
                                                                   columns=['Source', 'Target', 'Weight', 'AF_ROP', 'seen_in', 'AF', 'BG'])])
            
    if save_name:
        if save_name.endswith('.csv'):
            save_path = f'/Users/poppy/Dropbox/{save_name}'
        else:
            save_path = f'/Users/poppy/Dropbox/{save_name}.csv'
        with open(save_path, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['Source', 'Target', 'Weight', 'AF_ROP', 'seen_in', 'AF', 'BG'])
            for index, row in comparison_df.iterrows():
                writer.writerow([row['Source'], row['Target'], row['Weight'], row['AF_ROP'], row['seen_in'], row['AF'], row['BG']])

    return comparison_df

bg_af2_comparison = compare_interactions(interactions_df, af2, save_name='BioGRID_AF_2025.04.04_comparison.csv')

# Make df for comparison of all interactions with AF, and second for Galletta interactions with AF
bg_af_comparison = compare_interactions(interactions_df, af)
galletta_af_comparison = compare_interactions(galletta_interactions_df, af)

bg_data_filename = 'BioGRID_AF_comparison.csv'
with open(f'/Users/poppy/Dropbox/{bg_data_filename}', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['Source', 'Target', 'Weight', 'AF_ROP', 'seen_in', 'AF', 'BG'])
    for index, row in bg_af_comparison.iterrows():
        writer.writerow([row['Source'], row['Target'], row['Weight'], row['AF_ROP'], row['seen_in'], row['AF'], row['BG']])

galletta_data_filename = 'Galletta2016_AF_comparison.csv'
with open(f'/Users/poppy/Dropbox/{galletta_data_filename}', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['Source', 'Target', 'Weight', 'AF_ROP', 'seen_in', 'AF', 'BG'])
    for index, row in galletta_af_comparison.iterrows():
        writer.writerow([row['Source'], row['Target'], row['Weight'], row['AF_ROP'], row['seen_in'], row['AF'], row['BG']])
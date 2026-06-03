import pandas as pd
from ast import literal_eval

def prep_coconat_input(protein_data_csv, output_file):
    protein_data = pd.read_csv(protein_data_csv)
    fragment_dict = {}
    with open(output_file, 'w') as f:
        for index, row in protein_data.iterrows():
            if 'dimer' in row['name'].lower():
                continue
            name = row['name']
            fragment_sequences = literal_eval(row['fragment_sequences'])
            for i, seq in enumerate(fragment_sequences):
                fragment_dict[f'{name}_F{i+1}'] = seq
                if len(seq) > 700:  # sequences >700 res not allowed in coconat
                    parts = []
                    remaining = seq
                    part_num = 1
                    while len(remaining) > 700:
                        parts.append(remaining[:600])
                        remaining = remaining[600:]
                        part_num += 1
                    parts.append(remaining)  # add the final part
                    
                    for part_idx, part in enumerate(parts, 1):
                        f.write(f'>{name}_F{i+1}_part{part_idx}\n')
                        f.write(f'{part}\n')
                else:
                    f.write(f'>{name}_F{i+1}\n')
                    f.write(f'{seq}\n')
    # count number of sequences in output file
    with open(output_file, 'r') as f:    num_sequences = sum(1 for line in f if line.startswith('>'))
    print(f'Number of sequences in output file: {num_sequences}')
    return fragment_dict

drosophila_fragment_dict = prep_coconat_input('/Users/poppy/Dropbox/all_fragments.csv', '/Users/poppy/Desktop/coconat_input.fasta')
human_fragment_dict = prep_coconat_input('/Users/poppy/Dropbox/alphafragment_human_protein_fragments.csv', '/Users/poppy/Desktop/human_coconat_input.fasta')

def prep_interface_data(interface_data_csv):
    interface_data = pd.read_csv(interface_data_csv)
    # remove dimers as haven't implemented analysis for these
    interface_data = interface_data[~interface_data['Protein1'].str.contains('dimer', case=False) & ~interface_data['Protein2'].str.contains('dimer', case=False)]
    return interface_data

def find_interface_contributors(interface_data, fragment_dict):
    # get interface contributors from each side for each interface
    interface_data['residue_pairs'] = interface_data['residue_pairs'].apply(literal_eval)
    interface_data['chain_a_contributors'] = interface_data['residue_pairs'].apply(lambda pairs: set(pair[0] for pair in pairs))
    interface_data['chain_b_contributors'] = interface_data['residue_pairs'].apply(lambda pairs: set(pair[1] for pair in pairs))

    # adjust chain b to account for pdb wide residue numbering
    def adjust_chain_b_contributors(row):
        chain_a_length = len(fragment_dict[row['Protein1_Domain']])
        adjusted_contributors = set(res - chain_a_length for res in row['chain_b_contributors'])
        return adjusted_contributors
    interface_data['chain_b_contributors'] = interface_data.apply(adjust_chain_b_contributors, axis=1)
    return interface_data

# collapse interface data to one row per fragment pair, keeping pdockq and list of interface contributors from each chain
def collapse_interface_data(interface_data):
    collapsed_data = interface_data.groupby(['Protein1_Domain', 'Protein2_Domain']).agg({
        'chain_a_contributors': lambda x: set().union(*x),
        'chain_b_contributors': lambda x: set().union(*x),
        'pdockq': 'first',
        'min_pae': 'min',
        'avg_pae': 'min',
        'rop': 'max',
        'iptm': 'first',
        'evenness': 'max',
        'max_promiscuity': 'min'
    }).reset_index()

    return collapsed_data

#interface_data = prep_interface_data('/Users/Poppy/Dropbox/t7_interface_analysis_2026.04.11.csv')
interface_data = prep_interface_data('/Users/Poppy/Desktop/t7_interface_analysis_2026.04.29_absolute.csv')
interface_data = find_interface_contributors(interface_data, drosophila_fragment_dict)
collapsed_data = collapse_interface_data(interface_data)

human_interface_data = prep_interface_data('/Users/Poppy/Desktop/human_interface_analysis_2026.04.30.csv')
human_interface_data = find_interface_contributors(human_interface_data, human_fragment_dict)
human_collapsed_data = collapse_interface_data(human_interface_data)

def process_coil_prob_data(coil_data_json):
    coil_data = pd.read_json(coil_data_json)
    coil_data['fragment_id'] = coil_data['accession']
    coil_data['sequence'] = coil_data['res']
    coil_data['length'] = coil_data['length']
    coil_data['coil_probs'] = coil_data['prob']
    coil_df = coil_data[['fragment_id', 'sequence', 'length', 'coil_probs']]

    # combine fragments that had to be split into multiple parts for coil analysis
    combined_indices = []
    for index, row in coil_df.iterrows():
        if '_part1' in row['fragment_id']:
            base_id = row['fragment_id'].replace('_part1', '')
            combined_sequence = row['sequence']
            combined_coil_probs = row['coil_probs']
            combined_length = row['length']
            part_num = 2
            while True:
                part_id = f'{base_id}_part{part_num}'
                part_row = coil_df[coil_df['fragment_id'] == part_id]
                if part_row.empty:
                    break
                combined_sequence += part_row['sequence'].values[0]
                combined_coil_probs += part_row['coil_probs'].values[0]
                combined_length += part_row['length'].values[0]
                combined_indices.append(part_row.index[0])
                part_num += 1
            # Update the _part1 row with combined data
            coil_df.at[index, 'sequence'] = combined_sequence
            coil_df.at[index, 'coil_probs'] = combined_coil_probs
            coil_df.at[index, 'length'] = combined_length
            coil_df.at[index, 'fragment_id'] = base_id
    # Remove all part rows (including the updated _part1)
    coil_df = coil_df.drop(combined_indices)

    return coil_df

def validate_coil_data(coil_df, fragment_dict):
    # check that every fragment from fragment_dict is present in coil_df
    missing_fragments = set(fragment_dict.keys()) - set(coil_df['fragment_id'])
    if missing_fragments:
        print(f'Missing fragments in coil data: {missing_fragments}')
    # check there are no duplicate fragment ids in coil_df
    if coil_df['fragment_id'].duplicated().any():
        print('Duplicate fragment ids found in coil data:')
        print(coil_df[coil_df['fragment_id'].duplicated(keep=False)]['fragment_id'])
    # check there are no fragment ids in coil_df that are not in fragment_dict, and if there are, print which ones
    extra_fragments = set(coil_df['fragment_id']) - set(fragment_dict.keys())
    if extra_fragments:
        print(f'Extra fragments in coil data that are not in fragment_dict: {extra_fragments}')
    # check sequence lengths match between coil_df and fragment_dict
    length_mismatches = []
    for _, row in coil_df.iterrows():
        frag_id = row['fragment_id']
        if frag_id in fragment_dict:
            coil_len = len(row['sequence'])
            frag_len = len(fragment_dict[frag_id])
            if coil_len != frag_len:
                length_mismatches.append((frag_id, coil_len, frag_len))
    if length_mismatches:
        for frag_id, coil_len, frag_len in length_mismatches:
                print(f'Length mismatch for {frag_id}: coil_df length={coil_len}, fragment_dict length={frag_len}')

# for each row in the collapsed interface data, find the corresponding coil probabilities for the chain a and chain b contributors, and average these to get an average coil propensity for the interface contributors in chain a and chain b, and add these as new columns to the collapsed interface data
# first add list of propensities for each chain a and chain b contributor
def add_coil_prob_data(interface_data, coil_df):
    def get_chain_probs(row, column, coil_df):
        coil_match = coil_df.loc[coil_df['fragment_id'] == row['Protein1_Domain'], 'coil_probs']
        if coil_match.empty:
            return []
        probs = coil_match.values[0]
        chain_probs = []
        for res in row[column]:
            if res < len(probs):
                chain_probs.append(probs[res])
        return chain_probs
    
    interface_data['chain_a_coil_probs'] = interface_data.apply(lambda row: get_chain_probs(row, 'chain_a_contributors', coil_df), axis=1)
    interface_data['chain_b_coil_probs'] = interface_data.apply(lambda row: get_chain_probs(row, 'chain_b_contributors', coil_df), axis=1)

    # get averages
    interface_data['chain_a_avg_coil_prob'] = interface_data['chain_a_coil_probs'].apply(lambda probs: sum(probs) / len(probs) if len(probs) > 0 else 0)
    interface_data['chain_b_avg_coil_prob'] = interface_data['chain_b_coil_probs'].apply(lambda probs: sum(probs) / len(probs) if len(probs) > 0 else 0)
    return interface_data

coil_df = process_coil_prob_data('/Users/Poppy/Desktop/coconat_results_reformatted.json')
validate_coil_data(coil_df, drosophila_fragment_dict)
collapsed_data = add_coil_prob_data(collapsed_data, coil_df)
interface_data = add_coil_prob_data(interface_data, coil_df)  # add to non collapsed interface data as well

human_coil_df = process_coil_prob_data('/Users/Poppy/Desktop/human_coconat_results_reformatted.json')
validate_coil_data(human_coil_df, human_fragment_dict)
human_collapsed_data = add_coil_prob_data(human_collapsed_data, human_coil_df)
human_interface_data = add_coil_prob_data(human_interface_data, human_coil_df)

# save non collapsed interface data with coil propensities as csv
interface_data.to_csv('/Users/Poppy/Dropbox/t7_interface_analysis_coil_absolute_no_dimers_2026.04.29.csv', index=False)

#############################################################
#ANALYSIS AND PLOTTING
#############################################################
# plot spread of pdockq scores alone
import matplotlib.pyplot as plt
plt.hist(collapsed_data['pdockq'], bins=20, color='blue', alpha=0.7)
plt.xlim(0, 1)
plt.xlabel('pdockq')
plt.ylabel('Frequency')
plt.title('Distribution of pdockq Scores')
plt.show()

# filter out interface size < 15
filtered_interface_data = interface_data[interface_data['size'] >= 15]

# Function to plot scatter of chain A vs B coil propensity colored by a metric
def plot_coil_scatter_by_metric(data, metric, cmap='viridis', clim=(0, 1)):
    plt.scatter(data['chain_a_avg_coil_prob'], data['chain_b_avg_coil_prob'], c=data[metric], cmap=cmap, alpha=0.3)
    plt.colorbar(label=metric)
    if clim is not None:
        plt.clim(*clim)
    plt.xlabel('Chain A Average Coil Propensity')
    plt.ylabel('Chain B Average Coil Propensity')
    plt.title(f'Average Coil Propensity of Interface Contributors vs {metric}')
    plt.show()

# Plot for each metric
plot_coil_scatter_by_metric(collapsed_data, 'pdockq')
plot_coil_scatter_by_metric(collapsed_data, 'iptm')
plot_coil_scatter_by_metric(filtered_interface_data, 'min_pae', cmap='viridis_r', clim=None)  # Use full data
plot_coil_scatter_by_metric(filtered_interface_data, 'avg_pae', cmap='viridis_r', clim=None)  # Use full data
plot_coil_scatter_by_metric(filtered_interface_data, 'rop', clim=None)  # Use full data
plot_coil_scatter_by_metric(interface_data, 'max_promiscuity', clim=None)  # Use full data

#human
human_filtered_interface_data = human_interface_data[human_interface_data['size'] >= 15]
plot_coil_scatter_by_metric(human_collapsed_data, 'pdockq')
plot_coil_scatter_by_metric(human_collapsed_data, 'iptm')
plot_coil_scatter_by_metric(human_filtered_interface_data, 'min_pae', cmap='viridis_r', clim=None)  # Use full data
plot_coil_scatter_by_metric(human_filtered_interface_data, 'avg_pae', cmap='viridis_r', clim=None)  # Use full data
plot_coil_scatter_by_metric(human_filtered_interface_data, 'rop', clim=None)  # Use full data
plot_coil_scatter_by_metric(human_filtered_interface_data, 'max_promiscuity', clim=None)  # Use full data

# split into 3 interface types - high coil propensity in both chains, high coil propensity in one chain only, low coil propensity in both chains
def classify_interface(row, threshold=0.5):
    if row['chain_a_avg_coil_prob'] > threshold and row['chain_b_avg_coil_prob'] > threshold:
        return 'coil + coil'
    elif row['chain_a_avg_coil_prob'] > threshold or row['chain_b_avg_coil_prob'] > threshold:
        return 'coil + non-coil'
    else:
        return 'non-coil + non-coil'

collapsed_data['interface_type'] = collapsed_data.apply(classify_interface, axis=1)
interface_data['interface_type'] = interface_data.apply(classify_interface, axis=1)

human_collapsed_data['interface_type'] = human_collapsed_data.apply(classify_interface, axis=1)
human_filtered_interface_data['interface_type'] = human_filtered_interface_data.apply(classify_interface, axis=1)

# Function to plot violin plot of metric vs interface type
def plot_metric_vs_type(data, metric, ylim=(0, 1), plot_type='violin'):
    import seaborn as sns
    if plot_type == 'violin':
        sns.violinplot(x='interface_type', y=metric, data=data)
    else:
        sns.boxplot(x='interface_type', y=metric, data=data)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel('Interface Type')
    plt.ylabel(metric)
    plt.title(f'{metric} Distribution by Interface Type')
    plt.show()

# Plot for each metric
for metric in ['pdockq', 'iptm']:
    plot_metric_vs_type(collapsed_data, metric, plot_type='violin')
    plot_metric_vs_type(collapsed_data, metric, plot_type='box')
for metric in ['min_pae', 'avg_pae', 'rop']:
    plot_metric_vs_type(collapsed_data, metric, ylim=None, plot_type='violin')
    plot_metric_vs_type(collapsed_data, metric, ylim=None, plot_type='box')

# human
for metric in ['pdockq', 'iptm']:
    plot_metric_vs_type(human_collapsed_data, metric, plot_type='violin')
    plot_metric_vs_type(human_collapsed_data, metric, plot_type='box')
for metric in ['min_pae', 'avg_pae', 'rop']:
    plot_metric_vs_type(human_filtered_interface_data, metric, ylim=None, plot_type='violin')
    plot_metric_vs_type(human_filtered_interface_data, metric, ylim=None, plot_type='box')

# statistical test to see if there is a significant difference in pdockq between the 3 interface types
from scipy.stats import kruskal, mannwhitneyu

def perform_statistical_test(data, value_col, group_col='interface_type'):
    groups = data[group_col].unique()
    group_values = [data[data[group_col] == group][value_col] for group in groups]
    stat, p_value = kruskal(*group_values)
    print(f'For category {value_col}: Kruskal-Wallis H-statistic: {stat}, p-value: {p_value}')

    # pairwise comparisons if the Kruskal-Wallis test is significant
    if p_value < 0.05:
        print(f'    Pairwise comparisons for {value_col}:')
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                group1 = data[data[group_col] == groups[i]][value_col]
                group2 = data[data[group_col] == groups[j]][value_col]
                stat, p_value = mannwhitneyu(group1, group2)
                print(f'    {groups[i]} vs {groups[j]}: U-statistic: {stat}, p-value: {p_value}')
                if p_value < 0.05:
                    print(f'        The difference between {groups[i]} and {groups[j]} is statistically significant.')
                else:
                    print(f'        The difference between {groups[i]} and {groups[j]} is not statistically significant.')

# Perform for each metric
for metric in ['pdockq', 'iptm', 'min_pae', 'avg_pae', 'rop']:
    perform_statistical_test(collapsed_data, metric)

for metric in ['pdockq', 'iptm', 'min_pae', 'avg_pae', 'rop']:
    perform_statistical_test(human_collapsed_data, metric)

# want to see how good pdockq is as a classifer for coiled interfaces
from scipy.stats import spearmanr

# Create combined coiledness score: product of both chains' average coil probs (captures asymmetry)
def calculate_coiledness_score(a, b):
    if a + b == 0:
        return 0
    return a * b * (1 - abs(a - b) / (a + b))

collapsed_data['coiledness_score'] = collapsed_data.apply(lambda row: calculate_coiledness_score(row['chain_a_avg_coil_prob'], row['chain_b_avg_coil_prob']), axis=1)
# add to full data and save
interface_data['coiledness_score'] = interface_data.apply(lambda row: calculate_coiledness_score(row['chain_a_avg_coil_prob'], row['chain_b_avg_coil_prob']), axis=1)
# save without interface contributor columns
to_save = interface_data.drop(columns=['chain_a_contributors', 'chain_b_contributors', 'chain_a_coil_probs', 'chain_b_coil_probs'])
# drop any columns with no header
to_save = to_save.drop(columns=[col for col in to_save.columns if not col.strip()])
# change protein1 to Source, protein2 to Target, and add weight column with value 1
to_save = to_save.rename(columns={'Protein1': 'Source', 'Protein2': 'Target'})
to_save['weight'] = 1
to_save.to_csv('/Users/Poppy/Dropbox/t7_interface_analysis_with_coil_2026.04.29.csv', index=False)

# Drop rows with NA in coiledness_score or metrics
metrics = ['pdockq', 'iptm', 'min_pae', 'avg_pae', 'rop']
eval_data = collapsed_data.dropna(subset=['coiledness_score'] + metrics)

human_interface_data['coiledness_score'] = human_interface_data.apply(lambda row: calculate_coiledness_score(row['chain_a_avg_coil_prob'], row['chain_b_avg_coil_prob']), axis=1)
human_to_save = human_interface_data.drop(columns=['chain_a_contributors', 'chain_b_contributors', 'chain_a_coil_probs', 'chain_b_coil_probs'])
human_to_save = human_to_save.drop(columns=[col for col in human_to_save.columns if not col.strip()])
human_to_save = human_to_save.rename(columns={'Protein1': 'Source', 'Protein2': 'Target'})
human_to_save['weight'] = 1
human_to_save.to_csv('/Users/Poppy/Dropbox/human_interface_analysis_with_coil_2026.04.30.csv', index=False)
human_collapsed_data['coiledness_score'] = human_collapsed_data.apply(lambda row: calculate_coiledness_score(row['chain_a_avg_coil_prob'], row['chain_b_avg_coil_prob']), axis=1)

# Compute and print Spearman correlation for each metric vs coiledness_score
for metric in metrics:
    corr, p_value = spearmanr(eval_data['coiledness_score'], eval_data[metric])
    print(f'Spearman correlation between {metric} and coiledness score: {corr:.3f} (p-value: {p_value:.3e})')
    # interpret correlation
    if abs(corr) >= 0.7:
        interpretation = 'strong'
    elif abs(corr) >= 0.5:
        interpretation = 'moderate'
    elif abs(corr) >= 0.3:
        interpretation = 'weak'
    else:
        interpretation = 'negligible'
    print(f'    Interpretation: {interpretation} correlation')


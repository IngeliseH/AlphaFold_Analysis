import os
import pandas as pd
from tqdm import tqdm
import traceback
from iptm_visualisation import extract_iptm
from analysis_utility import find_rank_001_files, extract_pae, parse_structure_file
from process_pae import find_min_pae

def find_proteins(folder):
    """
    Extract protein names from the folder name.
    Expects the folder name to be in the format "fold_protein1_fl_protein2_fl", "fold_protein1_protein2_fl" or "fold_protein_fl_dimer".
    """
    # Get the folder name and split by underscores
    base_folder = os.path.basename(folder)
    parts = base_folder.split('_')
    if parts[-1] == "dimer":    # format = "fold_protein_fl_dimer"
        protein1, protein2 = parts[1], parts[1]
    elif len(parts) == 5:   # format = "fold_protein1_fl_protein2_fl"
        protein1, protein2 = parts[1], parts[3]
    elif len(parts) == 4 or len(parts) == 3:    # format = "fold_protein1_protein2_fl" or "fold_protein1_protein2"
        protein1, protein2 = parts[1], parts[2]
    else:
        raise ValueError(f"Unexpected folder name format: {base_folder}")
    return protein1, protein2

### RUN ANALYSIS ON FL SCREEN ###
def process_all_predictions(base_folder):
    """
    Process all AlphaFold predictions in the given base folder and return results as a DataFrame.
    Expects base folder to contain protein pair folders.

    Parameters:
        - base_folder (str): Path to the base folder containing all protein pair folders.
    
    Returns:
        - DataFrame containing protein pairs with their ipTM and min_PAE values
    """
    all_data = []
    protein_pair_folders = [d for d in os.listdir(base_folder) if os.path.isdir(os.path.join(base_folder, d))]

    with tqdm(total=len(protein_pair_folders), desc='Processing Protein Pairs') as pbar:
        # Find each protein pair folder
        for protein_pair_folder in protein_pair_folders:
            protein_pair_path = os.path.join(base_folder, protein_pair_folder)
            protein1, protein2 = find_proteins(protein_pair_folder)
            if protein2 < protein1:
                protein1, protein2 = protein2, protein1   # alphabetically order protein names
            
            # find ipTM and minPAE for each protein pair
            try:
                pbar.set_postfix_str(f"{protein1} - {protein2}")
                structure_file, json_file, log_file, _, _ = find_rank_001_files(protein_pair_path, af3=True)
                iptm = extract_iptm(log_file)
                pae = extract_pae(json_file)
                model = parse_structure_file(structure_file, is_pdb=False)
                min_pae = find_min_pae(model, pae)[0]
                
                # Append the results to the list
                all_data.append({
                    'Protein1': protein1,
                    'Protein2': protein2,
                    'fl_iptm': iptm,
                    'fl_min_pae': min_pae
                })
                    
            except Exception as e:
                print(f"Error processing protein pair {protein_pair_folder}: {e}")
                traceback.print_exc()
            pbar.update(1)

    if all_data:
        return pd.DataFrame(all_data)
    else:
        print("No data found.")

directory = "/Users/poppy/Dropbox/AF3_FL_Centriole_screen"
data = process_all_predictions(directory)
print(data[:5])

### READ IN FRAGMENT DATA AND REFORMAT TO MATCH FL DATA ###
def reformat_interface_data(input_path):
    """
    Read in the interface data and reformat it to match the full length data.
    
    Parameters:
        - input_path (str): Path to the CSV file containing the interface data.
    
    Returns:
        - DataFrame containing the reformatted interface data.
    """
    # read in csv
    fragment_data = pd.read_csv(input_path)
    def rename_protein(old_name, new_name):
        fragment_data['Protein1'] = fragment_data['Protein1'].str.replace(old_name, new_name)
        fragment_data['Protein2'] = fragment_data['Protein2'].str.replace(old_name, new_name)
    # rename protein names to match full length data
    rename_protein('Sas-4', 'sas4')
    rename_protein('Sas-6', 'sas6')
    rename_protein('beta-tubulin', 'btub')
    rename_protein('alpha-tubulin', 'atub')
    # for each row, lower case and alphabetically order the source and target columns
    fragment_data['Protein1'] = fragment_data['Protein1'].str.lower()
    fragment_data['Protein2'] = fragment_data['Protein2'].str.lower()
    for index, row in fragment_data.iterrows():
        source = row['Protein1']
        target = row['Protein2']
        if source > target:
            fragment_data.at[index, 'Protein1'] = target
            fragment_data.at[index, 'Protein2'] = source
    return fragment_data

fragment_input = '/Users/poppy/Dropbox/all_interface_analysis_2025.02.26.csv'
fragment_data = reformat_interface_data(fragment_input)
print(fragment_data[:5])

### CONDENSE FRAGMENT DATA TO ONE VALUE PER PROTEIN PAIR ###
def condense_to_one_val_per_protein_pair(dataframe, column_name, func='min'):
    """
    Condense the dataframe to one value per protein pair for the specified column using the specified function.
    
    Parameters:
        - dataframe (DataFrame): The input DataFrame.
        - column_name (str): The name of the column to condense.
        - func (str): The function to use for condensing ('min' or 'max').
    """
    condensed_data = dataframe[['Protein1', 'Protein2', column_name]].copy()
    # rename column to 'fragments_' + column_name
    condensed_data.rename(columns={column_name: 'fragments_' + column_name}, inplace=True)
    # find min or max for each protein pair
    if func == 'min':
        condensed_data = condensed_data.groupby(['Protein1', 'Protein2']).min().reset_index()
    elif func == 'max':
        condensed_data = condensed_data.groupby(['Protein1', 'Protein2']).max().reset_index() 
    else:
        raise ValueError("Prioritise must be either 'min' or 'max'")
    return condensed_data
# rop
fragment_rop = condense_to_one_val_per_protein_pair(fragment_data, 'rop', func='max')
# min_pae
fragment_min_pae = condense_to_one_val_per_protein_pair(fragment_data, 'min_pae', func='min')
# iptm
fragment_iptm = condense_to_one_val_per_protein_pair(fragment_data, 'iptm', func='max')
# merge all three dataframes - highest rop and iptm, lowest min pae for each protein pair (not necessarily all from same fragment pair)
fragment_scores = pd.merge(fragment_rop, fragment_min_pae, on=["Protein1", "Protein2"], how="outer")
fragment_scores = pd.merge(fragment_scores, fragment_iptm, on=["Protein1", "Protein2"], how="outer")

### MERGE FRAGMENT DATA AND FL DATA ###
full_data = pd.merge(data, fragment_scores, on=["Protein1", "Protein2"], how="inner") # using inner to exclude any protein pairs not run with FL
print(full_data[:5])

### ADD INPUT LENGTH COLUMN ###
# dictionary of protein lengths
lengths = {'rcd4': 199, 'poc1': 391, 'ana2': 420, 'btub': 447, 'atub': 462, 'sas6': 472,
           'centrobin': 689, 'plk4': 769, 'cep97': 806, 'sas4': 901, 'asl': 994,
           'cep135': 1059, 'ana1': 1729, 'ana3': 1977, 'plp': 2895, 'cdk1': 297,
           'cdk2': 314, 'aura': 411, 'cyca': 491, 'cycb1': 530, 'cycb3': 575, 'plk1': 576,
           'cyce': 709, 'mzt1': 82, 'gtub': 457, 'grip71': 646, 'gcp4': 650, 'gcp2': 852,
           'gcp3': 917, 'gcp5': 1092, 'spd2': 1146, 'cnn': 1148, 'tacc': 1226, 'gcp6': 1351,
           'chc': 1678, 'asp': 1954, 'msps': 2042, 'mud': 2567, 'cp110': 666, 'centrin': 182,
           'ctp': 89, 'cid': 225, 'aurb': 329, 'fzy': 526, 'bubr1': 1460, 'cenpc': 1411,
           'mad1': 730, 'mps1': 672, 'mad2': 207, 'bub1': 1099, 'zw10': 721, 'zwilch': 641,
           'rod': 2098}
# for each row in full_data, sum lengths of protein1 and protein2 to give input length
full_data['input_length'] = 0
full_data['min_length'] = 0
full_data['max_length'] = 0
for index, row in full_data.iterrows():
    length1 = lengths[row['Protein1']]
    length2 = lengths[row['Protein2']]
    full_data.at[index, 'input_length'] = length1 + length2
    full_data.at[index, 'min_length'] = min(length1, length2)
    full_data.at[index, 'max_length'] = max(length1, length2)
print(full_data[:5])

### MAKE SCATTERPLOTS ###
import matplotlib.pyplot as plt
import seaborn as sns
def make_scatterplot(data, x, y, colour, size, title=None, xlabel=None, ylabel=None, save_path=None):
    """
    Create a scatterplot with the given parameters.
    
    Parameters:
        - data (DataFrame): The input DataFrame.
        - x (str): The name of the column for the x-axis.
        - y (str): The name of the column for the y-axis.
        - colour (str): The name of the column for the colour.
        - size (str): The name of the column for the size of the points.
        - title (str): The title of the plot.
        - xlabel (str): The label for the x-axis.
        - ylabel (str): The label for the y-axis.
        - save_path (str): The path to save the plot.
    """
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=data, x=x, y=y, hue=colour, palette='viridis', size=size, sizes=(20, 200))
    if not title:
        title = f"{y}"
    if not xlabel:
        xlabel = x
    if not ylabel:
        ylabel = y
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if save_path:
        plt.savefig(save_path)
    plt.show()
# fl_min_pae vs fragments_min_pae
make_scatterplot(full_data, 'fragments_min_pae', 'fl_min_pae', colour='input_length', size='fragments_rop')
# fl_iptm vs fragments_iptm
make_scatterplot(full_data, 'fragments_ipTM', 'fl_ipTM', colour='input_length', size='fragments_rop')

### MAKE INTERACTIVE PLOTS ###
import plotly.express as px
def make_interactive_plot(data, x, y, colour, size, save_path, title=None, x_label=None, y_label=None):
    """
    Create an interactive scatterplot with the given parameters.
    
    Parameters:
        - data (DataFrame): The input DataFrame.
        - x (str): The name of the column for the x-axis.
        - y (str): The name of the column for the y-axis.
        - colour (str): The name of the column for the colour.
        - size (str): The name of the column for the size of the points.
        - title (str): The title of the plot.
        - save_path (str): The path to save the plot.
    """
    import plotly.express as px
    fig = px.scatter(data, x=x, y=y, hover_name='Protein1', hover_data=['Protein2'], color=colour, size=size)
    fig.update_traces(marker=dict(line=dict(width=2, color='DarkSlateGrey')), selector=dict(mode='markers'))
    if not title:
        title = f"{y} vs {x}"
    if not x_label:
        x_label = x
    if not y_label:
        y_label = y
    fig.update_layout(title=title, xaxis_title=x_label, yaxis_title=y_label)
    if save_path:
        fig.write_html(save_path)
    fig.show()

# fl_min_pae vs fragments_min_pae
make_interactive_plot(full_data, 'fragments_min_pae', 'fl_min_pae', colour='input_length', size='fragments_rop',
                      save_path='/Users/poppy/Dropbox/AF3_FL_Centriole_screen/fl_vs_fragments_min_pae.html')
# fl_iptm vs fragments_iptm
make_interactive_plot(full_data, 'fragments_iptm', 'fl_iptm', colour='input_length', size='fragments_rop',
                      save_path='/Users/poppy/Dropbox/AF3_FL_Centriole_screen/fl_vs_fragments_iptm.html')

### PROCESS FOR MATCHING PROTEIN NAMES BETWEEN FL AND FRAGMENT DATA ###
# read in fragment data
fragment_input = '/Users/poppy/Dropbox/all_interface_analysis_2025.02.26.csv'
fragment_data = pd.read_csv(fragment_input)
# get all protein names from fragment data
fragment_names = set()
for index, row in fragment_data.iterrows():
    p1 = row['Protein1']
    p2 = row['Protein2']
    fragment_names.add(p1.lower())
    fragment_names.add(p2.lower())

# for all folders in FL directory, find protein names
fl_names = set()
for folder_name in os.listdir(directory):
    # Check if the folder is a directory
    if os.path.isdir(os.path.join(directory, folder_name)):
        # Get the protein names from the folder name
        p1, p2 = find_proteins(os.path.join(directory, folder_name))
        fl_names.add(p1.lower())
        fl_names.add(p2.lower())
print(fl_names)

# find proteins that are the same in both sets
print(fl_names.intersection(fragment_names))
# find proteins in the proteins set that are not in the fragment names set
# (ie included in FL screen but named differently to in fragment data)
print(fl_names.difference(fragment_names))

### ASSESS ONE FL PREDICTION ###
input = "/Users/poppy/Dropbox/AF3_FL_Centriole_screen/fold_ana3_fl_ana2_fl"
structure_file, json_file, log_file, PAE_png, fasta_file = find_rank_001_files(input, af3=True)
print(extract_iptm(log_file))
pae = extract_pae(json_file)
model = parse_structure_file(structure_file, is_pdb=False)
print(find_min_pae(model, pae))
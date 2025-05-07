# Making Chimera scripts for colouring interfaces/pae
###################################################################################################
# import necessary functions
import os
from collections import defaultdict
from analysis_utility import extract_pae, parse_structure_file, modify_bfactors

def create_visualisation(structure_file, json_file, residue_pairs, output_folder):
    pae = extract_pae(json_file)
    model = parse_structure_file(structure_file)

    #Finding min PAE for each interface residue and saving to pdb
    pae_values = defaultdict(lambda: float('inf'))  # Initialize with infinity
    for pair in residue_pairs:
        pae1 = min(pae[pair], pae[pair[::-1]])  # Get the minimum PAE for the pair
        if pae1 < pae_values[pair[0]]:
            pae_values[pair[0]] = pae1
        if pae1 < pae_values[pair[1]]:
            pae_values[pair[1]] = pae1
    # find max key in pae_values
    max_res = max(pae_values.keys())
    # make list of values, with each either being pae value for that residue if residue involved in interface pair, or 30 if not
    pae_values = [pae_values[i] if i in pae_values else 30 for i in range(0, max_res+1)]
    # all values over 30 are set to 30
    pae_values = [30 if i > 30 else i for i in pae_values]
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_folder + 'pdbs/'):
        os.makedirs(output_folder + 'pdbs/')
    file_name = structure_file.split('/')[-1].split('_unrelaxed')[0]
    pdb_output = output_folder + 'pdbs/' + file_name + '.pdb'
    modify_bfactors(model, pae_values, output=pdb_output, default=30)

    #Creating Chimera script
    # define palettes
    cols = ["0,dark red:3,red:5,salmon:15,light pink",
            "0,dark blue:3,blue:5,royal blue:15,light sky blue",
            "0,rebecca purple:3,purple:5,medium purple:15,thistle",
            "0,tomato:3,dark orange:5,orange:15,light goldenrod yellow",
            "0,dark green:3,green:5,lime green:15,pale green"
    ]

    # find chain_ids in model
    chain_ids = [chain.id for chain in model.get_chains()]
    chain_colours = ""
    for i, chain in enumerate(chain_ids):
        chain_colours += f"color bfactor #1/{chain} palette \"{cols[i]}\" \n"

    # Chimera script content
    script_content = f"""
    close
    set bgColor white
    open {'pdbs/' + file_name + '.pdb'}
    {chain_colours}
    select @@bfactor<=7 & #1
    show sel atoms
    select clear
    """

    # Save the script
    cxc_output = output_folder + file_name + '.cxc'
    with open(cxc_output, 'w') as file:
        file.write(script_content)

########################
# # Example usage
# import pandas as pd
# interfaces_csv_path = "/Users/poppy/Dropbox/all_interface_analysis_2025.04.29.csv"
# interfaces_df = pd.read_csv(interfaces_csv_path)
# # initialise empty df to store data
# domain_pairs_df = pd.DataFrame(columns=['model_rank', 'log_file', 'residue_pairs'])
# # Filter dataframe to only create files for predicted interactions
# for index, row in interfaces_df.iterrows():
#     if row['min_pae'] < 5 and row['size'] > 5 and row['avg_pae'] < 15 and row['rop'] >= 2:
#         domain_pairs_df = pd.concat([domain_pairs_df, row.to_frame().T], ignore_index=True)
# domain_pairs_df['residue_pairs'] = domain_pairs_df['residue_pairs'].apply(eval)
# # combine all rows representing predicted interaction between same domain pairs
# domain_pairs_df = domain_pairs_df.groupby(['model_rank', 'log_file']).agg({'residue_pairs': lambda x: [item for sublist in x for item in sublist]}).reset_index()
# for index, row in domain_pairs_df.iterrows():
#     # get the structure file and json file from the row
#     model_rank = row['model_rank']
#     folder_path = os.path.dirname(row['log_file'])
#     # use model_rank to find structure file and json file
#     structure_file = None
#     json_file = None
#     for file in os.listdir(folder_path):
#         if f'rank_00{model_rank}' in file:
#             if file.endswith('.pdb'):
#                 structure_file = os.path.join(folder_path, file)
#             elif file.endswith('.json'):
#                 json_file = os.path.join(folder_path, file)
#         if structure_file and json_file:
#             break
#     create_visualisation(structure_file, json_file, row['residue_pairs'], output_folder="/Users/poppy/Dropbox/visualisation/")

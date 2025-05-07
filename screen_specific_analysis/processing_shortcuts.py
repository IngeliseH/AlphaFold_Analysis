import os
import pandas as pd
from interface_analysis import find_and_score_interfaces
from run_folders import process_all_predictions
from identify_run_failures import find_missing_structure_files
from folder_cleanup import tidy_folder_structure
from format_for_gephi import gephi_format
from promiscuity_analysis import find_interface_promiscuity

def tidy_and_run(protein_names, to_tidy=True, find_missing=True, run_analysis=True, gephi_formatting=True):
    # check if protein_names is a list of protein names or a string
    if isinstance(protein_names, str):
        protein_names = [protein_names]
    for protein_name in protein_names:
        base_folder = f'/Users/poppy/Dropbox/{protein_name}'
        if to_tidy:
            tidy_folder_structure(base_folder)
        if find_missing:
            find_missing_structure_files(base_folder, save_to_csv=True)
        if run_analysis:
            process_all_predictions(base_folder, find_and_score_interfaces, output_file=f'{protein_name}_interface_analysis.csv', ipTM_graphic=False)
        if gephi_formatting:
            criteria = {
                'rop': lambda x: x >= 2,
                'avg_pae': lambda x: x <= 15,
                'min_pae': lambda x: x <= 5,
                'size': lambda x: x >= 5
            }
            gephi_format(f'/Users/poppy/Dropbox/{protein_name}/{protein_name}_interface_analysis.csv', include=['avg_pae', 'iptm', 'pdockq', 'min_pae', 'size'], criteria=criteria)
# tidy_and_run(['centriole_screen'], to_tidy=False)

def combine_data(protein_names):
    import datetime
    missing_files = []
    analysis_results = []
    date = datetime.datetime.now().strftime('%Y.%m.%d')
    for protein_name in protein_names:
        base_folder = f'/Users/poppy/Dropbox/{protein_name}'
        # find missing_structure_files.csv and read in data (4 columns, Protein1, Protein2, Protein1_Domain, Protein2_Domain)
        missing_file = os.path.join(base_folder, 'missing_structure_files.csv')
        if os.path.exists(missing_file):
            missing_files.append(pd.read_csv(missing_file))
        else:
            print(f"No missing_structure_files.csv found in {base_folder}")
        analysis = os.path.join(base_folder, f'{protein_name}_interface_analysis.csv')
        if os.path.exists(analysis):
            analysis_results.append(pd.read_csv(analysis))
        else:
            print(f"No {protein_name}_interface_analysis.csv found in {base_folder}")
    if missing_files:
        all_missing_files = pd.concat(missing_files)
        # don't need date here as not important to know previously missing files, keeps clear which is most recent
        all_missing_files.to_csv('/Users/poppy/Dropbox/all_run_failures.csv', index=False)
    if analysis_results:
        # combine all analysis results from subfolders
        all_analysis_results = pd.concat(analysis_results)
        # add promiscuity analysis
        all_analysis_results = find_interface_promiscuity(all_analysis_results)
        # name file with todays date to avoid overwriting and distinguish versions
        all_analysis_results.to_csv(f'/Users/poppy/Dropbox/all_interface_analysis_{date}.csv', index=False)
        criteria = {
            'rop': lambda x: x >= 2,
            'avg_pae': lambda x: x <= 15,
            'min_pae': lambda x: x <= 5,
            'size': lambda x: x >= 5
        }
        # format all data for gephi
        gephi_format(input_csv=f'/Users/poppy/Dropbox/all_interface_analysis_{date}.csv', include=['avg_pae', 'iptm', 'pdockq', 'min_pae', 'size', 'p1d_promiscuity', 'p2d_promiscuity'], output_csv=f'all_gephi_input_{date}.csv', criteria=criteria)


#protein_list = ['PCM', 'centriole_screen', 'Cdk2', 'CENPC', 'Centrin', 'Centrobin', 'Cep97', 'Cep135', 'CHC', 'CID', 'Cnn', 'CP110', 'Ctp', 'CycA', 'CycB1', 'CycB3', 'CycE', 'FZY', 'GCP2', 'GCP3', 'GCP4', 'GCP5', 'Grip71', 'MAD1_MAD2', 'MPS1', 'Mud', 'Mzt1']
#combine_data(protein_list)
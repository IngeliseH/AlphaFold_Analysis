import pandas as pd
import os

# function for processing a single domain pair (to be used with process_all_predictions function from run_folders script)
def find_incomplete_jobs(single_folder_path):
    # get the number of complete models from the log file
    log_file_path = os.path.join(single_folder_path, 'log.txt')
    if os.path.exists(log_file_path):
        with open(log_file_path, 'r') as log_file:
            lines = log_file.readlines()
            num_complete = 0
            num_complete = sum(1 for line in lines if '(7 recycles)' in line)
            done_found = any('Done' in line for line in lines)
            failure = False
            if done_found and num_complete != 5: # job failure rather than slow runtime
                num_complete = 0
                failure = True
    else:
        print(f"Warning: Log file {log_file_path} does not exist. Assuming all models are incomplete.")
        num_complete = 0
    
    return {'folder_path': single_folder_path, 'num_complete': num_complete, 'failure': failure}

# find 'incomplete_jobs_status.csv' file for each folder and concatenate them into one df
def concatenate_csv_files(folder_paths, status_file='incomplete_jobs_status.csv'):
    all_dataframes = []
    for folder_path in folder_paths:
        csv_file_path = os.path.join(folder_path, status_file)
        if os.path.exists(csv_file_path):
            df = pd.read_csv(csv_file_path)
            all_dataframes.append(df)
        else:
            print(f"Warning: CSV file {csv_file_path} does not exist.")
    
    if all_dataframes:
        combined_df = pd.concat(all_dataframes, ignore_index=True)
        return combined_df

# split into 5 dataframes based on number of complete models
def split_dataframe(df):
    df_0 = df[df['num_complete'] == 0]
    df_0_failure = df_0[df_0['failure'] == True]
    df_0_slow = df_0[df_0['failure'] == False]
    print(f"Number of jobs with 0 runs complete and failure: {len(df_0_failure)}")
    print(f"Number of jobs with 0 runs complete without failure: {len(df_0[df_0['failure'] == False])}")
    df_1 = df[df['num_complete'] == 1]
    print(f"Number of jobs with 1 run complete: {len(df_1)}")
    df_2 = df[df['num_complete'] == 2]
    print(f"Number of jobs with 2 runs complete: {len(df_2)}")
    df_3 = df[df['num_complete'] == 3]
    print(f"Number of jobs with 3 runs complete: {len(df_3)}")
    df_4 = df[df['num_complete'] == 4]
    print(f"Number of jobs with 4 runs complete: {len(df_4)}")
    df_5 = df[df['num_complete'] == 5]
    print(f"Number of jobs with 5 runs complete: {len(df_5)}")
    return df_0_failure, df_0_slow, df_1, df_2, df_3, df_4, df_5

# make a new folder for reruns and copy the fasta files into it
def prep_rerun_folder(df, output_folder_name):
    # join output folder name with 'Users/poppy/Dropbox/' to create the full path
    output_folder = os.path.join('/Users/poppy/Dropbox/', output_folder_name)
    # create a new folder for the reruns
    os.makedirs(output_folder, exist_ok=True)
    
    # iterate through the dataframe and copy the fasta files to the new folder
    for index, row in df.iterrows():
        folder_path = row['folder_path']
        # replace 4th part of folder path with the new folder name
        folder_parts = folder_path.split('/')
        protein_pair_folder = '/'.join(folder_parts[:6])
        folder_parts[4] = output_folder_name
        new_folder_path = '/'.join(folder_parts[:6])
        # create the new folder if it doesn't exist
        os.makedirs(new_folder_path, exist_ok=True)
        
        # copy the fasta file to the new folder
        fasta_name = folder_parts[6].replace('_output', '.fasta')

        fasta_file_path = os.path.join(protein_pair_folder, fasta_name)
        second_fasta_file_path = os.path.join(folder_path, fasta_name)
        new_fasta_file_path = os.path.join(new_folder_path, fasta_name)
        if os.path.exists(fasta_file_path):
            with open(fasta_file_path, 'r') as fasta_file:
                with open(new_fasta_file_path, 'w') as new_fasta_file:
                    new_fasta_file.write(fasta_file.read())
        elif os.path.exists(second_fasta_file_path):
            with open(second_fasta_file_path, 'r') as fasta_file:
                with open(new_fasta_file_path, 'w') as new_fasta_file:
                    new_fasta_file.write(fasta_file.read())
        else:
            print(f"Warning: Fasta file {fasta_file_path} does not exist.")

################################################################################
from run_folders import process_all_predictions

process_all_predictions('/Users/poppy/Dropbox/Centriole_dimers', find_incomplete_jobs, output_file='incomplete_jobs_status.csv', ipTM_graphic=False)
process_all_predictions('/Users/poppy/Dropbox/PCM_dimers', find_incomplete_jobs, output_file='incomplete_jobs_status.csv', ipTM_graphic=False)
process_all_predictions('/Users/poppy/Dropbox/more_dimers', find_incomplete_jobs, output_file='incomplete_jobs_status.csv', ipTM_graphic=False)

total_status = concatenate_csv_files(['/Users/poppy/Dropbox/Centriole_dimers', '/Users/poppy/Dropbox/PCM_dimers', '/Users/poppy/Dropbox/more_dimers'])

df_0_failure, df_0_slow, df_1, df_2, df_3, df_4, df_5 = split_dataframe(total_status)
# print out the names of the folders with 0 runs complete
print("Folders with 0 runs complete without failure:")
for folder in df_0_slow['folder_path']:
    print(folder)

# prep the rerun folders
prep_rerun_folder(df_0_slow, 'dimer_reruns_0_slow') # to be run with colabfold_folders_long.sh
prep_rerun_folder(df_0_failure, 'dimer_reruns_0_failure') # to be run with colabfold_folders.sh
prep_rerun_folder(df_1, 'dimer_reruns_1') # to be run with colabfold_folders_4_models_long.sh
prep_rerun_folder(df_2, 'dimer_reruns_2') # to be run with colabfold_folders_3_models_long.sh
prep_rerun_folder(df_3, 'dimer_reruns_3') # to be run with colabfold_folders_2_models.sh
prep_rerun_folder(df_4, 'dimer_reruns_4') # to be run with colabfold_folders_1_model.sh
# no need to make a rerun folder for df_5 as these are all complete

##################################################################################
import regex as re
from datetime import datetime

def combine_new_and_old_models(new_folder_path, old_folder_path):
    
    def update_model_number(file_name, folder_path, old_count):
        file_path = os.path.join(folder_path, file_name)
        match = re.search(r'model_(\d+)', file_name)
        if not match:
            print(f"Error: No model number found in {file_name}.")
            return file_name
        model_number = int(match.group(1))
        new_file_name = file_name.replace(f'model_{model_number}', f'model_{model_number + old_count}')
        new_file_path = os.path.join(folder_path, new_file_name)
        if os.path.exists(new_file_path):
            print(f"Warning: {new_file_name} already exists in {folder_path}.")
            return file_name
        os.rename(file_path, new_file_path)
    
    # count number of pdb files in old folder
    old_pdb_files = [f for f in os.listdir(old_folder_path) if f.endswith('.pdb')]
    old_pdb_count = len(old_pdb_files)
    
    new_pdb_files = []
    new_json_files = []
    # rename/delete files duplicate files in new folder
    for f in os.listdir(new_folder_path):
        # split name from extension
        name, ext = os.path.splitext(f)
        if (ext in ['.fasta', '.bibtex']) or ('.done' in name):
            os.remove(os.path.join(new_folder_path, f))
        elif (ext in ['.png', '.a3m']) or (name in ['.DS_Store', 'log', 'config']) or ('_env' in name) or ('_pairgreedy' in name) or ('predicted_aligned_error_v1' in name):
            new_name = name + '_2' + ext
            os.rename(os.path.join(new_folder_path, f), os.path.join(new_folder_path, new_name))
        elif ext == '.pdb':
            new_pdb_files.append(f)
        elif ext == '.json' and '_model' in name:
            new_json_files.append(f)
    print('new_pdb_files:', new_pdb_files)
    print('new_json_files:', new_json_files)
    
    # should iterate from high model number to low to avoid overwriting
    new_pdb_files.sort(key=lambda x: int(re.search(r'model_(\d+)', x).group(1)), reverse=True)
    new_json_files.sort(key=lambda x: int(re.search(r'model_(\d+)', x).group(1)), reverse=True)
    for f in new_pdb_files:
        update_model_number(f, new_folder_path, old_pdb_count)
    for f in new_json_files:
        update_model_number(f, new_folder_path, old_pdb_count)
    
    model_dict = {1:{'iptm': None, 'ptm':None, 'plddt': None, 'rank': None}, 
                  2:{'iptm': None, 'ptm':None, 'plddt': None, 'rank': None}, 
                  3:{'iptm': None, 'ptm':None, 'plddt': None, 'rank': None}, 
                  4:{'iptm': None, 'ptm':None, 'plddt': None, 'rank': None}, 
                  5:{'iptm': None, 'ptm':None, 'plddt': None, 'rank': None}}

    new_log_file = os.path.join(new_folder_path, 'log_2.txt')
    with open(new_log_file, 'r') as log_file:
        lines = log_file.readlines()
        for i, line in enumerate(lines):
            model_match = re.search(r'alphafold2_multimer_v3_model_(\d)_seed_000', line)
            final_match = re.search(r'model_(\d+)_seed_000 recycle=7 pLDDT=(\d+\.\d+) pTM=(\d+\.\d+) ipTM=(\d+\.\d+) tol=(\d+\.\d+)', line)
            if model_match:
                original_number = int(model_match.group(1))
                model_number = original_number + old_pdb_count
                lines[i] = line.replace(f'model_{original_number}', f'model_{model_number}')
                if final_match:
                    model_dict[model_number]['iptm'] = float(final_match.group(4))
                    model_dict[model_number]['ptm'] = float(final_match.group(3))
                    model_dict[model_number]['plddt'] = float(final_match.group(2))
                
        with open(new_log_file, 'w') as log_file:
            log_file.writelines(lines)
    print(model_dict)

    # move all files and folders from new folder to old folder
    for f in os.listdir(new_folder_path):
        current_file_path = os.path.join(new_folder_path, f)
        moving_to_file_path = os.path.join(old_folder_path, f)
        if os.path.exists(moving_to_file_path):
            print(f"Warning: {f} already exists in {old_folder_path}.")
        else:
            os.rename(current_file_path, moving_to_file_path)
    
    # go through the log file and find the lines with the new model numbers
    original_log_file = os.path.join(old_folder_path, 'log.txt')
    with open(original_log_file, 'r') as log_file:
        lines = log_file.readlines()
        for line in lines:
           final_match = re.search(r'model_(\d+)_seed_000 recycle=7 pLDDT=(\d+\.\d+) pTM=(\d+\.\d+) ipTM=(\d+\.\d+) tol=(\d+\.\d+)', line)
           if final_match:
                model_number = int(final_match.group(1))
                model_dict[model_number]['iptm'] = float(final_match.group(4))
                model_dict[model_number]['ptm'] = float(final_match.group(3))
                model_dict[model_number]['plddt'] = float(final_match.group(2))
    model_dict = {k: v for k, v in sorted(model_dict.items(), key=lambda item: item[1]['iptm'], reverse=True)}
    print(model_dict)
    # add ranks to dict
    for i, (model_number, model_data) in enumerate(model_dict.items()):
        model_data['rank'] = i + 1

    # add ranks to pdb and json file names
    for f in os.listdir(old_folder_path):
        if f.endswith('.pdb'):
            match = re.search(r'_unrelaxed(.*?)_multimer_v3_model_(\d+)_seed', f)
            if match:
                model_number = int(match.group(2))
                rank = model_dict[model_number]['rank']
                new_name = f.replace(match.group(1), f'_00{rank}_alphafold2')
                os.rename(os.path.join(old_folder_path, f), os.path.join(old_folder_path, new_name))
        elif f.endswith('.json') and 'predicted_aligned_error' not in f:
            match = re.search(r'_scores(.*?)_multimer_v3_model_(\d+)_seed', f)
            if match:
                model_number = int(match.group(2))
                rank = model_dict[model_number]['rank']
                new_name = f.replace(match.group(1), f'_00{rank}_alphafold2')
                os.rename(os.path.join(old_folder_path, f), os.path.join(old_folder_path, new_name))
        
    # generate ranking paragraph
    now = datetime.now()
    date_time = now.strftime("%Y-%m-%d %H:%M:%S")
    ranking_paragraph = f"{date_time} reranking models by 'ipTM' metric\n{date_time} ranking has been added on {date_time}\n{date_time} this is due to merging of an initial incomplete prediction run generating {old_pdb_count} models with a second run generating {5-old_pdb_count} models\n"
    for i, (model_number, model_data) in enumerate(model_dict.items()):
        ranking_paragraph += f"{date_time} rank_00{model_data['rank']}_alphafold2_multimer_v3_model_{model_number}_seed_000 pLDDT={model_data['plddt']} pTM={model_data['ptm']} ipTM={model_data['iptm']}\n"
    # add ranking paragraph to end of log file
    new_log_file = os.path.join(old_folder_path, 'log_2.txt')
    with open(original_log_file, 'a') as log_file:
        # add contents of log_2.txt to end of log.txt
        with open(new_log_file, 'r') as new_log_file:
            log_file.write(new_log_file.read())
        # add ranking paragraph to end of log file
        log_file.write(ranking_paragraph)

old_folder_path = "/Users/poppy/Dropbox/test/old/Ana1_F2+Asl_dimer_F1_output"
new_folder_path = "/Users/poppy/Dropbox/test/new/Ana1_F2+Asl_dimer_F1_output"
combine_new_and_old_models(new_folder_path, old_folder_path)
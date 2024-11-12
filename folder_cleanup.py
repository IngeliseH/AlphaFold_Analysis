"""
Scripts to clean up folder structure by moving files from subdirectories to the parent directory and removing
empty directories, or removing unwanted files based on a regex pattern to save space.

Functions:
    - tidy_folder_structure
    - delete_files_matching_pattern
"""

import os
import re
from pathlib import Path
import shutil

def tidy_folder_structure(base_folder):
    """
    Move files from subdirectories to the parent directory and remove empty directories.

    Parameters:
        - base_folder (str): The path to the base directory containing subdirectories to be tidied.

    Returns:
        None

    Note:
        Once a folder is correctly formatted, rerunning the function will not move any files or
        remove any directories.
    """
    # Go through each protein pair folder in the base folder
    for protein_pair_folder in os.listdir(base_folder):
        protein_pair_path = os.path.join(base_folder, protein_pair_folder)
        if os.path.isdir(protein_pair_path):
            # Go through each fragment pair folder in the protein pair folder
            for fragment_pair_folder in os.listdir(protein_pair_path):
                fragment_pair_path = os.path.join(protein_pair_path, fragment_pair_folder)
                if os.path.isdir(fragment_pair_path):
                    # Now traverse down to find the colabfold_*_tmp folder
                    for root, dirs, files in os.walk(fragment_pair_path):
                        if os.path.basename(root).startswith('colabfold_') and root.endswith('_tmp'):
                            # Move each file from the colabfold_*_tmp folder to the fragment pair folder
                            for file in files:
                                file_path = os.path.join(root, file)
                                destination_path = os.path.join(fragment_pair_path, file)
                                
                                shutil.move(file_path, destination_path)
                                print(f'Moved file {file_path} to {destination_path}')
                            
                            # Move each subdirectory from the colabfold_*_tmp folder to the fragment pair folder
                            for dir_name in dirs:
                                dir_path = os.path.join(root, dir_name)
                                destination_path = os.path.join(fragment_pair_path, dir_name)
                                
                                shutil.move(dir_path, destination_path)
                    
                    # Remove empty directories within the fragment pair path
                    for root, dirs, files in os.walk(fragment_pair_path, topdown=False):
                        if not os.listdir(root):
                            try:
                                os.rmdir(root)
                                print(f'Removed empty directory {root}')
                            except OSError as e:
                                print(f"Error removing directory {root}: {e}")

def delete_files_matching_pattern(folder_path, pattern):
    """
    Delete files within a folder that match a given regex pattern.

    Parameters
        - folder_path (str): The path to the directory containing files to be deleted
        - pattern (str): The regex pattern to match the file names against

    Returns
        None
    """
    folder = Path(folder_path)
    total_size = 0
    deleted_files = 0

    # Compile the regex pattern
    regex = re.compile(pattern)

    # Recursively search for files matching the regex pattern
    for file in folder.rglob('*'):
        # Check if the file name matches the regex
        if regex.search(str(file)):
            file_size = file.stat().st_size
            file.unlink()
            total_size += file_size
            deleted_files += 1

    print(f"Deleted {deleted_files} files. Total space cleared: {total_size} bytes")

    return None

# Example usage
#folder_path = "../../../../../Dropbox/2022.10.20_Drosophila_Version_1"

# to tidy folder structure
#tidy_folder_structure(folder_path)


# to delete all .fasta files in the folder
#fasta_pattern = r".*\.fasta$"  # Regex to match any file that ends with '.fasta'
#delete_files_matching_pattern(folder_path, fasta_pattern)

# to delete all .json files in the folder that do not contain 'rank_001'
#json_pattern = r"^(?!.*rank_001).*\.json$"  # Regex to match '.json' files not containing 'rank_001'
#delete_files_matching_pattern(folder_path, json_pattern)


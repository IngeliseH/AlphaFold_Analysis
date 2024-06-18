import re
from pathlib import Path

def delete_files_matching_pattern(folder_path, pattern):
    """
    Delete files within a folder that match a given regex pattern.

    Parameters
        - folder_path (str): The path to the directory containing files to be deleted
        - pattern (str): The regex pattern to match the file names against
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


# to delete all .fasta files in the folder
#fasta_pattern = r".*\.fasta$"  # Regex to match any file that ends with '.fasta'
#delete_files_matching_pattern(folder_path, fasta_pattern)

# to delete all .json files in the folder that do not contain 'rank_001'
#json_pattern = r"^(?!.*rank_001).*\.json$"  # Regex to match '.json' files not containing 'rank_001'
#delete_files_matching_pattern(folder_path, json_pattern)


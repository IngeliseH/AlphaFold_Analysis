"""
Functions to extract data from the files generated by AlphaFold

Functions:
find_rank_001_files
extract_fasta_protein_lengths
parse_structure_file
map_chains_and_residues
determine_chain_lengths
extract_pae
"""
import json
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.PDB import MMCIFParser

def find_rank_001_files(folder_path, af3=False):
    """
    Identifies the highest ranked model files in the given folder, as well as the PAE
    PNG and FASTA files if using AF2.

    Parameters:
        - folder_path (str): Path to the folder containing the model files.
        - af3 (bool): Whether the files are from AlphaFold3 (default: False).
    
    Returns:
        - tuple: Path to the PDB file, path to the JSON file, path to the log file, path to the PAE PNG file,
          path to the FASTA file.
    
    Note:
        - The af3 argument can be unset even if the files are from AlphaFold3 - if no
            PDB file is found, the function will automatically look for a CIF file. If
            files are from AF3, setting af3 to True will marginally increase speed.
    """
    folder = Path(folder_path)
    structure_file, json_file, log_file, PAE_png, fasta_file = None, None, None, None, None
    if not af3:
        # AlphaFold2 file handling
        for file in folder.glob('*'):
            if 'rank_001' in file.stem or 'rank_1' in file.stem:
                if file.suffix == '.pdb':
                    structure_file = file
                elif file.suffix == '.json':
                    json_file = file
            if 'log' in file.stem and file.suffix == '.txt':
                log_file = file
            if ('pae' in file.stem or 'PAE' in file.stem) and file.suffix == '.png':
                PAE_png = file
            if file.suffix == '.fasta':
                fasta_file = file
    # AF3 will have cif file instead of pdb file
    if af3 == True or not structure_file:
        for file in folder.glob('*'):
            if 'model_0' in file.stem and file.suffix == '.cif':
                structure_file = file
            elif 'full_data_0' in file.stem and file.suffix == '.json':
                json_file = file
            if 'summary_confidences_0' in file.stem and file.suffix == '.json':
                log_file = file
    return structure_file, json_file, log_file, PAE_png, fasta_file

def extract_fasta_protein_lengths(fasta_file):
    """
    Extracts the lengths of the two proteins from the FASTA file.
    Assumes the fasta file contains 2 sequences separated by ':'.

    Parameters:
        - fasta_file (str): Path to the FASTA file.

    Returns:
        - tuple: Length of the first protein, length of the second protein.
    """
    with open(fasta_file, 'r') as file:
        lines = file.readlines()  # Read all lines in the file
        fasta_content = ''.join(lines[1:])  # Join lines starting from the second line
    sequences = fasta_content.strip().split(':')
    return len(sequences[0]), len(sequences[1])

def parse_structure_file(file_path, is_pdb=True):
    """
    Parses the structure file (PDB or CIF) and returns the model.
    
    Parameters:
        - file_path (str): Path to the structure file.
        - is_pdb (bool): Whether the file is in PDB format (default: True).
    """
    parser = PDBParser() if is_pdb else MMCIFParser()
    structure = parser.get_structure('protein', file_path)
    model = structure[0]
    return model

def map_chains_and_residues(model):
    """
    Maps the chains and residues in the structure model to their absolute residue IDs.
    Absolute residue IDs are unique across all chains.
    
    Parameters:
        - model (Bio.PDB.Model.Model): The structure model.
        
    Returns:
        - chain_residue_map (list): List of tuples representing chain, residue ID, residue
          name, and absolute residue ID.
    """
    chain_residue_map = []
    unique_residues = set()
    abs_res_id = 0

    for chain in model.get_chains():
        chain_id = chain.id
        for residue in chain.get_residues():
            res_id = residue.id[1]
            res_name = residue.resname
            residue_tuple = (chain_id, res_id, res_name)
            if residue_tuple not in unique_residues:
                chain_residue_map.append((chain_id, res_id, res_name, abs_res_id))
                unique_residues.add(residue_tuple)
                abs_res_id += 1
    return sorted(chain_residue_map)

def determine_chain_lengths(model):
    """
    Determines the length of each chain in the given Biopython Model object.

    Parameters:
        - model (Bio.PDB.Model.Model): A model object from Biopython containing the chains.

    Returns:
        - list: A list of chain lengths in the order they appear.
    """
    chain_lengths = []
    # Iterate through each chain in the model
    for chain in model:
        count = 0
        for _ in chain.get_residues():
            count += 1
        chain_lengths.append(count)  # Append the count of residues in the chain

    return chain_lengths

def extract_pae(json_file):
    """
    Parses the PAE JSON file to extract the PAE matrix.

    Parameters:
        - pae_file (str): Path to the PAE JSON file.

    Returns:
        - list: PAE matrix.
    """
    with open(json_file, 'r') as file:
        data = json.load(file)
    if isinstance(data, list):
        # If data is a list, extract the first element
        pae = data[0].get('pae', data[0].get('predicted_aligned_error', 'Error: PAE not found'))
    elif isinstance(data, dict):
        # If data is a dictionary
        pae = data.get('pae', data.get('predicted_aligned_error', 'Error: PAE not found'))
    else:
        raise ValueError("Unexpected data format in JSON file")
        
    pae_matrix = np.array(pae)
    return pae_matrix

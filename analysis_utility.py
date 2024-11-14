"""
Functions to extract data from the files generated by AlphaFold

Functions:
find_rank_001_files
extract_fasta_protein_lengths
parse_structure_file
map_chains_and_residues
determine_chain_lengths
extract_pae
modify_bfactors
readable_ranges
"""
import json
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser, MMCIFParser, PDBIO

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
        - file_path (str or Path): Path to the structure file.
        - is_pdb (bool): Whether the file is in PDB format (default: True).

    Returns:
        - model: The first model in the structure.
    """
    parser = PDBParser() if is_pdb else MMCIFParser()
    # Convert file_path to string if it's a Path object
    if isinstance(file_path, Path):
        file_path = str(file_path)
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
    chain_lengths = [len(list(chain.get_residues())) for chain in model]
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

def modify_bfactors(protein_model, residue_scores_list, output, default=0.0):
    """
    Modifies the B-factor column of a PDB file to store custom scores.

    Parameters:
        - protein_model (Model): The protein structure model.
        - residue_scores_list (list): List of scores for each residue.
        - output (str or Path): Name + file path for the output PDB file.
        - default (float): Default score to use if the list is shorter than the number of residues.
    Returns:
        None, but saves the modified PDB file.
    """
    # Flatten the residue similarity scores into a single list
    all_scores = []
    for scores in residue_scores_list:
        all_scores.extend(scores)
    
    residue_counter = 0
    for chain in protein_model.get_chains():
        for residue in chain.get_residues():
            for atom in residue.get_atoms():
                if residue_counter < len(all_scores):
                    atom.bfactor = all_scores[residue_counter]
                else:
                    atom.bfactor = default  # Default B-factor if out of scores
            residue_counter += 1
    
    # Save the modified structure
    io = PDBIO()
    io.set_structure(protein_model)

    # Ensure output_name is a string
    if isinstance(output, Path):
        output = str(output)
    io.save(output)

def readable_ranges(absolute_pairs, structure_model, abs_res_lookup_dict=None):
    """
    Convert absolute residue IDs to chain-specific IDs and format ranges by protein in an easily readable way.

    Parameters:
        - absolute_pairs (list of tuples): List of residue pairs in absolute numbering.
        - structure_model (Bio.PDB.Model.Model): Protein structure model containing chain IDs.
        - abs_res_lookup_dict (dict, optional): Precomputed dictionary mapping chain-specific
          to absolute residue IDs. If None, this dictionary is computed within the function.

    Returns:
        - dict: Dictionary containing formatted residue ranges for each protein.
    """
    # Compute the reverse lookup dictionary if not provided
    if abs_res_lookup_dict is None:
        chain_residue_map = map_chains_and_residues(structure_model)
        abs_res_lookup_dict = {(chain_id, res_id): abs_res_id for chain_id, res_id, _, abs_res_id in chain_residue_map}

    # Reverse lookup to map absolute IDs to chain-specific IDs
    inv_abs_res_lookup_dict = {v: k for k, v in abs_res_lookup_dict.items()}
    chain_ids = [chain.id for chain in structure_model.get_chains()]

    # Convert each pair to chain-specific IDs
    chain_specific_pairs = []
    residues_by_protein = {'Protein1': set(), 'Protein2': set()}

    for res1, res2 in absolute_pairs:
        chain1, res_num1 = inv_abs_res_lookup_dict[res1]
        chain2, res_num2 = inv_abs_res_lookup_dict[res2]
        chain_specific_pairs.append((f"{chain1} {res_num1}", f"{chain2} {res_num2}"))

        # Group residues by protein based on chain
        if chain1 == chain_ids[0]:  # Assuming chain1 corresponds to Protein1
            residues_by_protein['Protein1'].add(res_num1)
            residues_by_protein['Protein2'].add(res_num2)
        else:
            residues_by_protein['Protein1'].add(res_num2)
            residues_by_protein['Protein2'].add(res_num1)

    # Format residue ranges for simplified output
    def format_ranges(numbers):
        if not numbers:
            return "None"
        sorted_nums = sorted(numbers)
        ranges = []
        start = end = sorted_nums[0]
        for n in sorted_nums[1:]:
            if n <= end + 2:
                end = n
            else:
                ranges.append(f"{start}-{end}" if start != end else str(start))
                start = end = n
        ranges.append(f"{start}-{end}" if start != end else str(start))
        return ", ".join(map(str, ranges))

    # Add chain-specific numbering information as an attribute
    return {
        'Protein1': format_ranges(residues_by_protein['Protein1']),
        'Protein2': format_ranges(residues_by_protein['Protein2'])
    }

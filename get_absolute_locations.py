"""
use protein csv
for each interface, find chain groupings, protein names (minus dimer)
find fasta file (or recreate from pdb if not found)
make dict matching sequences from fasta to protein names, and to 'locations'
find full length sequences from uniprot or from proteins csv (use actual sequences due to issues with numbering in some fragmentations)
find fasta subsection in full length sequence, use this to get position shift
for each 'location', add corresponding position shift to get absolute location in full length protein
"""

import os
import requests
import pandas as pd
from ast import literal_eval

def fetch_sequence(uniprot_id):
    """
    Shorter version of the AlphaFragment Uniprot fetch function - this one has less error handling
    and specifically returns sequences rather than all data.
    """
    url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}"
    try:
        r = requests.get(url, headers={"Accept":"application/json"}, timeout=30)
        r.raise_for_status()
        data = r.json()
        return data.get('sequence', None)
    except Exception as e:
        print("Error fetching", uniprot_id, e)
        return None

def recreate_fasta(pdb_path, fasta_path):
    """
    Reads a PDB or cif file and writes a FASTA file of its amino acid sequence(s).
    Standard amino acids only; multiple chains are joined by ':' in the output sequence.

    Parameters:
      - pdb_path (str): Path to the input PDB/cif file.
      - fasta_path (str): Path to the output FASTA file.
    """
    # Mapping 3-letter amino acid codes to 1-letter codes
    _THREE_TO_ONE = {
        'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
        'GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I',
        'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
        'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'
    }

    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    chains = {}  # chain_id => list of (residue_id, resname)

    if pdb_path.endswith('.pdb'):
        with open(pdb_path, 'r') as f:
            for line in f:
                # Parse ATOM records for standard residues
                if line.startswith(('ATOM  ', 'HETATM')):
                    resname = line[17:20].strip()
                    chain_id = line[21].strip()
                    res_seq = line[22:26].strip()
                    key = (res_seq, line[26].strip())  # include insertion code
                    if resname in _THREE_TO_ONE:
                        chains.setdefault(chain_id, []).append((key, resname))

    elif pdb_path.endswith('.cif'):
        from gemmi import cif
        doc = cif.read(pdb_path)
        block = doc.sole_block()
        atom_site = block.find(['_atom_site.label_asym_id', '_atom_site.label_atom_id', '_atom_site.label_comp_id',
                                '_atom_site.label_seq_id', '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z'])
        for row in atom_site:
            resname = row[2].strip()
            chain_id = row[0].strip()
            res_seq = row[3].strip()
            key = (res_seq, '')
            if resname in _THREE_TO_ONE:
                chains.setdefault(chain_id, []).append((key, resname))

    # Build sequences, ensuring unique ordering per chain
    seqs = []
    for chain_id, residues in sorted(chains.items()):
        seen = set()
        seq = []
        for (res_seq, ins_code), resname in residues:
            identifier = (res_seq, ins_code)
            if identifier in seen:
                continue
            seen.add(identifier)
            seq.append(_THREE_TO_ONE[resname])
        seqs.append(''.join(seq))

    joined = ':'.join(seqs)

    header = os.path.splitext(os.path.basename(pdb_path))[0]

    with open(fasta_path, 'w') as out:
        out.write(f">{header}\n")
        out.write(joined + '\n')

def add_absolute_interface_locations(interface_df, protein_df):
    interface_data = pd.read_csv(interface_df)
    protein_data = pd.read_csv(protein_df)

    proteins_dict = {}
    for index, row in protein_data.iterrows():
        protein_name = row['name']
        uniprot_id = row['accession_id']
        sequence = None
        if uniprot_id and pd.notna(uniprot_id):
            sequence = fetch_sequence(uniprot_id)
        if not sequence:
            sequence = row['sequence']
        if not protein_name.endswith('_dimer'):
            proteins_dict[protein_name] = sequence

    for index, row in interface_data.iterrows():
        protein1 = row['Protein1']
        protein2 = row['Protein2']
        dimer1 = False
        if '_dimer' in protein1:
            protein1 = protein1.replace('_dimer', '')
            group1 = ('A', 'B')
        else:
            group1 = ('A')
        if '_dimer' in protein2:
            protein2 = protein2.replace('_dimer', '')
            group2 = ('B', 'C') if not dimer1 else ('C', 'D')
        else:
            group2 = ('B') if not dimer1 else ('C')
        chain_grouping = [group1, group2]

        folder = row['log_file'].split('/')
        file_name = folder[-2].replace('_output', '')
        fasta_file = '/'.join(folder[0:-1]) + '/' + file_name + '.fasta'
        if not os.path.exists(fasta_file):
            fasta_file = '/'.join(folder[0:-2]) + file_name + '.fasta'
        if not os.path.exists(fasta_file):
            print(f"Fasta file not found for {fasta_file}, trying to create from pdb.")
            pdb_folder = '/'.join(folder[0:-1])
            # look for any .pdb or .cif extension file in pdb_folder (any will do)
            pdb_file = None
            for file in os.listdir(pdb_folder):
                if file.endswith('.pdb') or file.endswith('.cif'):
                    pdb_file = os.path.join(pdb_folder, file)
                    break
            if pdb_file:
                fasta_file = '/'.join(folder[0:-2]) + file_name + '.fasta'
                recreate_fasta(pdb_file, fasta_file)
                print(f"Fasta file created at {fasta_file}.")
            else:
                print(f"No pdb or cif file found in {pdb_folder}, skipping.")            
            continue

        try:
            locations = literal_eval(row['location'])
        except (ValueError, SyntaxError):
            locations = None
        if not isinstance(locations, dict):
            continue
        # locations is a dict - get values as list
        locations = list(locations.values())

        # Get fragment sequences from fasta file, use chain grouping to assign to protein names
        with open(fasta_file, 'r') as f:
            lines = f.readlines()
            sequence = lines[1].strip()
            fragments = sequence.split(':')
            protein_names = []
            chains = []

            for group, protein in zip(chain_grouping, [protein1, protein2]):
                for chain in group:
                    chains.append(chain)
                    protein_names.append(protein)
            class Fragment:
                def __init__(self, chain, sequence, protein_name, location):
                    self.chain = chain
                    self.sequence = sequence
                    self.protein_name = protein_name
                    self.location = location
                    self.shift = None  # to be determined
                    self.absolute_location = None  # to be determined
            fragment_objects = []
            for chain, frag, prot, location in zip(chains, fragments, protein_names, locations):
                fragment_objects.append(Fragment(chain, frag, prot, location))

        # find fasta subsection in full length sequence, use this to get position shift
        for fragment in fragment_objects:
            full_sequence = proteins_dict.get(fragment.protein_name, None)
            if not full_sequence:
                print(f"Full length sequence not found for {fragment.protein_name}, skipping.")
                continue
            start_index = full_sequence.find(fragment.sequence)
            if start_index == -1:
                print(f"Fragment sequence not found in full length sequence for {fragment.protein_name}, skipping.")
                continue
            fragment.shift = start_index

            if fragment.location is None or fragment.location == 'None':
                fragment.absolute_location = 'None'
                continue
            loc_parts = fragment.location.split(',')
            abs_locs = []
            for part in loc_parts:
                part = part.strip()
                if '-' in part:
                    start, end = part.split('-')
                    start = int(start) + fragment.shift
                    end = int(end) + fragment.shift
                    abs_locs.append(f"{start}-{end}")
                else:
                    pos = int(part) + fragment.shift
                    abs_locs.append(f"{pos}")
            fragment.absolute_location = (', '.join(abs_locs))

        absolute_location_dict = {}
        for fragment in fragment_objects:
            if fragment.absolute_location:
                absolute_location_dict[f"Chain {fragment.chain}"] = fragment.absolute_location
        interface_data.at[index, 'absolute_location'] = str(absolute_location_dict)

    interface_data.to_csv(interface_df.replace('.csv', '_absolute.csv'), index=False)

#protein_df = "/Volumes/T7/screen_results/all_fragments_2025.06.04.csv"
#interface_df = "/Volumes/T7/screen_results/all_interface_analysis_2025.11.26.csv"
#add_absolute_interface_locations(interface_df, protein_df)

# Mapping interface locations onto FL sequences
#################################################################
"""
• provide mappings of protein names => uniprot ids/FL sequences
• if not got sequences, get from UniProt
• store in dict to avoid repeated queries
• define protein class
    • attributes = uniprot id, protein name, FL sequence, interfaces
    * methods = partners (finds interactions without having to determine if protein1 or 2)
• define interface class
    • attributes = protein1, protein2, res in protein1, res in protein2, interface_id, some measure of quality (eg minpae, avgpae, size), screen (FL, fragment, Isaac's)
• make protein class objects for each protein

• condense all_int dataframe to only good predictions (via criteria that can be varied in the future)
• condense further to one row per domain pair, with all interface positions (currently dicts in 'location' column) assembled into a list of dicts
• for each prediction:
    • use log_file location to navigate to fasta
    • get fragment sequences from fasta, splitting by ":" (don't make domains/fragments as objects, as may be small seq variations)
    • use Protein1 and Protein2 columns to get the FL sequence for each protein
    • get start and end positions of each fragment in the FL sequence by searching for fragment seq in the FL seq
    • use this to get the positions of each interface in the FL sequence - add shift based on the start position of the fragment, then use 'expand_loc' function to get the full range of residues
    • make 'interface' class objects for each interface
    • add interface to an overall list, and as an attribute of each protein
"""
import os
import pandas as pd
import requests
from analysis_utility import expand_loc_range

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Polygon
from matplotlib.colors import Normalize, is_color_like
from matplotlib.cm import ScalarMappable
import time

class Protein:
    def __init__(self, name, accession_id, sequence, category=None, screens=[]):
        self.name = name
        self.accession_id = accession_id
        self.sequence = sequence
        self.interfaces = []
        self.category = category
        # literal eval screens as list of strings
        if isinstance(screens, str):
            try:
                screens = eval(screens)
            except Exception as e:
                print(f"Error evaluating screens: {e}")
                screens = []
        self.screens = screens
    
    def add_interface(self, interface):
        # check interface is Interface object
        if not isinstance(interface, Interface):
            raise ValueError(f"Expected Interface object, got {type(interface)}")
        else:
            self.interfaces.append(interface)

    def partners(self):
        """
        Yield one tuple per interface:
          (partner_protein,
           my_residue_list,
           their_residue_list,
           min_pae,
           avg_pae,
           size,
           interface_id,
           category_id,
           screen)
        """
        for iface in self.interfaces:
            if iface.protein1 is self:
                yield (
                    iface.protein2,
                    iface.res1,       # residues on self
                    iface.res2,       # residues on partner
                    iface.min_pae,
                    iface.avg_pae,
                    iface.size,
                    iface.interface_id,
                    iface.category_id,
                    iface.screen
                )
            elif iface.protein2 is self:
                yield (
                    iface.protein1,
                    iface.res2,       # residues on self
                    iface.res1,       # residues on partner
                    iface.min_pae,
                    iface.avg_pae,
                    iface.size,
                    iface.interface_id,
                    iface.category_id,
                    iface.screen
                )
            else:
                # this should never happen if you only add valid interfaces
                raise ValueError(f"{self.name} not found in interface {iface.interface_id}")

class Interface:
    def __init__(self, protein1, protein2, res1, res2, min_pae=None, avg_pae=None, size=None, screen=None):
        # make sure protein1 and 2 are alphabetically ordered
        if protein1.name > protein2.name:
            protein1, protein2 = protein2, protein1
            res1, res2 = res2, res1
        # check protein1 and 2 are Protein objects
        if not isinstance(protein1, Protein):
            raise ValueError(f"Expected Protein object, got {type(protein1)}")
        if not isinstance(protein2, Protein):
            raise ValueError(f"Expected Protein object, got {type(protein2)}")
        # check res1 and 2 are lists
        if not isinstance(res1, list) and not (isinstance(r, int) for r in res1):
            if isinstance(res1, int):
                res1 = [res1]
            elif isinstance(res1, list) and all(isinstance(r, str) for r in res1):
                res1 = [int(r) for r in res1]
            else:
                raise ValueError(f"Expected list of ints, got {type(res1)}")
        if not isinstance(res2, list) and not (isinstance(r, int) for r in res2):
            if isinstance(res2, int):
                res2 = [res2]
            elif isinstance(res2, list) and all(isinstance(r, str) for r in res2):
                res2 = [int(r) for r in res2]
            else:
                raise ValueError(f"Expected list of ints, got {type(res2)}")
        # check min_pae, avg_pae are floats and size is int
        if min_pae is not None and not isinstance(min_pae, (float, int)):
            if isinstance(min_pae, str):
                try:
                    min_pae = float(min_pae)
                except ValueError:
                    raise ValueError(f"Expected float or int, got {type(min_pae)}, with value {min_pae}")
            elif isinstance(min_pae, list) and len(min_pae) == 1 and isinstance(min_pae[0], (float, int)):
                min_pae = float(min_pae[0])
            else:
                raise ValueError(f"Expected float or int, got {type(min_pae)}, with value {min_pae}")
        if avg_pae is not None and not isinstance(avg_pae, (float, int)):
            if isinstance(avg_pae, str):
                try:
                    avg_pae = float(avg_pae)
                except ValueError:
                    raise ValueError(f"Expected float or int, got {type(avg_pae)}, with value {avg_pae}")
            elif isinstance(avg_pae, list) and len(avg_pae) == 1 and isinstance(avg_pae[0], (float, int)):
                avg_pae = float(avg_pae[0])
            else:
                raise ValueError(f"Expected float or int, got {type(avg_pae)}, with value {avg_pae}")
        if size is not None and not isinstance(size, (float, int)):
            if isinstance(size, (float, str)):
                try:
                    size = int(size)
                except ValueError:
                    raise ValueError(f"Expected int, got {type(size)}, with value {size}")
            elif isinstance(size, list) and len(size) == 1 and isinstance(size[0], (float, int)):
                size = int(size[0])
            else:
                raise ValueError(f"Expected int, got {type(size)}, with value {size}")
        self.protein1 = protein1
        self.protein2 = protein2
        self.res1 = res1
        self.res2 = res2
        self.min_pae = min_pae
        self.avg_pae = avg_pae
        self.size = size
        self.interface_id = f"{protein1.name}_{protein2.name}"
        if protein1.category < protein2.category:
            self.category_id = f"{protein1.category}_{protein2.category}"
        else:
            self.category_id = f"{protein2.category}_{protein1.category}"
        self.screen = screen

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

def condense_interfaces(df, criteria):
    """
    Condense the interface dataframe to one row per domain pair, with all interface positions
    assembled into a list of dicts.

    Parameters:
      - df (pd.DataFrame): The dataframe containing interface data.
    """
    # convert location column to dict (instead of string representation)
    df['location'] = df['location'].apply(lambda x: eval(x) if isinstance(x, str) else x)
    # group by domain pairs and aggregate interface positions
    condensed_df = df.groupby(['Protein1', 'Protein2', 'Protein1_Domain', 'Protein2_Domain', 'log_file']).agg({
        'location': lambda x: list(x),
        'min_pae': lambda x: list(x),
        'avg_pae': lambda x: list(x),
        'size': lambda x: list(x)
    }).reset_index()
    return condensed_df

def find_pos_shift(protein1, protein2, log_file):
    """
    Find the start and end positions of a fragment in a full-length sequence.
    """
    # check if protein1 and 2 are strings or protein objects
    if isinstance(protein1, str):
        protein1 = next((p for p in proteins if p.name == protein1), None)
    if isinstance(protein2, str):
        protein2 = next((p for p in proteins if p.name == protein2), None)
    
    # Use log file path to construct fasta file path
    directory = os.path.dirname(log_file)
    file_name = os.path.basename(directory)
    # if '_output' in file_name, remove it
    if '_output' in file_name:
        file_name = file_name.replace('_output', '')
    fasta_file = os.path.join(directory, f"{file_name}.fasta")
    # Check if fasta file exists here, if not look in parent directory
    if not os.path.exists(fasta_file):
        # look for fasta file in the parent directory
        parent_directory = os.path.dirname(directory)
        fasta_file = os.path.join(parent_directory, f"{file_name}.fasta")
    if not os.path.exists(fasta_file):
        # create fasta file from pdb file
        # search for any pdb file in the same directory as the log file - rank doesn't matter
        pdb_file = None
        for file in os.listdir(directory):
            if file.endswith('.pdb') or file.endswith('.cif'):
                pdb_file = os.path.join(directory, file)
                break
        if pdb_file:
            # create fasta file from pdb file
            recreate_fasta(pdb_file, fasta_file)

    with open(fasta_file, 'r') as f:
        fasta_data = f.read().splitlines()
        # Get separate sequences by splitting by :
        sequences = fasta_data[1].split(':')
        # if there are not 2 sequences, raise an error
        if len(sequences) != 2:
            raise ValueError(f"Error: Expected 2 sequences in {fasta_file}, found {len(sequences)}")
    
    # Find the fragment sequence positions within the full-length sequence for each protein
    # try, and report if not found
    try:
        protein1_start = protein1.sequence.find(sequences[0])
    except IndexError:
        print(f"Error: Unable to find fragment sequence in FL sequence for {protein1.name}, using fragment from {fasta_file}")
        protein1_start = None
    try:
        protein2_start = protein2.sequence.find(sequences[1])
    except IndexError:
        print(f"Error: Unable to find fragment sequence in FL sequence for {protein2.name}, using fragment from {fasta_file}")
        protein2_start = None
    #protein1_start = protein1.sequence.find(sequences[0])
    #protein2_start = protein2.sequence.find(sequences[1])
    return protein1_start, protein2_start

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

    # Join chains by ':'
    joined = ':'.join(seqs)

    # Determine header from filename
    header = os.path.splitext(os.path.basename(pdb_path))[0]

    # Write FASTA
    with open(fasta_path, 'w') as out:
        out.write(f">{header}\n")
        # write sequence
        out.write(joined + '\n')

def compare_interface_screens(ifaces_a, ifaces_b, min_overlap=1):
    """
    Compare two screens of Interface objects.
    Returns a list of dicts for each protein-pair seen in both screens
    with at least `min_overlap` residues in common on each side.
    """
    def build_map(ifaces):
        # Map interface_id -> {protein1, protein2, res1:set, res2:set}
        m = {}
        for iface in ifaces:
            key = iface.interface_id
            # make sure res1/res2 always correspond to protein1/protein2
            res1, res2 = set(iface.res1), set(iface.res2)
            if key not in m:
                m[key] = {
                    'protein1': iface.protein1,
                    'protein2': iface.protein2,
                    'res1': set(res1),
                    'res2': set(res2)
                }
            else:
                m[key]['res1'].update(res1)
                m[key]['res2'].update(res2)
        return m

    map_a = build_map(ifaces_a)
    map_b = build_map(ifaces_b)

    overlaps = []
    for key, a in map_a.items():
        if key not in map_b:
            continue
        b = map_b[key]
        # compute overlaps on each chain
        ov1 = a['res1'] & b['res1']
        ov2 = a['res2'] & b['res2']
        if len(ov1) >= min_overlap and len(ov2) >= min_overlap:
            overlaps.append({
                'interface_id': key,
                'protein1': a['protein1'].name,
                'protein2': a['protein2'].name,
                'screenA_residues': sorted(a['res1']),
                'screenB_residues': sorted(b['res1']),
                'overlap_residues_1': sorted(ov1),
                'screenA_residues_2': sorted(a['res2']),
                'screenB_residues_2': sorted(b['res2']),
                'overlap_residues_2': sorted(ov2),
                'min_pae_A': min(iface.min_pae for iface in ifaces_a if iface.interface_id==key),
                'min_pae_B': min(iface.min_pae for iface in ifaces_b if iface.interface_id==key),
            })

    return overlaps

def map_interfaces(interface_df, proteins, criteria, alternate_protein_names, screen):
    """
    Map interfaces from a DataFrame onto Protein objects.

    Parameters:
      - interface_df (pd.DataFrame): DataFrame containing interface data.
      - proteins (list): List of Protein objects. These objects will be modified
                         by adding interfaces to them.
      - criteria (dict): Dictionary of criteria to filter the DataFrame.
      - alternate_protein_names (dict): Dictionary mapping alternate protein names
                                        to standard names.

    Returns:
      - A list of created Interface objects.
    """
    # filter interface dataframe based on criteria
    for key, func in criteria.items():
        if key in interface_df.columns:
            interface_df = interface_df[interface_df[key].apply(func)]
    
    # condense interface dataframe if domains and not fl proteins
    if 'Protein1_Domain' in interface_df.columns and 'Protein2_Domain' in interface_df.columns:
        interface_df = condense_interfaces(interface_df, criteria)

    def find_protein(protein_name, proteins, alternate_protein_names, screen):
        """
        Find a protein in the list of proteins by name.
        alternate_protein_names = {'aTub': 'alpha-tubulin', 'bTub': 'beta-tubulin', 'gTub': 'gamma-tubulin', 'sak': 'PLK4', 'polo':'PLK1', 'cnb': 'Centrobin', 'sas-4': 'Sas4', 'sas-6': 'Sas6'}

        """
        protein = next((p for p in proteins if p.name.lower() == protein_name.lower()), None)
        # Adjust for cases where different protein screens use different variants of the same protein
        if protein and (type(protein.screens) == list):
          if (len(protein.screens) >= 1) and (screen not in protein.screens):
                protein = next((p for p in proteins if (p.name.lower() == protein_name.lower()) and (screen in p.screens)), None)
        
        if not protein:
            # Look in alternate_protein_names dict for an alternative name
            protein_name_from_dict = alternate_protein_names.get(protein_name.lower(), protein_name)
            # try looking in alternate names dict, give error if not found
            try:
                protein_name_from_dict = alternate_protein_names[protein_name.lower()]
            except KeyError:
                print(f"Error: Protein '{protein_name}' not found in protein list.")
                return None
            protein = next((p for p in proteins if p.name.lower() == protein_name_from_dict.lower()), None)
            if protein and (type(protein.screens) == list):
                if (len(protein.screens) >= 1) and (screen not in protein.screens):
                    protein = next((p for p in proteins if (p.name.lower() == protein_name_from_dict.lower()) and (screen in p.screens)), None)
        return protein

    interfaces = []
    for index, row in interface_df.iterrows():
        protein1_name = row['Protein1']
        protein2_name = row['Protein2']
        protein1 = find_protein(protein1_name, proteins, alternate_protein_names, screen)
        protein2 = find_protein(protein2_name, proteins, alternate_protein_names, screen)
        
        log_file = row['log_file']
        shift1, shift2 = find_pos_shift(protein1, protein2, log_file)
        if shift1 is None or shift2 is None or shift1 < 0 or shift2 < 0:
            print(f"Error: Unable to find position shift for {protein1.name} and {protein2.name}. - got shift1: {shift1}, shift2: {shift2}")
            print("screen = ", screen, "log file = ", log_file)
            continue
        
        # expand location column to get full range of residues involved in interface
        try:
            locations = expand_loc_range(row['location'])
        except Exception as e:
            try:
                from ast import literal_eval
                locations = expand_loc_range(literal_eval(row['location']))
            except Exception as e:
                print(f"Error: Unable to expand location for {protein1.name} and {protein2.name}. - got {row['location']}")
                continue
        # for every dict in the list of dicts, make interface objects
        if isinstance(locations, list) and len(locations) > 1: # multiple locations, should be separated out
            for i, location in enumerate(locations):
                # Standardize chain keys (e.g., "Chain A", "Chain B", ... -> "Protein1", "Protein2", ...)
                # Collect all keys starting with "Chain "
                chain_keys = [key for key in location if key.startswith("Chain ")]
                if chain_keys:
                    chain_keys.sort() # Sort alphabetically (e.g., "Chain A", "Chain B", "Chain C")
                    
                    # Store values associated with original chain keys and remove them from location
                    chain_values = {}
                    for key in chain_keys:
                        chain_values[key] = location.pop(key)
                    
                    # Add new keys "Protein1", "Protein2", ... using the sorted order
                    for j, old_key in enumerate(chain_keys):
                        new_key = f"Protein{j+1}"
                        location[new_key] = chain_values[old_key]
                # shift every interface res by shift value
                location['Protein1'] = [(x + shift1) for x in location['Protein1']]
                location['Protein2'] = [(x + shift2) for x in location['Protein2']]
                # get the min_pae, avg_pae and size for this interface
                min_pae = row['min_pae'][i]
                avg_pae = row['avg_pae'][i]
                size = row['size'][i]
                # make interface object
                interface = Interface(protein1, protein2, location['Protein1'], location['Protein2'], min_pae, avg_pae, size, screen)
                # add interface to protein objects
                interfaces.append(interface)
                protein1.add_interface(interface)
                protein2.add_interface(interface)
        else: # single location
            if isinstance(locations, list):
                locations = locations[0]
            # Standardize chain keys (e.g., "Chain A", "Chain B", ... -> "Protein1", "Protein2", ...)
            # Collect all keys starting with "Chain "
            chain_keys = [key for key in locations if key.startswith("Chain ")]
            if chain_keys:
                chain_keys.sort() # Sort alphabetically (e.g., "Chain A", "Chain B", "Chain C") 
                # Store values associated with original chain keys and remove them from location
                chain_values = {}
                for key in chain_keys:
                    chain_values[key] = locations.pop(key)
                    
                # Add new keys "Protein1", "Protein2", ... using the sorted order
                for j, old_key in enumerate(chain_keys):
                    new_key = f"Protein{j+1}"
                    locations[new_key] = chain_values[old_key]
            # shift every interface res by shift value
            locations['Protein1'] = [(x + shift1) for x in locations['Protein1']]
            locations['Protein2'] = [(x + shift2) for x in locations['Protein2']]
            # get the min_pae, avg_pae and size for this interface
            min_pae = row['min_pae']
            avg_pae = row['avg_pae']
            size = row['size']
            # make interface object
            interface = Interface(protein1, protein2, locations['Protein1'], locations['Protein2'], min_pae, avg_pae, size, screen)
            # add interface to protein objects
            interfaces.append(interface)
            protein1.add_interface(interface)
            protein2.add_interface(interface)
            
    return interfaces

##########################################################################
def merge_interfaces(proteins, interfaces):
    names = [p.name for p in proteins]
    idx = {n: i for i, n in enumerate(names)}
    iface_data = {}
    # first pass - make dict with all info for each residue pair
    for iface in interfaces:
        n1, n2 = iface.protein1.name, iface.protein2.name
        if n1 in idx and n2 in idx:
            i, j = idx[n1], idx[n2]
            key = tuple(sorted((i, j)))
            entry = iface_data.get(key, {'res1': [], 'res2': [], 'min_pae': np.inf, 'screen': None})
            # collect residues in orientation
            if i < j:
                i_res, j_res = iface.res1, iface.res2
            else:
                i_res, j_res = iface.res2, iface.res1
            entry['res1'].extend(i_res)
            entry['res2'].extend(j_res)
            if n1 == n2: # if same protein, add residues to both sides
                entry['res1'].extend(j_res)
                entry['res2'].extend(i_res)
            entry['min_pae'] = min(entry['min_pae'], iface.min_pae)
            if entry['screen'] is None:
                entry['screen'] = iface.screen
            elif entry['screen'] != iface.screen and iface.screen and iface.screen not in entry['screen']:
                entry['screen'] = entry['screen'] + '_' + iface.screen
            iface_data[key] = entry
    # convert dict back to list of interface objects
    merged_interfaces = []
    for (i, j), data in iface_data.items():
        # create a new Interface object
        merged_interface = Interface(proteins[i], proteins[j], data['res1'], data['res2'], min_pae=data['min_pae'], screen=data['screen'])
        merged_interfaces.append(merged_interface)
    return merged_interfaces

def draw_chord_by_position(proteins, interfaces, chord_colour='protein', protein_colour='alternating', merge=False, smoothness=750, text_direction='perpendicular', text_flip=False, pad=2, figsize=8, cmap_name='viridis', title='Chord Diagram showing predicted interfaces'):
    """
    Draws a chord diagram where:
      - sector widths ∝ sequence length
      - chord endpoints and ribbons positioned by full interface residue spans
      - sequence ends labelled with length
      - depending on colour scheme, either:
        - alternate wedge backgrounds for sectors and ribbons colored by set variable
        - protein-specific wedge colors and gradient ribbons
    
    Parameters:
      - proteins (list): List of Protein objects.
      - interfaces (list): List of Interface objects.
      - chord_colour (str): Variable to color the chords by. Options: 'protein', any variable in interfaces.
      - protein_colour (str): Variable to color proteins by. Options: 'protein' - picks a specific colour for each protein.
                                                                      'alternating' - alternating shades of grey
                                                                       colour name - any colour name from matplotlib
            Note - if chord_colour is 'protein', protein_colour will always be shown as 'protein'
      - smoothness (int): Curve smoothness. Larger values give smoother curves, but slower plotting. Recommended values are between 100 and 1000.
      - text_direction (str): Angle for text rotation - 'perpendicular' (default) or 'parallel'.
      - text_flip (bool): If True, no text will be upside down
      - pad (float): Padding between sectors.
      - figsize (int): Size of the figure.
      - cmap_name (str): Name of the colormap to use.
      - title (str): Plot title.
    """
    # Map names to sequences and lengths
    names = [p.name for p in proteins]
    seq_map = {p.name: p.sequence for p in proteins if p.name in names}
    seq_lens = {n: len(seq_map[n]) for n in names}
    total_len = sum(seq_lens.values())

    if merge:
        interfaces = merge_interfaces(proteins, interfaces)

    # Prepare interface data
    idx = {n: i for i, n in enumerate(names)}
    iface_data = []
    for iface in interfaces:
        n1, n2 = iface.protein1.name, iface.protein2.name
        if n1 in idx and n2 in idx:
            idx1, idx2 = idx[n1], idx[n2]
            # orient residues so order is consistent
            if idx1 < idx2:
                res_1, res_2 = iface.res1, iface.res2
            else:
                res_1, res_2 = iface.res2, iface.res1
            idx1, idx2 = min(idx1, idx2), max(idx1, idx2)
            data = {
                'res1': res_1,
                'res2': res_2,
                'min_pae': iface.min_pae,
                'avg_pae': iface.avg_pae,
                'size': iface.size,
                'interface_id': iface.interface_id,
                'category_id': iface.category_id,
                'screen': iface.screen
            }
            iface_data.append((idx1, idx2, data))

    colour_map = None
    norm = None
    continuous = False
    chord_legend = False
    prot_legend = False

    combined_colours = []
    for map in ['Set2', 'Set3', 'Accent', 'tab20', 'Pastel1', 'Pastel2', 'Dark2', 'Paired', 'tab20b', 'tab20c', 'tab10']:
        plt_map = plt.get_cmap(map)
        map_colours = plt_map(np.linspace(0, 1, plt_map.N))
        combined_colours.extend(map_colours)
    custom_discrete_colormap = np.array(combined_colours)
    
    # assign colours for chords
    if chord_colour != 'protein':
        # check if all items have the colour attribute
        if all (chord_colour in data for _, _, data in iface_data):
            # check if attribute is numerical
            if isinstance(iface_data[0][2][chord_colour], (int, float)):
                # normalise and use cmap
                norm_vals = [e[2][chord_colour] for e in iface_data] or [0,1]
                norm = Normalize(vmin=min(norm_vals), vmax=max(norm_vals))
                colour_map = plt.get_cmap(cmap_name)
                chord_legend = True
                continuous = True # determines legend type
            else:
                # if variable is categorical, use discrete colormap
                unique_vals = set(data[chord_colour] for _, _, data in iface_data)
                if len(unique_vals) <= 4:
                    if len(unique_vals) <= 2:
                        cmap = plt.get_cmap('Set1')
                        colour_map = {val: cmap(i) for i, val in enumerate(unique_vals)}
                    else:
                        cmap1 = plt.get_cmap('Dark2')
                        cmap2 = plt.get_cmap('Set2')
                        cmap = [cmap1(0), cmap1(3), cmap2(5), cmap1(1)]
                        colour_map = {val: cmap[i] for i, val in enumerate(unique_vals)}
                elif protein_colour == 'protein': # use from end of map (so not making false links)
                    colour_map = {val: custom_discrete_colormap[-(i+1)] for i, val in enumerate(unique_vals)}
                else:
                    colour_map = {val: custom_discrete_colormap[i] for i, val in enumerate(unique_vals)}
                chord_legend = True
                continuous = False # determines legend type
        elif is_color_like(chord_colour): # check if chord_colour is a valid colour name
            colour_map = {name: chord_colour for name in names}
        else:
            print(f"Colour variable '{chord_colour}' not found in interfaces. Changing to default.")
            chord_colour = 'protein'
    
    # Assign colors for proteins/wedges
    if protein_colour == 'protein' or (chord_colour == 'protein' and not all(hasattr(p, protein_colour) for p in proteins)):
        prot_colors = {name: custom_discrete_colormap[i] for i, name in enumerate(names)}
    elif protein_colour == 'alternating': # alternate wedge colour
        prot_colors = {name: ('whitesmoke' if i%2==0 else 'grey') for i, name in enumerate(names)}
    elif is_color_like(protein_colour): # check if protein_colour is a valid colour name
        prot_colors = {name: protein_colour for name in names}
    # check if protein colour is attribute of each protein object in list
    elif all(hasattr(p, protein_colour) for p in proteins):
        unique_vals = set([getattr(p, protein_colour) for p in proteins])
        if len(unique_vals) <= 4:
            if len(unique_vals) <= 2:
                cmap = plt.get_cmap('Set1')
                prot_value_colours = {val: cmap(i) for i, val in enumerate(unique_vals)}
            else:
                # custom colormap
                cmap1 = plt.get_cmap('Dark2')
                cmap2 = plt.get_cmap('Set2')
                cmap = [cmap1(0), cmap1(3), cmap2(5), cmap1(1)]
                prot_value_colours = {val: cmap[i] for i, val in enumerate(unique_vals)}
            prot_colors = {protein.name: prot_value_colours[getattr(protein, protein_colour)] for protein in proteins}
        elif chord_legend and not continuous: # if chords already using colourmap, use from end to avoid confusion
            prot_value_colours = {val: custom_discrete_colormap[-(i+1)] for i, val in enumerate(unique_vals)}
            prot_colors = {protein.name: prot_value_colours[getattr(protein, protein_colour)] for protein in proteins}
        else: # if chords are not using discrete colormap, can use from start for proteins
            prot_value_colours = {val: custom_discrete_colormap[i] for i, val in enumerate(unique_vals)}
            prot_colors = {protein.name: prot_value_colours[getattr(protein, protein_colour)] for protein in proteins}
        prot_legend = True
    else: # give error and reset to alternating
        print(f"Invalid protein colour '{protein_colour}'. Using default alternating colours.")
        prot_colors = {name: ('whitesmoke' if i%2==0 else 'grey') for i, name in enumerate(names)}

    # Compute angular span of the circle for each protein
    total_angle = 360.0 - pad * len(names)
    angles = []
    start = 0.0
    for name in names:
        span = seq_lens[name] / total_len * total_angle
        angles.append((360 - start, 360 - (start + span)))
        start += span + pad

    def polar(theta, radius=1.0):
        rad = np.deg2rad(theta)
        return np.cos(rad)*radius, np.sin(rad)*radius

    def quad_bezier(P0, C, P2, t):
        return (1-t)**2 * P0 + 2*(1-t)*t * C + t**2 * P2

    # Draw wedges with protein colors
    def draw_wedge(colour, a0, a1, inner_radius, ax):
        wedge = Wedge((0,0), inner_radius*1.1, a1, a0,
                      width=inner_radius*0.1,
                      facecolor=colour,
                      edgecolor='none', zorder=0)
        ax.add_patch(wedge)
    
    def draw_arc(i, a0, a1, inner_radius, ax):
        thetas = np.linspace(a0, a1, 200)
        xs, ys = zip(*[polar(theta, inner_radius) for theta in thetas])
        ax.plot(xs, ys, color='black', lw=0.8, zorder=2)
    
    def add_protein_labels(i, a0, a1, inner_radius, names, ax, text_angle, text_flip):
        mid = (a0 + a1) / 2
        lx, ly = polar(mid, inner_radius*1.2)
        rot = mid - text_angle
        if text_flip and ((text_angle == 0 and (270 >= mid > 90)) or (text_angle == 90 and mid > 180)):
            rot += 180
        ax.text(lx, ly, names[i], ha='center', va='center', rotation=rot,
                rotation_mode='anchor', zorder=4)
        
    def add_sequence_labels(i, a0, a1, inner_radius, names, ax):
        # ticks at start and end
        for frac in (0, 1):
            angle = a0 + frac*(a1 - a0)
            x_out, y_out = polar(angle, inner_radius + 0.05)
            x_in, y_in = polar(angle, inner_radius)
            ax.plot([x_out, x_in], [y_out, y_in], color='gray', lw=0.8, zorder=2)
        label = str(seq_lens[names[i]])
        tx, ty = polar(a1, inner_radius + 0.06)
        ax.text(tx, ty, label, ha='center', va='center', fontsize=6,
                rotation=a1-90, rotation_mode='anchor', zorder=4)
    
    def draw_chord(i, j, data, angles, seq_lens, names, ax, inner_radius, chord_frac,
                    smoothness, prot_colors, chord_colour, norm, colour_map=None):
        # compute angles for endpoints
        a0_i, a1_i = angles[i]
        a0_j, a1_j = angles[j]
        min_i, max_i = min(data['res1']), max(data['res1'])
        min_j, max_j = min(data['res2']), max(data['res2'])
        theta_i0 = a0_i + (min_i/seq_lens[names[i]])*(a1_i - a0_i)
        theta_i1 = a0_i + (max_i/seq_lens[names[i]])*(a1_i - a0_i)
        theta_j0 = a0_j + (min_j/seq_lens[names[j]])*(a1_j - a0_j)
        theta_j1 = a0_j + (max_j/seq_lens[names[j]])*(a1_j - a0_j)
        # control point
        mid_angle = ( (theta_i0+theta_i1)/2 + (theta_j0+theta_j1)/2 ) / 2
        C = np.array(polar(mid_angle, inner_radius*chord_frac))
        # end points
        P0_top = np.array(polar(theta_i0, inner_radius))
        P2_top = np.array(polar(theta_j1, inner_radius))
        P0_bot = np.array(polar(theta_i1, inner_radius))
        P2_bot = np.array(polar(theta_j0, inner_radius))

        if chord_colour == 'protein':
            # protein end colors
            col_i = np.array(prot_colors[names[i]])
            col_j = np.array(prot_colors[names[j]])
            # sample and draw small quads
            ts = np.linspace(0,1,smoothness)
            top_pts = np.array([quad_bezier(P0_top, C, P2_top, t) for t in ts])
            bot_pts = np.array([quad_bezier(P0_bot, C, P2_bot, t) for t in ts])
            for k in range(len(ts)-1):
                t_mid = 0.5*(ts[k] + ts[k+1])
                color = (1-t_mid)*col_j + t_mid*col_i
                quad = [top_pts[k], top_pts[k+1], bot_pts[k+1], bot_pts[k]]
                poly = Polygon(quad, facecolor=color, edgecolor=None, zorder=1, alpha=0.8)
                ax.add_patch(poly)
        elif is_color_like(chord_colour):
            # if colour is a valid colour name, use it for the whole chord
            ts = np.linspace(0,1,smoothness)
            top_pts = np.array([quad_bezier(P0_top, C, P2_top, t) for t in ts])
            bot_pts = np.array([quad_bezier(P0_bot, C, P2_bot, t) for t in ts])
            verts = np.vstack([top_pts, bot_pts[::-1]])
            poly = Polygon(verts, facecolor=chord_colour,
                        alpha=0.4, edgecolor=None, zorder=1)
            ax.add_patch(poly)
        # if colour is based on a variable, use colour scheme for this
        else:
            ts = np.linspace(0,1,100)
            top_pts = np.array([quad_bezier(P0_top, C, P2_top, t) for t in ts])
            bot_pts = np.array([quad_bezier(P0_bot, C, P2_bot, t) for t in ts])
            verts = np.vstack([top_pts, bot_pts[::-1]])
            if continuous:
                poly = Polygon(verts, facecolor=colour_map(norm(data[chord_colour])),
                            alpha=0.4, edgecolor=None, zorder=1)
            else:
                poly = Polygon(verts, facecolor=colour_map[data[chord_colour]],
                            alpha=0.4, edgecolor=None, zorder=1)
            ax.add_patch(poly)

    # Set parameters for plot

    text_angle = 90 if text_direction == 'parallel' else 0

    inner_radius = 1.0
    chord_frac = 0.4

    figsize_y = figsize * 1.2 if (chord_legend or prot_legend) else figsize
    #figsize_y = figsize

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(figsize, figsize_y), constrained_layout=True)
    ax.set_aspect('equal')
    ax.axis('off')

    for i, (a0, a1) in enumerate(angles):
        draw_wedge(prot_colors[names[i]], a0, a1, inner_radius, ax)
        draw_arc(i, a0, a1, inner_radius, ax)
        add_protein_labels(i, a0, a1, inner_radius, names, ax, text_angle, text_flip)
        add_sequence_labels(i, a0, a1, inner_radius, names, ax)

    # Draw chord for each interface
    for i, j, data in iface_data:
        draw_chord(
            i, j, data,
            angles, seq_lens, names, ax,
            inner_radius, chord_frac,
            smoothness, prot_colors,
            chord_colour, norm, colour_map
        )
 
    # Adding legend
    if chord_legend:
        # if chord_colour variable is continuous
        if continuous:
            sm = ScalarMappable(norm=norm, cmap=colour_map)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, orientation='horizontal', pad=0.05)
            cbar.set_label(chord_colour)
        # if chord_colour variable is categorical
        else:
            # create legend handles
            chord_legend_handles = []
            for val, color in colour_map.items():
                chord_legend_handles.append(plt.Line2D([0], [0], marker='o', color='w', label=val,
                                                  markerfacecolor=color, markersize=10))
            ax.legend(handles=chord_legend_handles, loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.15), fontsize=8)
        if prot_legend:
            # if prot_legend and chord_legend, move prot_legend to top
            prot_legend_handles = []
            for name, color in prot_value_colours.items():
                prot_legend_handles.append(plt.Line2D([0], [0], marker='o', color='w', label=name,
                                                  markerfacecolor=color, markersize=10))
            fig.legend(handles=prot_legend_handles, loc='upper right', bbox_to_anchor=(1.02, 0.85), frameon=False, fontsize=8)
    elif prot_legend: # if only prot_legend, move to bottom
        prot_legend_handles = []
        for name, color in prot_value_colours.items():
            prot_legend_handles.append(plt.Line2D([0], [0], marker='o', color='w', label=name,
                                              markerfacecolor=color, markersize=10))
        ax.legend(handles=prot_legend_handles, loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.15), fontsize=8)

    plt.title(title, y=1.05)
    plt.show()
    return fig


### Making protein objects with FL sequences
protein_csv = "/Users/poppy/Dropbox/protein_ids.csv"
protein_df = pd.read_csv(protein_csv)
# Find missing sequences from UniProt accession IDs
mask = protein_df['sequence'].isna() & protein_df['accession_id'].notna()
protein_df.loc[mask, 'sequence'] = protein_df.loc[mask, 'accession_id'].apply(fetch_sequence)

# Make list of protein objects
proteins = [
    Protein(r.name, r.accession_id, r.sequence, r.category, r.screens)
    for r in protein_df.itertuples(index=False)
]
print("Protein objects created: {}".format(len(proteins)))
print([protein.name for protein in proteins])

### Making interface objects
# read in interface dataframe
interface_csv_screen1 = "/Users/poppy/Dropbox/2022.10.20_Drosophila_Version_1/2022.10.20_Drosophila_Version_1_interface_analysis.csv"
interface_csv_screen2 = "/Users/poppy/Dropbox/all_interface_analysis_2025.05.15.csv"
interface_csv_fl_screen = "/Users/poppy/Dropbox/AF3_FL_Centriole_screen/AF3_FL_Centriole_screen_interface_analysis.csv"
interface_df_screen1 = pd.read_csv(interface_csv_screen1)
interface_df_screen2 = pd.read_csv(interface_csv_screen2)
interface_df_fl_screen = pd.read_csv(interface_csv_fl_screen)

criteria = {
    'rop': lambda x: x >= 3,
    'avg_pae': lambda x: x <= 10,
    'min_pae': lambda x: x <= 3.5,
    'size': lambda x: x >= 5
}
alternate_protein_names = {'alpha-tubulin': 'aTub', 'beta-tubulin': 'bTub', 'gamma-tubulin': 'gTub', 'sak': 'PLK4', 'polo':'PLK1', 'cnb': 'Centrobin', 'sas-4': 'Sas4', 'sas-6': 'Sas6'}

interfaces_screen1 = map_interfaces(interface_df_screen1, proteins, criteria, alternate_protein_names, 'screen1')
interfaces_screen2 = map_interfaces(interface_df_screen2, proteins, criteria, alternate_protein_names, 'screen2')
interfaces_fl_screen = map_interfaces(interface_df_fl_screen, proteins, criteria, alternate_protein_names, 'fl_screen')
print([(interface.protein1.name, interface.protein2.name) for interface in interfaces_fl_screen])
overlaps_1_2 = compare_interface_screens(interfaces_screen1, interfaces_screen2, min_overlap=10)
print(f"Screen1 size: {len(interfaces_screen1)}, Screen2 size: {len(interfaces_screen2)}, Number of overlaps: {len(overlaps_1_2)}")
print(f"Screen FL size: {len(interfaces_fl_screen)}")
print([i['interface_id'] for i in overlaps_1_2])

merged_interfaces1 = merge_interfaces(proteins, interfaces_screen1)
merged_interfaces2 = merge_interfaces(proteins, interfaces_screen2)
merged_interfaces_fl = merge_interfaces(proteins, interfaces_fl_screen)
merged_interfaces = merged_interfaces1 + merged_interfaces2 + merged_interfaces_fl


#######
# Draw chord diagram
start_time = time.time()
names = [p.name for p in proteins]
interfaces = interfaces_screen1 + interfaces_screen2 + interfaces_fl_screen
centriole_proteins = [p for p in proteins if p.category=='Centriole core']
names = ['Ana3', 'Sas4', 'Ana2']
names = ['Ana1', 'ERD10', 'FUS', 'MFP1', 'Ana2', 'Ana3', 'Asl', 'AurA', 'Cep135', 'CHC', 'Cnn', 'Msps', 'PIK3CA', 'Spd2', 'TACC', 'PLK1', 'bamC', 'cpcB', 'PLK4', 'Rcd4', 'Sas4', 'Sas6', ]
centriole_proteins = [p for p in proteins if p.name in names]
core_interfaces = [iface for iface in merged_interfaces
                    if iface.protein1.name in names and iface.protein2.name in names]
draw_chord_by_position(
    centriole_proteins,
    core_interfaces,
    chord_colour='screen',
    text_flip=True
    #protein_colour='category',
    #pad=3,
    #figsize=8
)
end_time = time.time()
print(f"draw_chord_by_position for my screen (centriole) took {end_time - start_time:.2f} seconds")

names = ['Ana3', 'Sas4', 'Ana2']
#centriole_proteins = [p for p in proteins if (p.category=='Centriole core' and p.name not in ['Ctp', 'Centrin', 'CP110', ])]
#centriole_proteins = [p for p in proteins if p.name in names]
centriole_proteins = [p for p in proteins if p.category in ['Centriole core', 'PCM core', 'Kinase'] and p.name not in ['gTub', 'GCP2', 'GCP3', 'GCP4', 'GCP5', 'Grip71', 'Mzt1', 'Mud', 'Gcp6', 'Asp', 'Ctp', 'Centrin', 'Cdk2', 'CycE', 'CycA', 'CycB1', 'CycB3']]
#for screen in ['screen1', 'screen2', 'fl_screen']:
for screen in ['screen1', 'screen2']:
    screen_proteins = []
    for p in centriole_proteins:
        #print(p.name, p.screens)
        if type(p.screens) == list and len(p.screens) > 0:
            if screen in p.screens:
                screen_proteins.append(p)
        else:
            screen_proteins.append(p)
    iface_subset = [i for i in merged_interfaces if i.screen == screen]
    draw_chord_by_position(
        screen_proteins,
        iface_subset,
        #protein_colour='category',
        chord_colour='red' if screen == 'screen1' else 'blue',
        text_flip = True,
        smoothness=100,
        title=screen
    )

draw_chord_by_position(
    centriole_proteins,
    merged_interfaces,
    chord_colour='screen',
    text_flip=True,
    #protein_colour='alternating',
    pad=3,
    figsize=8
)

draw_chord_by_position(
    proteins,
    merged_interfaces2,
    protein_colour='category',
    text_flip=True
)

cep135 = [p for p in proteins if p.name == 'Cep135'][0]
print(cep135.sequence)

# Usage:
centriole_proteins = [p for p in proteins if p.category=='Centriole core']
names = [p.name for p in centriole_proteins]
core_interfaces = [iface for iface in interfaces
                    if iface.protein1.name in names and iface.protein2.name in names]
#test_proteins = [p for p in proteins if p.name in ['Ana2', 'Sas-4', 'Ana3', 'Sas-6', 'Plk4']]
#names = [p.name for p in test_proteins]
#core_interfaces = [iface for iface in interfaces_s2
#                   if iface.protein1.name in names and iface.protein2.name in names]

draw_chord_by_position(
    names,
    centriole_proteins,
    core_interfaces,
    pad=3,
    figsize=8
)

# isaacs screen
centriole_proteins = [p for p in proteins if p.category=='Centriole core']
centriole_names = [p.name for p in centriole_proteins]
draw_chord_by_position(
    centriole_names,
    centriole_proteins,
    interfaces_s1,
    pad=3,
    figsize=8
)
# my screen
draw_chord_by_position(
    centriole_names,
    centriole_proteins,
    interfaces_s2,
    pad=3,
    figsize=8
)
pcm_proteins = [p for p in proteins if p.category=='PCM core']
pcm_names = [p.name for p in pcm_proteins]
draw_chord_by_position(
    pcm_names,
    pcm_proteins,
    interfaces_s1,
    pad=3,
    figsize=8
)
# my screen
draw_chord_by_position(
    pcm_names,
    pcm_proteins,
    interfaces_s2,
    pad=3,
    figsize=8
)



# d1 - d5 of Spd-2+Spd-2
# d1
recreate_fasta('/Users/poppy/Dropbox/2022.10.20_Drosophila_Version_1/Spd-2+Spd-2/Results/Spd-2_D1+Spd-2_D1/Spd-2_D1_Spd-2_D1_unrelaxed_rank_1_model_4.pdb',
                '/Users/poppy/Dropbox/2022.10.20_Drosophila_Version_1/Spd-2+Spd-2/Results/Spd-2_D1+Spd-2_D1.fasta')
# d2
recreate_fasta('/Users/poppy/Dropbox/2022.10.20_Drosophila_Version_1/Spd-2+Spd-2/Results/Spd-2_D2+Spd-2_D2/Spd-2_D2_Spd-2_D2_unrelaxed_rank_1_model_4.pdb',
                '/Users/poppy/Dropbox/2022.10.20_Drosophila_Version_1/Spd-2+Spd-2/Results/Spd-2_D2+Spd-2_D2.fasta')
# d3
recreate_fasta('/Users/poppy/Dropbox/2022.10.20_Drosophila_Version_1/Spd-2+Spd-2/Results/Spd-2_D3+Spd-2_D3/Spd-2_D3_Spd-2_D3_unrelaxed_rank_1_model_1.pdb',
                '/Users/poppy/Dropbox/2022.10.20_Drosophila_Version_1/Spd-2+Spd-2/Results/Spd-2_D3+Spd-2_D3.fasta')
# d4
recreate_fasta('/Users/poppy/Dropbox/2022.10.20_Drosophila_Version_1/Spd-2+Spd-2/Results/Spd-2_D4+Spd-2_D4/Spd-2_D4_Spd-2_D4_unrelaxed_rank_1_model_1.pdb',
                '/Users/poppy/Dropbox/2022.10.20_Drosophila_Version_1/Spd-2+Spd-2/Results/Spd-2_D4+Spd-2_D4.fasta')
# d5
recreate_fasta('/Users/poppy/Dropbox/2022.10.20_Drosophila_Version_1/Spd-2+Spd-2/Results/Spd-2_D5+Spd-2_D5/Spd-2_D5_Spd-2_D5_unrelaxed_rank_1_model_5.pdb',
                '/Users/poppy/Dropbox/2022.10.20_Drosophila_Version_1/Spd-2+Spd-2/Results/Spd-2_D5+Spd-2_D5.fasta')

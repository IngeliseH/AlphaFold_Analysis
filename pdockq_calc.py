"""
pdockq calculation (for AF2 or AF3 output)
Adapted from Elofsson Lab's code
Changes made:
    - Runs as functions rather than through command line
    - Adapted to be compatible with AF3 output (cif files) as well as AF2 output (pdb files)

Functions:
convert_json_to_pkl
parse_atm_record
read_pdb
parse_cif_atm_record
average_plddt_per_residue
calc_pdockq
compute_pdockq
"""
import json
import numpy as np
from collections import defaultdict
import pickle
import gemmi

def convert_json_to_pkl(json_filename):
    """
    Convert JSON file to PKL file

    Parameters:
        - json_filename (str): JSON file to convert
    """
    with open(json_filename, 'r') as json_file:
        data = json.load(json_file)
        if 'plddt' in data:
            plddt_scores = np.array(data['plddt'])
        elif 'atom_plddts' in data:
            plddt_scores = np.array(data['atom_plddts'])
        else:
            raise ValueError("Neither 'plddt' nor 'atom_plddts' found in the JSON file.")
    pkl_filename = str(json_filename).rsplit('.', 1)[0] + '.pkl'
    with open(pkl_filename, 'wb') as pkl_file:
        pickle.dump({'plddt': plddt_scores}, pkl_file)
    return pkl_filename

def parse_atm_record(line):
    """
    Parse an ATOM record from a PDB file
    
    Parameters:
        - line (str): ATOM record from a PDB file
    
    Returns:
        - record (dict): Dictionary containing the parsed ATOM record
    """
    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['atm_alt'] = line[17]
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])
    return record

def read_pdb(pdbfile):
    """
    Read a PDB file and extract the coordinates of the CB atoms

    Parameters:
        - pdbfile (str): Path to the PDB file to read

    Returns:
        - chain_coords (dict): Dictionary containing the coordinates of the CB atoms for each chain
    """
    chain_coords = {}
    with open(pdbfile, 'r') as file:
        for line in file:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            if record['atm_name']=='CB' or (record['atm_name']=='CA' and record['res_name']=='GLY'):
                if record['chain'] in chain_coords:
                    chain_coords[record['chain']].append([record['x'],record['y'],record['z']])
                else:
                    chain_coords[record['chain']] = [[record['x'],record['y'],record['z']]]
    return chain_coords

def parse_cif_atm_record(cif_file):
    """
    Parse an ATOM record from a CIF file
    
    Parameters:
        - cif_file (str): Path to the CIF file to read
    
    Returns:
        - chains (dict): Dictionary containing the coordinates of the CB atoms for each chain
    """
    doc = gemmi.cif.read(str(cif_file))
    block = doc.sole_block()

    #print("Block categories:", block.get_mmcif_category_names())

    atom_site = block.find(['_atom_site.label_asym_id', '_atom_site.label_atom_id', '_atom_site.label_comp_id',
                            '_atom_site.label_seq_id', '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z'])
    if not atom_site:
        print("No _atom_site found")
        return defaultdict(list)

    chains = defaultdict(list)
    residues = defaultdict(list)
    
    for row in atom_site:
        atm_name = row['_atom_site.label_atom_id']
        res_name = row['_atom_site.label_comp_id']
        chain = row['_atom_site.label_asym_id']
        res_id = int(row['_atom_site.label_seq_id'])
        x = float(row['_atom_site.Cartn_x'])
        y = float(row['_atom_site.Cartn_y'])
        z = float(row['_atom_site.Cartn_z'])

        residues[(chain, res_id)].append((atm_name, x, y, z))
        
        if atm_name == 'CB' or (atm_name == 'CA' and res_name == 'GLY'):
            chains[chain].append([x, y, z])
    
    return chains, residues

def average_plddt_per_residue(plddt, residues):
    """
    Calculate the average pLDDT score for each residue

    Parameters:
        - plddt (np.array): Array of pLDDT scores for each atom
        - residues (dict): Dictionary containing the atoms for each residue

    Returns:
        - residue_plddt (np.array): Array of average pLDDT scores for each residue
    """
    residue_plddt = []
    for residue_atoms in residues.values():
        atom_indices = [i for i, (atm_name, x, y, z) in enumerate(residue_atoms)]
        avg_plddt = np.mean(plddt[atom_indices])
        residue_plddt.append(avg_plddt)
    return np.array(residue_plddt)

def calc_pdockq(chain_coords, plddt, t=8):
    """
    Calculate pDockQ and PPV

    Parameters:
        - chain_coords (dict): Dictionary containing the coordinates of the CB atoms for each chain
        - plddt (np.array): Array of pLDDT scores for each residue
        - t (int): Distance threshold for contact calculation

    Returns:
    pdockq (float): pDockQ score
    ppv (float): PPV score
    """
    ch1, ch2 = list(chain_coords.keys())
    coords1, coords2 = np.array(chain_coords[ch1]), np.array(chain_coords[ch2])
    if coords1.shape[0] + coords2.shape[0] != plddt.shape[0]:
        print('Length mismatch with plDDT:', coords1.shape[0] + coords2.shape[0], plddt.shape[0])
    mat = np.append(coords1, coords2, axis=0)
    a_min_b = mat[:, np.newaxis, :] - mat[np.newaxis, :, :]
    dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T
    contact_dists = dists[:len(coords1), len(coords1):]
    contacts = np.argwhere(contact_dists <= t)
    pdockq, ppv = 0, 0
    if contacts.shape[0] > 0:
        avg_if_plddt = np.average(np.concatenate([plddt[np.unique(contacts[:,0])], plddt[np.unique(contacts[:,1]) + len(coords1)]]))
        x = avg_if_plddt * np.log10(len(contacts) + 1)
        pdockq = 0.724 / (1 + np.exp(-0.052 * (x - 152.611))) + 0.018
        
        #PPV
        PPV = np.array([0.98128027, 0.96322524, 0.95333044, 0.9400192 ,
            0.93172991, 0.92420274, 0.91629946, 0.90952562, 0.90043139,
            0.8919553 , 0.88570037, 0.87822061, 0.87116417, 0.86040801,
            0.85453785, 0.84294946, 0.83367787, 0.82238224, 0.81190228,
            0.80223507, 0.78549007, 0.77766077, 0.75941223, 0.74006263,
            0.73044282, 0.71391784, 0.70615739, 0.68635536, 0.66728511,
            0.63555449, 0.55890174])

        pdockq_thresholds = np.array([0.67333079, 0.65666073, 0.63254566, 0.62604391,
            0.60150931, 0.58313803, 0.5647381 , 0.54122438, 0.52314392,
            0.49659878, 0.4774676 , 0.44661346, 0.42628389, 0.39990988,
            0.38479715, 0.3649393 , 0.34526004, 0.3262589 , 0.31475668,
            0.29750023, 0.26673725, 0.24561247, 0.21882689, 0.19651314,
            0.17606258, 0.15398168, 0.13927677, 0.12024131, 0.09996019,
            0.06968505, 0.02946438])
        inds = np.argwhere(pdockq_thresholds>=pdockq)
        if len(inds)>0:
            ppv = PPV[inds[-1]][0]
        else:
            ppv = PPV[0]

    return pdockq, ppv

def compute_pdockq(structure_file, json_file):
    """
    Compute pDockQ and PPV scores for a given structure and JSON file
    
    Parameters:
        - structure_file (str): Path to the structure file (CIF or PDB)
        - json_file (str): Path to the JSON file
    
    Returns:
        - pdockq (float): pDockQ score
        - ppv (float): PPV score
    """
    # Determine file type
    if structure_file.suffix == '.cif':
        chain_coords, residues = parse_cif_atm_record(structure_file)
    elif structure_file.suffix == '.pdb':
        chain_coords = read_pdb(structure_file)
        residues = None  # Not used for PDB files
    else:
        raise ValueError("Unsupported file format. Please provide a CIF or PDB file.")

    pickle_file = convert_json_to_pkl(json_file)
    if len(chain_coords.keys()) < 2:
        print(f'Only one chain in {structure_file}')
        return None, None
    with open(pickle_file, 'rb') as pkl_file:
        pickle_data = pickle.load(pkl_file)
    plddt = pickle_data['plddt']

    # if input is cif, calculate plddt per residue (rather than per atom)
    if residues:
        residue_plddt = average_plddt_per_residue(plddt, residues)
    else:
        residue_plddt = plddt
    return calc_pdockq(chain_coords, residue_plddt)

####################################################################################################
# Example usage
#from analysis_utility import find_rank_001_files

# Example usage for AF2 files
#structure_file, json_file, PAE_png, fasta_file = find_rank_001_files('Sak_Sas6/Sak_D3+Sas6_D1')
#pdockq, ppv = compute_pdockq(structure_file, json_file)
#print(f"Calculated pDockQ = {pdockq}, PPV = {ppv}")

# Example usage for AF3 files
#structure_file, json_file, PAE_png, fasta_file = find_rank_001_files('fold_ana2_flp_sak_fl')
#pdockq, ppv = compute_pdockq(structure_file, json_file)
#print(f"Calculated pDockQ = {pdockq}, PPV = {ppv}")

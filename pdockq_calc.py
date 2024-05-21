"""
pdockq calculation (for AF2 output, currently not working for AF3 output)
includes attempt at adapting for AF3 output

Functions:
convert_json_to_pkl
parse_atm_record
read_pdb
calc_pdockq
compute_pdockq

af3_convert_json_to_pkl
af3_parse_cif_file
af3_calc_pdockq
af3_compute_pdockq
"""
import json
import numpy as np
from collections import defaultdict
import pickle

def convert_json_to_pkl(json_filename):
    with open(json_filename, 'r') as json_file:
        data = json.load(json_file)
        plddt_scores = np.array(data['plddt'])
    pkl_filename = json_filename.rsplit('.', 1)[0] + '.pkl'
    with open(pkl_filename, 'wb') as pkl_file:
        pickle.dump({'plddt': plddt_scores}, pkl_file)
    return pkl_filename

def parse_atm_record(line):
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

def calc_pdockq(chain_coords, plddt, t=8):
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
        ###################
        # Placeholder for actual PPV calculation logic
        ppv = None
    return pdockq, ppv

def compute_pdockq(pdb_file, json_file):
    pickle_file = convert_json_to_pkl(json_file)
    chain_coords = read_pdb(pdb_file)
    if len(chain_coords.keys()) < 2:
        print('Only one chain in pdbfile', pdb_file)
        return None, None
    with open(pickle_file, 'rb') as pkl_file:
        pickle_data = pickle.load(pkl_file)
    plddt = pickle_data['plddt']
    return calc_pdockq(chain_coords, plddt)


from analysis_utility import find_rank_001_files
structure_file, json_file, PAE_png, fasta_file = find_rank_001_files('Sak_Sas6/Sak_D3+Sas6_D1')
structure_file, json_file, PAE_png, fasta_file = find_rank_001_files('fold_ana2p_cterm_dimer_sas6_dimer', af3=True)


if structure_file is not None  and json_file is not None:  # Check if files were found
    pdockq, ppv = compute_pdockq(structure_file, json_file)
    print(f"Calculated pDockQ = {pdockq}, PPV = {ppv}")





############################################################################################################
#attempt to calculate pdockq from AF3 output (cif file instead of pdb and atom_plddts in json instead of plddt)
from analysis_utility import find_rank_001_files
import gemmi

def af3_convert_json_to_pkl(json_filename):
    with open(json_filename, 'r') as json_file:
        data = json.load(json_file)
        # Extract 'atom_plddts' instead of 'plddt'
        plddt_scores = np.array(data['atom_plddts'])
    pkl_filename = json_filename.rsplit('.', 1)[0] + '.pkl'
    with open(pkl_filename, 'wb') as pkl_file:
        pickle.dump({'plddt': plddt_scores}, pkl_file)
    return pkl_filename

def af3_parse_cif_atm_record(cif_file):
    doc = gemmi.cif.read(cif_file)
    block = doc.sole_block()

    print("Block categories:", block.get_mmcif_category_names())

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

def af3_average_plddt_per_residue(plddt, residues):
    residue_plddt = []
    for residue_atoms in residues.values():
        atom_indices = [i for i, (atm_name, x, y, z) in enumerate(residue_atoms)]
        avg_plddt = np.mean(plddt[atom_indices])
        residue_plddt.append(avg_plddt)
    return np.array(residue_plddt)

def af3_calc_pdockq(chain_coords, plddt, t=8):
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
        ###################
        # Placeholder for actual PPV calculation logic
        ppv = None
    return pdockq, ppv

def af3_compute_pdockq(cif_file, json_file):
    pickle_file = af3_convert_json_to_pkl(json_file)
    chain_coords, residues = af3_parse_cif_atm_record(cif_file)
    print('chain coords from cif parse', chain_coords)
    if len(chain_coords.keys()) < 2:
        print('Only one chain in cif file', cif_file)
        return None, None
    with open(pickle_file, 'rb') as pkl_file:
        pickle_data = pickle.load(pkl_file)
    plddt = pickle_data['plddt']

    residue_plddt = af3_average_plddt_per_residue(plddt, residues)
    return af3_calc_pdockq(chain_coords, residue_plddt)

# Example usage
structure_file, json_file, PAE_png, fasta_file = find_rank_001_files('fold_ana2_flp_sak_fl')
print(structure_file, json_file)

if structure_file is not None and json_file is not None:  # Check if files were found
    pdockq, ppv = af3_compute_pdockq(structure_file, json_file)
    print(f"Calculated pDockQ = {pdockq}, PPV = {ppv}")

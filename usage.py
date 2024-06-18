import os
from analysis_utility import find_rank_001_files, parse_structure_file, determine_chain_lengths, extract_pae
from phosphorylation_handling import correct_cif_pae, output_chimerax_readable_json
from process_pae import visualize_pae_matrix, find_min_pae
from repeatability_from_pdb import measure_repeatability
from iptm_visualisation import create_iptm_matrix, visualize_iptm_matrix, extract_iptm
from pdockq_calc import compute_pdockq

# Example usage:
# AF2 files
af2_folder_path = 'Sak_Sas6/Sak_D3+Sas6_D1'
# Process folder to find relevant files and read in data
pdb_file, json_file, log_file, _, _ = find_rank_001_files(af2_folder_path)
structure_model = parse_structure_file(pdb_file)
pae_matrix = extract_pae(json_file)
chain_lengths = determine_chain_lengths(structure_model)

# Find ipTM score
iptm = extract_iptm(log_file)
print(f"ipTM = {iptm}")

# Plot pae
visualize_pae_matrix(pae_matrix, chain_lengths)
# Find minimum interface pae(s)
min_pae = find_min_pae(structure_model, pae_matrix)
print(f"Minimum interface PAE = {min_pae}")

# Calculate pdockq
pdockq, ppv = compute_pdockq(pdb_file, json_file)
print(f"Calculated pDockQ = {pdockq}, PPV = {ppv}")

# Measure repeatability
num_consistent, level_consistent = measure_repeatability(af2_folder_path, distance_cutoff=5, pae_cutoff=10, all_atom=True)

# For use only on protein pair folder containing multiple domain pair folders
# Create and visualize iptm matrix (AF2 only, as not planning on running multiple domain pairs in AF3)
iptm_matrix = create_iptm_matrix(af2_folder_path)
iptm_png_file_path = os.path.join(af2_folder_path, 'iptm_matrix.png')
visualize_iptm_matrix(iptm_matrix, iptm_png_file_path)


# AF3 files
af3_folder_path = 'fold_ana2p_cterm_dimer_sas6_dimer'
af3_folder_path = 'fold_ana2p_cterm_monomer_sas6_dimer'
af3_folder_path = 'fold_ana2_flp_sak_fl'
af3_folder_path = 'fold_ana2flp_all_dimer_sas6_dimer'
# Process folder to find relevant files and read in data
cif_file, json_file, _, _, _ = find_rank_001_files(af3_folder_path)
structure_model = parse_structure_file(cif_file, is_pdb=False)
initial_pae_matrix = extract_pae(json_file)
pae_matrix = correct_cif_pae(structure_model, initial_pae_matrix)
chain_lengths = determine_chain_lengths(structure_model)

# Find ipTM score
iptm = extract_iptm(log_file)
print(f"ipTM = {iptm}")

# Plot pae
visualize_pae_matrix(pae_matrix, chain_lengths, pae_cutoff=30)
# Save pae data in a format readable by ChimeraX
output_chimerax_readable_json(pae_matrix, save_location=os.path.join(af3_folder_path, 'pae_data.json'))
# Find minimum interface pae(s)
min_pae = find_min_pae(structure_model, pae_matrix)
print(f"Minimum interface PAE = {min_pae}")

# Calculate pdockq
pdockq, ppv = compute_pdockq(cif_file, json_file)
print(f"Calculated pDockQ = {pdockq}, PPV = {ppv}")

# Measure repeatability
num_consistent, level_consistent = measure_repeatability(af3_folder_path, distance_cutoff=5, pae_cutoff=10, all_atom=True)

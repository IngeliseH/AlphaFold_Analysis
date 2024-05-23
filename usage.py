from analysis_utility import parse_cif, parse_structure_file, extract_pae, determine_cif_chain_lengths, extract_fasta_protein_lengths, find_rank_001_files
from phosphorylation_handling import map_res_to_pos, condense_pae
from visualize_pae import visualize_pae_matrix

# Example usage:
# AF2 files
pdb_file, json_file, PAE_png, fasta_file = find_rank_001_files('Sak_Sas6/Sak_D3+Sas6_D1')

pae_data = extract_pae(json_file)

chain_lengths = extract_fasta_protein_lengths(fasta_file)
visualize_pae_matrix(pae_data, chain_lengths=chain_lengths)

# AF3 files
cif_file, json_file, PAE_png, fasta_file = find_rank_001_files('fold_ana2p_cterm_dimer_sas6_dimer')

chain_residue_map = parse_cif(cif_file)
res_to_pos = map_res_to_pos(chain_residue_map)
pae_matrix = extract_pae(json_file)
condensed_pae_matrix = condense_pae(pae_matrix, res_to_pos)

chain_lengths = determine_cif_chain_lengths(chain_residue_map)
visualize_pae_matrix(condensed_pae_matrix, chain_lengths=chain_lengths)


cif_file, json_file, PAE_png, fasta_file = find_rank_001_files('fold_ana2_flp_sak_fl')
print(cif_file)

chain_residue_map = parse_structure_file(cif_file)
print(chain_residue_map)
pae_matrix = extract_pae(json_file)

chain_lengths = determine_cif_chain_lengths(chain_residue_map)
visualize_pae_matrix(pae_matrix, chain_lengths=chain_lengths)

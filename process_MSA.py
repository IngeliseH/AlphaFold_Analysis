"""
Module to handle MSAs pre-converted from m3a to fasta form, calculating BLOSUM
similarity scores for each residue of a 2 protein prediction, and modifying PDB
file b factor column to record this. Calculates sequence similarity to decide
which BLOSUM matrix to use.

Functions:
    - parse_structure_file(PDB_file)
    - determine_chain_lengths(protein_model)
    - remove_gap_columns(aligned_sequences, query_seq)
    - calculate_percent_similarity(alignment)
    - calculate_residue_similarity_to_query(alignment, query_seq, blosum)
    - modify_bfactors(protein_model, residue_similarity_scores_list, protein_lengths, output_name)
    - process_alignment(MSA_file, PDB_file)
"""
from Bio import AlignIO
from Bio.PDB import PDBIO, Model
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
from pathlib import Path
import numpy as np
from analysis_utility import parse_structure_file, determine_chain_lengths

def remove_gap_columns(aligned_sequences, query_seq):
    """
    Removes columns where the query sequence has a gap from all sequences in the alignment.

    Parameters:
        - alignment (MultipleSeqAlignment): The full alignment.
        - query_seq (Seq): The query sequence in the alignment.

    Returns:
        - new_alignment (MultipleSeqAlignment): The alignment with gap columns removed.
    """
    # Convert alignment to a list of lists
    alignment_array = [list(record.seq) for record in aligned_sequences]
    alignment_length = len(alignment_array[0])
    
    # Identify columns that are not gaps in the query sequence
    non_gap_positions = [i for i in range(alignment_length) if query_seq[i] != '-']
    
    # Filter out gap columns
    gapless_alignment = []
    for record in aligned_sequences:
        seq = ''.join([record.seq[i] for i in non_gap_positions])
        gapless_alignment.append(SeqRecord(Seq(seq), id=record.id))
    
    return MultipleSeqAlignment(gapless_alignment)

def calculate_percent_similarity(alignment):
    """
    Calculates the percentage similarity of each sequence in the alignment to the query sequence.

    Parameters:
        - alignment (MultipleSeqAlignment): The alignment with sequences.
    
    Returns:
        - similarities (list): List of similarity scores for each sequence
    """
    query_seq = alignment[0].seq
    similarities = []
    for record in alignment[1:]:
        seq = record.seq
        matches = sum(1 for a, b in zip(query_seq, seq) if a == b and a != '-')
        total = sum(1 for a in query_seq if a != '-')
        similarity = (matches / total) * 100 if total > 0 else 0
        similarities.append(similarity)
    return similarities

def calculate_residue_similarity_to_query(alignment, query_seq, blosum):
    """
    Calculates average similarity scores for each residue in the query sequence across all aligned sequences.

    Parameters:
        - alignment (MultipleSeqAlignment): The alignment with sequences.
        - query_seq (Seq): The query sequence.
        - blosum (dict): Substitution matrix (e.g., blosum62).

    Returns:
        - residue_similarity_scores (list): List of average similarity scores for each residue in the query sequence.
    """
    seq_length = len(query_seq)
    scores = [0.0] * seq_length
    #num_sequences = len(alignment) - 1  # Exclude the query sequence

    for i in range(seq_length):
        aa_query = query_seq[i]
        if aa_query == '-':
            continue  # Skip gap positions in the query sequence

        total_score = 0
        valid_counts = 0
        for record in alignment[1:]:  # Exclude the query sequence
            aa_subject = record.seq[i]
            if aa_subject == '-':
                continue  # Skip gaps in subject sequences

            # Get the substitution score; default to 0 if not found
            score = blosum.get((aa_query, aa_subject), blosum.get((aa_subject, aa_query), 0))
            total_score += score
            valid_counts += 1

        if valid_counts > 0:
            scores[i] = total_score / valid_counts
        else:
            scores[i] = 0.0  # No valid comparisons at this position

    return scores

def modify_bfactors(protein_model, residue_similarity_scores_list, output_name):
    """
    Modifies the B-factor column of the PDB file to store the residue similarity scores.

    Parameters:
        - protein_model (Model): The protein structure model.
        - residue_similarity_scores_list (list): List of similarity scores for each residue.
        - output_name (str or Path): The name of the output PDB file.
    Returns:
        None, but saves the modified PDB file.
    """
    # Flatten the residue similarity scores into a single list
    all_scores = []
    for scores in residue_similarity_scores_list:
        all_scores.extend(scores)
    
    residue_counter = 0
    for chain in protein_model.get_chains():
        for residue in chain.get_residues():
            for atom in residue.get_atoms():
                if residue_counter < len(all_scores):
                    atom.bfactor = all_scores[residue_counter]
                else:
                    atom.bfactor = 0.0  # Default B-factor if out of scores
            residue_counter += 1
    
    # Save the modified structure
    io = PDBIO()
    io.set_structure(protein_model)

    # Ensure output_name is a string
    if isinstance(output_name, Path):
        output_name = str(output_name)
    io.save(output_name)

def process_alignment(MSA_file, structure_input):
    """
    Processes an alignment file and calculates residue conservation scores.
    The structure input can be either a file path or a protein model object.
    
    Parameters:
        - MSA_file (str): The path to the alignment file in fasta format.
        - structure_input (str, Path  or Bio.PDB.Model): Either a file path to the PDB file or a protein model object.

    Returns:
        - residue_similarity_scores_list (list): List of conservation scores for each protein.
        - start_positions (list): Starting positions of each protein in the concatenated sequence.
    """
    # Read the aligned sequences
    aligned_sequences = AlignIO.read(MSA_file, "fasta")
    
    # Query sequence is the first one in the alignment
    query_seq_record = aligned_sequences[0]
    query_seq = query_seq_record.seq
    if isinstance(structure_input, str) or isinstance(structure_input, Path):
        # Treat structure_input as a file path and parse the structure file
        protein_model = parse_structure_file(structure_input)
    elif isinstance(structure_input, Model.Model):
        # Treat structure_input as an already loaded protein model object
        protein_model = structure_input
    else:
        raise TypeError("structure_input must be either a file path (str or Path) or a protein model object (Bio.PDB.Model.Model).")


    # Remove gap columns based on the query sequence
    gapless_alignment = remove_gap_columns(aligned_sequences, query_seq)

    # Determine chain lengths and starting positions
    protein_lengths = determine_chain_lengths(protein_model)  # Returns a list of lengths
    num_proteins = len(protein_lengths)
    start_positions = [0]
    for length in protein_lengths:
        start_positions.append(start_positions[-1] + length)

    # Split the alignment into separate proteins
    protein_seqs = [[] for _ in range(num_proteins)]
    query_seqs = [None for _ in range(num_proteins)]
    first_sequence = True

    for seq_record in gapless_alignment:
        seq = seq_record.seq

        regions = []
        for i in range(num_proteins):
            start = start_positions[i]
            end = start_positions[i + 1]
            region = seq[start:end]
            regions.append(region)

        if first_sequence:
            query_seqs = regions
            first_sequence = False

        for i, region in enumerate(regions):
            if set(region) - set('-'):
                protein_seqs[i].append(SeqRecord(Seq(str(region)), id=seq_record.id))

    # Calculate conservation scores for each protein
    residue_similarity_scores_list = []
    for i in range(num_proteins):
        alignment = MultipleSeqAlignment(protein_seqs[i])
        avg_similarity = np.mean(calculate_percent_similarity(alignment))
        
        # Select appropriate BLOSUM matrix
        if avg_similarity >= 80:
            blosum = substitution_matrices.load("BLOSUM80")
        elif avg_similarity >= 50:
            blosum = substitution_matrices.load("BLOSUM62")
        else:
            blosum = substitution_matrices.load("BLOSUM45")

        # Calculate residue similarity scores
        scores = calculate_residue_similarity_to_query(alignment, query_seqs[i], blosum)
        residue_similarity_scores_list.append(scores)

    return residue_similarity_scores_list, start_positions

def calculate_interface_conservation_score(residue_similarity_scores_list, interface_residues, abs_res_reverse_lookup, start_positions):
    """
    Calculates the overall interface conservation score by taking the minimum conservation
    score for each interface residue pair and averaging these minimums.

    Parameters:
        - residue_similarity_scores_list (list): List of conservation scores for each protein.
        - interface_residues (list): List of residue pairs (tuples) representing interface contacts.
                                     Each tuple contains absolute residue indices (starting from 1).
        - abs_res_reverse_lookup (dict): Mapping of absolute residue indices to (chain_id, res_num).
        - start_positions (list): Starting positions of each protein in the concatenated sequence.

    Returns:
        - overall_conservation_score (float): The average of minimum conservation scores for interface residue pairs.
    """
    # Create a mapping from (chain_id, res_num) to conservation score index
    residue_index_mapping = {}
    pos_counter = 0
    for chain_idx, scores in enumerate(residue_similarity_scores_list):
        chain_length = len(scores)
        for res_idx in range(chain_length):
            abs_res_idx = start_positions[chain_idx] + res_idx + 1  # +1 to adjust to 1-based indexing
            chain_id, res_num = abs_res_reverse_lookup.get(abs_res_idx, (None, None))
            if chain_id is not None:
                residue_index_mapping[(chain_id, res_num)] = (chain_idx, res_idx)
        pos_counter += chain_length

    # Collect minimum conservation scores for each interface residue pair
    min_conservation_scores = []
    for res1_abs, res2_abs in interface_residues:
        # Map absolute residue indices back to (chain_id, res_num)
        res1_info = abs_res_reverse_lookup.get(res1_abs)
        res2_info = abs_res_reverse_lookup.get(res2_abs)

        if res1_info and res2_info:
            res1_chain_id, res1_res_num = res1_info
            res2_chain_id, res2_res_num = res2_info

            res1_mapping = residue_index_mapping.get((res1_chain_id, res1_res_num))
            res2_mapping = residue_index_mapping.get((res2_chain_id, res2_res_num))

            if res1_mapping and res2_mapping:
                res1_chain_idx, res1_idx = res1_mapping
                res2_chain_idx, res2_idx = res2_mapping

                # Get conservation scores
                score1 = residue_similarity_scores_list[res1_chain_idx][res1_idx]
                score2 = residue_similarity_scores_list[res2_chain_idx][res2_idx]

                # Take the minimum score
                min_score = min(score1, score2)
                min_conservation_scores.append(min_score)

    # Calculate the overall conservation score
    if min_conservation_scores:
        overall_conservation_score = np.mean(min_conservation_scores)
    else:
        overall_conservation_score = None  # No interface residues found

    return overall_conservation_score

# example usage
#process_alignment("data/MSA_analysis/GCP5_F1_Mzt1_F1.fasta", "data/MSA_analysis/GCP5_F1_Mzt1_F1_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_000.pdb")

#from analysis_utility import map_chains_and_residues
#from repeatability_from_pdb import find_confident_interface_residues
#structure_model = parse_structure_file("data/MSA_analysis/GCP5_F1_Mzt1_F1_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_000.pdb")
#residue_similarity_scores_list, start_positions = process_alignment("data/MSA_analysis/GCP5_F1_Mzt1_F1.fasta", structure_model)
#chain_residue_map = map_chains_and_residues(structure_model)
# Create mappings between absolute residue indices and (chain_id, res_num)
#abs_res_lookup_dict = {(entry[0], entry[1]): entry[3] for entry in chain_residue_map}
#abs_res_reverse_lookup = {entry[3]: (entry[0], entry[1]) for entry in chain_residue_map}
#residue_pairs = find_confident_interface_residues(structure_model, "data/MSA_analysis/GCP5_F1_Mzt1_F1_scores_rank_001_alphafold2_multimer_v3_model_4_seed_000.json", distance_cutoff=5, pae_cutoff=15, abs_res_lookup_dict=abs_res_lookup_dict, is_pdb=True, all_atom=True)
#interface_size = len(residue_pairs)
# Compute conservation score using the new function
#conservation_score = calculate_interface_conservation_score(residue_similarity_scores_list, residue_pairs, abs_res_reverse_lookup, start_positions)

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
from Bio.PDB import PDBParser, PDBIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
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
    num_sequences = len(alignment) - 1  # Exclude the query sequence

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
        - output_name (str): The name of the output PDB file.
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
    io.save(output_name)

def process_alignment(MSA_file, PDB_file):
    """
    Takes an alignment file in fasta form, using the first sequence as the query
    sequence. Calculates average percent similarity of sequences to query for each protein, 
    then calculates similarity scores for each residue in the query sequence across all
    aligned sequences, using the best matching BLOSUM substitution matrix. Modifies the
    B-factor column of the PDB file to store the residue conservation scores.

    Parameters:
        - MSA_file (str): The path to the alignment file in fasta format.
        - PDB_file (str): The path to the PDB file

    Returns:
        None, but saves the modified PDB file with conservation scores in the B-factor column.
    """
    aligned_sequences = AlignIO.read(MSA_file, "fasta")
    
    # Query sequence is the first one in the alignment
    query_seq_record = aligned_sequences[0]
    query_seq = query_seq_record.seq

    # Parse the structure file
    protein_model = parse_structure_file(PDB_file)

    gapless_alignment = remove_gap_columns(aligned_sequences, query_seq)

    protein_lengths = determine_chain_lengths(protein_model)  # Returns a list of lengths

    num_proteins = len(protein_lengths)
    protein_seqs = [[] for _ in range(num_proteins)]
    query_seqs = [None for _ in range(num_proteins)]
    first_sequence = True

    # Calculate start and end positions for each protein region
    start_positions = [0]
    for length in protein_lengths:
        # Adjust the index for zero-based indexing in Python
        start_positions.append(start_positions[-1] + length)

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

    # Assign B-factors (similarity scores)
    output_name = PDB_file.split(".")[0] + "_conservation.pdb"
    modify_bfactors(protein_model, residue_similarity_scores_list, output_name)

# Example usage
process_alignment("data/MSA_analysis/Ana2_D2_Sak_D1.fasta", "data/MSA_analysis/Ana2_D2_Sak_D1_unrelaxed_rank_1_model_2.pdb")

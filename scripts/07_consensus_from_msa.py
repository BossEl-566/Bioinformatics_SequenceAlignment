from Bio import AlignIO
from collections import Counter
import os

MSA_FILE = "results/msa/muscle_alignment.fasta"
OUTPUT_FILE = "results/msa/consensus_sequence.txt"

os.makedirs("results/msa", exist_ok=True)

alignment = AlignIO.read(MSA_FILE, "fasta")
alignment_length = alignment.get_alignment_length()

def consensus_from_alignment(aln, threshold=0.7, ambiguous="X"):
    """
    Build a simple consensus sequence by column voting.
    Ignores gaps when choosing the most common residue.
    If no residue reaches threshold, uses ambiguous character.
    """
    consensus_chars = []
    nseq = len(aln)

    for col_idx in range(alignment_length):
        column = aln[:, col_idx]  # string of residues at this column
        residues = [aa for aa in column if aa != "-"]

        if not residues:
            consensus_chars.append("-")
            continue

        counts = Counter(residues)
        top_residue, top_count = counts.most_common(1)[0]

        # threshold based on non-gap residues (more biologically sensible)
        if (top_count / len(residues)) >= threshold:
            consensus_chars.append(top_residue)
        else:
            consensus_chars.append(ambiguous)

    return "".join(consensus_chars)

consensus_seq = consensus_from_alignment(alignment, threshold=0.7, ambiguous="X")

with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
    f.write("Consensus sequence (from MUSCLE MSA)\n")
    f.write("Threshold: 0.7\n")
    f.write(f"Number of sequences: {len(alignment)}\n")
    f.write(f"Alignment length: {alignment_length}\n\n")
    f.write(consensus_seq + "\n")

print("Consensus generation complete.")
print(f"Alignment length: {alignment_length}")
print("Consensus sequence:")
print(consensus_seq)
print(f"\nSaved: {OUTPUT_FILE}")
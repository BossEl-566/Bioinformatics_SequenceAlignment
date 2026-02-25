from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
import os
import csv

INPUT_FASTA = "data/cleaned/hemoglobin_beta_cleaned.fasta"
PAIRWISE_DIR = "results/pairwise"
TABLES_DIR = "results/tables"

os.makedirs(PAIRWISE_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

# Load cleaned sequences
records = list(SeqIO.parse(INPUT_FASTA, "fasta"))

if len(records) < 2:
    raise ValueError("Need at least 2 sequences for pairwise alignment.")

# For now, use first two sequences (Human and Chimp if your order is unchanged)
rec1 = records[0]
rec2 = records[1]

seq1 = str(rec1.seq)
seq2 = str(rec2.seq)

# Load BLOSUM62 substitution matrix (for proteins)
blosum62 = substitution_matrices.load("BLOSUM62")

def run_alignment(mode, open_gap=-10, extend_gap=-0.5):
    aligner = Align.PairwiseAligner()
    aligner.mode = mode  # "global" or "local"
    aligner.substitution_matrix = blosum62
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap

    alignments = aligner.align(seq1, seq2)
    best = alignments[0]
    return best

# Run global and local alignments
global_alignment = run_alignment("global")
local_alignment = run_alignment("local")

# Save alignment outputs as text files
with open(os.path.join(PAIRWISE_DIR, "global_alignment.txt"), "w", encoding="utf-8") as f:
    f.write(f"Sequence 1: {rec1.id}\n")
    f.write(f"Sequence 2: {rec2.id}\n")
    f.write(f"Mode: global\n")
    f.write(f"Score: {global_alignment.score}\n\n")
    f.write(str(global_alignment))

with open(os.path.join(PAIRWISE_DIR, "local_alignment.txt"), "w", encoding="utf-8") as f:
    f.write(f"Sequence 1: {rec1.id}\n")
    f.write(f"Sequence 2: {rec2.id}\n")
    f.write(f"Mode: local\n")
    f.write(f"Score: {local_alignment.score}\n\n")
    f.write(str(local_alignment))

# Save summary scores to CSV
csv_path = os.path.join(TABLES_DIR, "pairwise_scores.csv")
with open(csv_path, "w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["sequence_1", "sequence_2", "mode", "gap_open", "gap_extend", "score"])
    writer.writerow([rec1.id, rec2.id, "global", -10, -0.5, global_alignment.score])
    writer.writerow([rec1.id, rec2.id, "local", -10, -0.5, local_alignment.score])

print("Pairwise alignment complete.")
print(f"Sequence 1: {rec1.id}")
print(f"Sequence 2: {rec2.id}")
print(f"Global score: {global_alignment.score}")
print(f"Local score: {local_alignment.score}")
print(f"Saved: {os.path.join(PAIRWISE_DIR, 'global_alignment.txt')}")
print(f"Saved: {os.path.join(PAIRWISE_DIR, 'local_alignment.txt')}")
print(f"Saved: {csv_path}")
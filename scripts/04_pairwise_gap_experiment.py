from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
import os
import csv

INPUT_FASTA = "data/cleaned/hemoglobin_beta_cleaned.fasta"
TABLES_DIR = "results/tables"
PAIRWISE_DIR = "results/pairwise"

os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(PAIRWISE_DIR, exist_ok=True)

records = list(SeqIO.parse(INPUT_FASTA, "fasta"))

# Find human and zebrafish
human = None
zebrafish = None

for r in records:
    rid = r.id.upper()
    if "HBB_HUMAN" in rid:
        human = r
    if "HBB1_DANRE" in rid or "DANRE" in rid:
        zebrafish = r

if human is None or zebrafish is None:
    raise ValueError("Could not find Human or Zebrafish sequences.")

seq1 = str(human.seq)
seq2 = str(zebrafish.seq)

blosum62 = substitution_matrices.load("BLOSUM62")

gap_settings = [
    (-8, -0.2),
    (-10, -0.5),
    (-12, -1.0),
]

def run_alignment(mode, open_gap, extend_gap):
    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = blosum62
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap
    return aligner.align(seq1, seq2)[0]

rows = []

for open_gap, extend_gap in gap_settings:
    for mode in ["global", "local"]:
        best = run_alignment(mode, open_gap, extend_gap)
        rows.append([
            human.id,
            zebrafish.id,
            mode,
            open_gap,
            extend_gap,
            best.score
        ])

csv_path = os.path.join(TABLES_DIR, "pairwise_gap_experiments_human_vs_zebrafish.csv")
with open(csv_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["sequence_1", "sequence_2", "mode", "gap_open", "gap_extend", "score"])
    writer.writerows(rows)

print("Gap penalty experiment complete.")
print(f"Saved: {csv_path}")
print("\nResults:")
for row in rows:
    print(row)
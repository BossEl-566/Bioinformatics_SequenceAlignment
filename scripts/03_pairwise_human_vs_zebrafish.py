from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
import os

INPUT_FASTA = "data/cleaned/hemoglobin_beta_cleaned.fasta"
PAIRWISE_DIR = "results/pairwise"
os.makedirs(PAIRWISE_DIR, exist_ok=True)

records = list(SeqIO.parse(INPUT_FASTA, "fasta"))

# Find human and zebrafish by ID keywords
human = None
zebrafish = None

for r in records:
    rid = r.id.upper()
    if "HBB_HUMAN" in rid:
        human = r
    if "DANRE" in rid or "ZEBRAFISH" in rid or "HBB1_DANRE" in rid:
        zebrafish = r

if human is None:
    raise ValueError("Human HBB sequence not found.")
if zebrafish is None:
    raise ValueError("Zebrafish HBB sequence not found.")

seq1 = str(human.seq)
seq2 = str(zebrafish.seq)

blosum62 = substitution_matrices.load("BLOSUM62")

def run_alignment(mode, open_gap=-10, extend_gap=-0.5):
    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = blosum62
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap
    return aligner.align(seq1, seq2)[0]

global_alignment = run_alignment("global")
local_alignment = run_alignment("local")

with open(os.path.join(PAIRWISE_DIR, "human_vs_zebrafish_global.txt"), "w", encoding="utf-8") as f:
    f.write(f"Sequence 1: {human.id}\n")
    f.write(f"Sequence 2: {zebrafish.id}\n")
    f.write("Mode: global\n")
    f.write(f"Score: {global_alignment.score}\n\n")
    f.write(str(global_alignment))

with open(os.path.join(PAIRWISE_DIR, "human_vs_zebrafish_local.txt"), "w", encoding="utf-8") as f:
    f.write(f"Sequence 1: {human.id}\n")
    f.write(f"Sequence 2: {zebrafish.id}\n")
    f.write("Mode: local\n")
    f.write(f"Score: {local_alignment.score}\n\n")
    f.write(str(local_alignment))

print("Human vs Zebrafish alignment complete.")
print("Sequence 1:", human.id)
print("Sequence 2:", zebrafish.id)
print("Global score:", global_alignment.score)
print("Local score:", local_alignment.score)
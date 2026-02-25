from Bio import AlignIO
import os
import csv

MSA_FILE = "results/msa/muscle_alignment.fasta"
TABLES_DIR = "results/tables"

os.makedirs(TABLES_DIR, exist_ok=True)

# Read MUSCLE output as FASTA alignment
alignment = AlignIO.read(MSA_FILE, "fasta")

num_sequences = len(alignment)
alignment_length = alignment.get_alignment_length()

print("MUSCLE MSA summary")
print(f"File: {MSA_FILE}")
print(f"Number of sequences: {num_sequences}")
print(f"Alignment length: {alignment_length}")
print("\nSequence IDs:")
for rec in alignment:
    print("-", rec.id)

# Save a summary CSV (first version, one tool only)
csv_path = os.path.join(TABLES_DIR, "msa_summary.csv")
with open(csv_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["tool", "num_sequences", "alignment_length"])
    writer.writerow(["MUSCLE", num_sequences, alignment_length])

print(f"\nSaved summary table: {csv_path}")
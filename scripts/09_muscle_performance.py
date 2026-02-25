import subprocess
import os
import csv
import time
import psutil
from Bio import SeqIO, AlignIO

INPUT_FASTA = "data/cleaned/hemoglobin_beta_cleaned.fasta"
OUTPUT_FASTA = "results/msa/muscle_alignment_perf.fasta"
MUSCLE_EXE = r"tools\muscle\muscle3.8.31_i86win32.exe"
CSV_PATH = "results/tables/msa_performance.csv"

os.makedirs("results/msa", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

# Basic dataset stats
records = list(SeqIO.parse(INPUT_FASTA, "fasta"))
num_sequences = len(records)
lengths = [len(r.seq) for r in records]
avg_length = sum(lengths) / len(lengths) if lengths else 0

cmd = [MUSCLE_EXE, "-in", INPUT_FASTA, "-out", OUTPUT_FASTA]

print("Running MUSCLE for performance measurement...")
print("Command:", " ".join(cmd))

# Measure this Python process memory before/after (approximate)
proc = psutil.Process(os.getpid())
mem_before_mb = proc.memory_info().rss / (1024 ** 2)

t0 = time.perf_counter()
result = subprocess.run(cmd, capture_output=True, text=True)
t1 = time.perf_counter()

mem_after_mb = proc.memory_info().rss / (1024 ** 2)
runtime_sec = round(t1 - t0, 4)
mem_delta_mb = round(mem_after_mb - mem_before_mb, 4)

if result.returncode != 0:
    print(result.stdout)
    print(result.stderr)
    raise RuntimeError("MUSCLE failed during performance run.")

# Read alignment to get aligned length
alignment = AlignIO.read(OUTPUT_FASTA, "fasta")
alignment_length = alignment.get_alignment_length()

row = [
    "MUSCLE",
    num_sequences,
    round(avg_length, 2),
    runtime_sec,
    mem_delta_mb,
    alignment_length
]

# Write CSV (create if not exists, otherwise overwrite for now)
with open(CSV_PATH, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["tool", "num_sequences", "avg_seq_length", "runtime_sec", "mem_delta_mb", "alignment_length"])
    writer.writerow(row)

print("\nPerformance measurement complete.")
print("Result row:", row)
print(f"Saved alignment: {OUTPUT_FASTA}")
print(f"Saved table: {CSV_PATH}")
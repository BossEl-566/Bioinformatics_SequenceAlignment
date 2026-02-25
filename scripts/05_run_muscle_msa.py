import subprocess
import os

INPUT_FASTA = "data/cleaned/hemoglobin_beta_cleaned.fasta"
OUTPUT_FASTA = "results/msa/muscle_alignment.fasta"

# Your MUSCLE executable path (works from your project root)
MUSCLE_EXE = r"tools\muscle\muscle3.8.31_i86win32.exe"

os.makedirs("results/msa", exist_ok=True)

# MUSCLE v3 syntax
cmd = [
    MUSCLE_EXE,
    "-in", INPUT_FASTA,
    "-out", OUTPUT_FASTA
]

print("Running MUSCLE MSA...")
print("Command:", " ".join(cmd))

result = subprocess.run(cmd, capture_output=True, text=True)

print("\nReturn code:", result.returncode)

if result.stdout.strip():
    print("\nSTDOUT:")
    print(result.stdout)

if result.stderr.strip():
    print("\nSTDERR:")
    print(result.stderr)

if result.returncode != 0:
    raise RuntimeError("MUSCLE failed. Check the messages above.")

print(f"\nMSA complete. Output saved to: {OUTPUT_FASTA}")
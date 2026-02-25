from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

INPUT_FASTA = "data/raw/hemoglobin_beta_all.fasta"
OUTPUT_FASTA = "data/cleaned/hemoglobin_beta_cleaned.fasta"

os.makedirs("data/cleaned", exist_ok=True)

def clean_sequence(seq):
    # Remove spaces, line breaks, and existing gaps if any
    return str(seq).upper().replace("-", "").replace(" ", "").replace("\n", "")

cleaned_records = []
lengths = []

for rec in SeqIO.parse(INPUT_FASTA, "fasta"):
    cleaned_seq = clean_sequence(rec.seq)

    # Skip empty sequences
    if not cleaned_seq:
        continue

    # Save cleaned record
    new_record = SeqRecord(
        Seq(cleaned_seq),
        id=rec.id,
        description=rec.description
    )
    cleaned_records.append(new_record)
    lengths.append(len(cleaned_seq))

# Write cleaned FASTA
SeqIO.write(cleaned_records, OUTPUT_FASTA, "fasta")

# Print summary
print(f"Input file: {INPUT_FASTA}")
print(f"Output file: {OUTPUT_FASTA}")
print(f"Sequences loaded and cleaned: {len(cleaned_records)}")

if lengths:
    print(f"Min length: {min(lengths)}")
    print(f"Max length: {max(lengths)}")
    print(f"Average length: {sum(lengths)/len(lengths):.2f}")

print("\nSequence IDs:")
for r in cleaned_records:
    print("-", r.id)
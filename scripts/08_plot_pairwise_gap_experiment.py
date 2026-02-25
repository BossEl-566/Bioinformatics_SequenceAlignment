import csv
import os
import matplotlib.pyplot as plt # type: ignore

INPUT_CSV = "results/tables/pairwise_gap_experiments_human_vs_zebrafish.csv"
OUTPUT_PNG = "results/figures/pairwise_gap_experiment_scores.png"

os.makedirs("results/figures", exist_ok=True)

global_labels = []
global_scores = []
local_labels = []
local_scores = []

with open(INPUT_CSV, "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        label = f"open={row['gap_open']}, ext={row['gap_extend']}"
        score = float(row["score"])
        mode = row["mode"].strip().lower()

        if mode == "global":
            global_labels.append(label)
            global_scores.append(score)
        elif mode == "local":
            local_labels.append(label)
            local_scores.append(score)

# Plot global scores
plt.figure(figsize=(8, 5))
plt.bar(global_labels, global_scores)
plt.title("Global Alignment Score vs Gap Penalties (Human vs Zebrafish HBB)")
plt.xlabel("Gap penalty settings")
plt.ylabel("Alignment score")
plt.xticks(rotation=20, ha="right")
plt.ylim(min(global_scores) - 2, max(global_scores) + 2)
plt.tight_layout()
plt.savefig("results/figures/global_gap_scores.png", dpi=300)
plt.close()

# Plot local scores
plt.figure(figsize=(8, 5))
plt.bar(local_labels, local_scores)
plt.title("Local Alignment Score vs Gap Penalties (Human vs Zebrafish HBB)")
plt.xlabel("Gap penalty settings")
plt.ylabel("Alignment score")
plt.xticks(rotation=20, ha="right")
plt.ylim(min(local_scores) - 2, max(local_scores) + 2)
plt.tight_layout()
plt.savefig("results/figures/local_gap_scores.png", dpi=300)
plt.close()

print("Plots created successfully.")
print("Saved: results/figures/global_gap_scores.png")
print("Saved: results/figures/local_gap_scores.png")
### **Slide 1 — Title Slide**

**Title:** Sequence Alignment with Biopython
**Course:** DCIT 411: Bioinformatics
**Project Title:** Sequence Alignment with Biopython (Hemoglobin Beta Case Study)
**Your Name**
**Submission Date:** February 27, 2026



### **Slide 2 — Project Overview**

* What is sequence alignment?
* Why it matters in bioinformatics
* Project focus: Hemoglobin beta (HBB) protein across vertebrates
* Tools used: Biopython, MUSCLE, UniProt, Matplotlib



### **Slide 3 — Objectives**

Use a shortened version of your report objectives:

* Collect and preprocess HBB protein sequences
* Perform pairwise alignment (global/local)
* Test gap penalties
* Run multiple sequence alignment (MSA)
* Generate consensus sequence
* Measure runtime and memory performance



### **Slide 4 — Dataset (Sequences Used)**

* **Protein family:** Hemoglobin beta (HBB)
* **Source:** UniProt (FASTA format)
* **Total sequences:** 11
* **Species included:** Human, Chimp, Mouse, Rat, Cow, Sheep, Goat, Pig, Chicken, Zebrafish, Frog
* Mention cleaned dataset:

  * `hemoglobin_beta_cleaned.fasta`


### **Slide 5 — Methodology Workflow**

Show as a flow (boxes/arrows):

1. Sequence collection (UniProt FASTA)
2. Preprocessing (clean FASTA using Biopython)
3. Pairwise alignment (Biopython PairwiseAligner + BLOSUM62)
4. Gap penalty experiment
5. MSA with MUSCLE
6. Consensus sequence generation
7. Performance profiling + plotting
8. Report and interpretation




### **Slide 6 — Pairwise Alignment Results (Human vs Chimp)**

* Global score: **780.0**
* Local score: **780.0**
* Observation: complete identity in selected HBB sequences
* Show a small screenshot/snippet of the alignment with `||||||||`

**Key message:** confirms pipeline correctness and high conservation in closely related species.



### **Slide 7 — Pairwise Alignment Results (Human vs Zebrafish)**

* Global score: **419.0**
* Local score: **419.0**
* More mismatches and gap(s)
* Conserved regions still present

**Key message:** shows evolutionary divergence but retained functional conservation.



### **Slide 8 — Gap Penalty Experiment (Human vs Zebrafish)**

Use your chart(s) here.

* `(-8, -0.2)` → **422.0**
* `(-10, -0.5)` → **419.0**
* `(-12, -1.0)` → **417.0**

**Observation:** stricter gap penalties reduce alignment score.

(Insert your global/local plot image(s))

---

### **Slide 9 — MSA + Consensus + Performance**

* MUSCLE aligned **11 sequences**
* Alignment length: **148**
* Consensus generated (threshold = 0.7)
* Runtime: **0.2146 s**
* Memory delta: **0.3125 MB**

Show:

* a short MSA snippet (with gaps `-`)
* first part of consensus sequence:

  * `MVXLTXXEKXAVXXXWGKV...`

---

### **Slide 10 — Conclusion and Future Work**

**Conclusion**

* Successfully implemented pairwise and multiple sequence alignment using Biopython and MUSCLE
* Identified conserved and variable HBB regions across vertebrates
* Demonstrated gap penalty effects and generated consensus sequence

**Future Work**

* Compare MUSCLE with MAFFT / ClustalW
* Use larger datasets
* Add structural alignment / HMM-based methods





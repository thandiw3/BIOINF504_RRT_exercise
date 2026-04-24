### The Reproducibility Audit: Find the Flaws, Fix the Pipeline

A ~60-minute practical exercise (created with the use of claude.ai) for graduate bioinformatics students. Students audit a deliberately flawed RNA-seq analysis for violations of statistical rigor, software engineering best practices, and FAIR data principles, then propose concrete fixes.

Companion repository (has the data files and code to audit): https://github.com/mitreacristina/RNAseq_analysis
---

### Files

| File | Description |
|---|---|
| `RTT_exercise_description.md` | **Student-facing assignment** — start here |
| `counts_final_FINAL_v3.csv` | Simulated RNA-seq count matrix (200 genes × 20 samples) |
| `meta.csv` | Sample metadata (condition, cohort, batch) |
| `fastq_samples.tar.gz` | 20 simulated FASTQ files, one per sample (500 reads × 75 bp, gzipped) |

---

### Simulated Data

Data were generated with `numpy.random.seed(42)` using a negative binomial model (dispersion 0.1–0.3).

**Ground truth DE genes:** `GENE_42`, `GENE_107`, `GENE_156` — 4× upregulated in tumor vs. normal.  
**Batch effect:** All 2022 samples (Patient_11–Patient_20) carry a 1.2× multiplicative expression increase across all genes.  
**Cohorts:** Cohort A = Patient_01–12 (n = 6 tumor, n = 6 normal, batch 2018/2022 mixed); Cohort B = Patient_13–20 (n = 4 tumor, n = 4 normal, batch 2022).

Running the broken script recovers ~5–15 genes at p < 0.05, including false positives driven by the batch effect. The corrected pipeline (DESeq2 / pydeseq2 with batch covariate and BH-FDR) recovers `GENE_42` and `GENE_107` reliably at n = 10 per group.

---

### Topics Covered

- Statistical power and multiple testing correction
- Appropriate models for count data (negative binomial vs. t-test)
- p-hacking and HARKing
- Batch effects and confounding
- Biological vs. technical replicates
- Code reproducibility: version pinning, containerization, workflow managers
- FAIR data principles
- Sampling bias and team diversity

---

### Prerequisites

Students should be familiar with R or Python and have completed course modules on statistics, software engineering, and FAIR principles before attempting this exercise.

**R packages used in solutions:** `pwr`, `DESeq2`, `RnaSeqSampleSize`, `renv`  
**Python packages used in solutions:** `pandas`, `numpy`, `scipy`, `pydeseq2`, `statsmodels`, `rpy2`
 

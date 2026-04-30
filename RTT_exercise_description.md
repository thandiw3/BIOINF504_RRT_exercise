## The Reproducibility Audit: Find the Flaws, Fix the Pipeline

**Estimated time:** ~60 minutes  
**Format:** Individual  
**Deliverable:** Completed tasks submitted as a pull request or a single PDF or `.md` file

#### Created with the use of claude.ai

---

### Learning Objectives

By the end of this exercise you should be able to:

- Identify p-hacking and underpowered study design
- Spot missing randomization and confounded replicates
- Evaluate code for versioning and reproducibility gaps
- Apply FAIR principles to data and metadata
- Recognize sampling bias and demographic confounders
- Propose concrete, actionable remediation steps
- Prioritize which flaws matter most in a real-world setting

---

### Scenario

A collaborator shares a GitHub repository for a published RNA-seq differential expression study comparing gene expression in tumor vs. normal tissue across two cohorts. They are excited — several genes hit *p* < 0.05 and they want to submit a follow-up grant. Before you endorse the work, you agree to do a quick reproducibility audit.

**Your job: find everything that could go wrong, and recommend fixes.**

---

### Background: Key Terms

Before you begin, make sure you are comfortable with the following terms used throughout this assignment.

- **HARKing** (Hypothesizing After Results are Known): reporting a finding as if it were an a priori hypothesis when it was actually identified after inspecting the data. This inflates the apparent significance of results.
- **Multiple testing correction**: when testing many hypotheses simultaneously (e.g., 18,000 genes), the probability of at least one false positive at α = 0.05 approaches certainty. Procedures such as Benjamini–Hochberg FDR control the expected proportion of false discoveries.
- **Negative binomial model**: RNA-seq counts are discrete, non-negative integers. They are overdispersed relative to a Poisson distribution (variance > mean), which makes the negative binomial distribution a better fit. Tools like DESeq2 (R) and pydeseq2 (Python) use this model.
- **CV (coefficient of variation)**: CV = σ/μ, the ratio of standard deviation to mean. A CV of 0.4 means the standard deviation is 40% of the mean — a typical level of biological variability in RNA-seq data.
- **FAIR principles**: a framework for data sharing — data should be Findable, Accessible, Interoperable, and Reusable.
- **Batch effect**: systematic technical variation introduced by differences in sample processing across time, lab, reagent lot, or sequencing run, which can masquerade as biological signal.

---

### Step 1 — Read the Methods Snapshot *(~7 min)*

#### Provided methods excerpt

> "We analyzed RNA-seq data from Cohort A (n = 6 tumor, n = 6 normal) and Cohort B (n = 4 tumor, n = 4 normal). All patients were male, aged 55–70, recruited at a single institution between 2018–2022. Library preparation and sequencing were performed in two runs: samples collected in 2018 were processed together, and samples collected in 2022 were processed together. Three tumor samples from Patient_01 were re-sequenced at higher depth and included as separate entries. We tested 18,000 genes for differential expression using a t-test and reported all genes with p < 0.05 without multiple testing correction. We selected the top 12 candidate genes based on fold-change ranking after seeing the data. The analysis was run in [Python/R]; the version is not recorded and no environment file was committed. Raw FASTQ files are available upon reasonable request. The pipeline script is included in the repo but requires manual setup with no version guidance. All samples are labeled 'Patient_01' through 'Patient_20' with no additional metadata file. Both authors are from the same lab and share the same clinical specialty (oncology)."

#### Task 1A — Annotation checklist

For each item below, mark whether a problem is present and write one sentence describing it. If a problem is not present or cannot be determined from the text, say so and explain.

| # | Category | Problem present? | Description |
|---|---|---|---|
| 1 | Sample size / statistical power | Yes |The study uses only 6 tumor and 6 normal samples in Cohort A and 4 tumor and 4 normal samples in Cohort B, which is likely underpowered for transcriptome-wide differential expression across 18,000 genes |
| 2 | Multiple testing correction |Yes | The authors tested 18,000 genes and reported genes with raw p<0.05, which makes false positives almost guaranteed without FDR correction. |
| 3 | P-hacking / HARKing |Yes |The top 12 candidate genes were selected after seeing the data, which means the final candidates may reflect post-hoc ranking rather than a pre-specified hypothesis |
| 4 | Appropriate statistical test for count data |Yes |A t-test is not appropriate for raw RNA-seq count data because RNA-seq counts are discrete and overdispersed, so a negative binomial model such as DESeq2 or pydeseq2 is more appropriate.|
| 5 | Randomization and batch effects |Yes |Samples from 2018 and 2022 were processed in separate sequencing runs, so year, batch, and biological condition may be confounded. |
| 6 | Biological vs. technical replicates |Yes |Three tumor samples from Patient_01 were re-sequenced and included as separate entries, which treats technical replicates as independent biological samples and inflates the apparent sample size.|
| 7 | Version pinning and dependency management |Yes |The Python/R version and package versions were not recorded, making it difficult to reproduce the analysis later.|
| 8 | Code portability / containerization |Yes |The pipeline requires manual setup with no environment guidance or container, so it may only run on the original author’s machine.|
| 9 | FAIR — Findability and Accessibility of raw data |Yes |Raw FASTQ files are only “available upon reasonable request,” which makes the data less findable and less reliably accessible than deposition in GEO, SRA, Zenodo, or another repository. |
| 10 | FAIR — Interoperability and metadata completeness |Yes |Samples are only labeled Patient_01 through Patient_20 with no complete metadata file, so key variables such as condition, batch, year, replicate status, and sequencing run are not machine-readable. |
| 11 | FAIR — Reusability, licensing, and provenance |Yes |No environment file, license, full metadata, or provenance information is provided, so another researcher cannot easily reuse or validate the data and pipeline. |
| 12 | Single-institution / demographic sampling bias |Yes |All patients were male, aged 55–70, and recruited from one institution, so the results may not generalize to other sexes, ages, institutions, or patient populations. |
| 13 | Confirmation bias from homogeneous team |Maybe |Both authors are from the same lab and specialty, which increases the chance that assumptions about the analysis go unchallenged. |
| 14 | Domain blind spots from homogeneous team |Yes |The team doesnt have a dedicated statistical genomics or computational expert, which explains the inappropriate t-test, missing FDR correction, and poor pipeline documentation. |

> **Hint:** Not every category has a single flaw — some have multiple layered problems. The goal is systematic thinking, not a perfect list.

---

### Step 2 — Inspect the Code *(~12 min)*

Chose **R** . Review the script for software engineering and reproducibility issues.

---

#### Option B — R

```r
# Use DESeq2 because RNA-seq counts are discrete and overdispersed,
# so a negative binomial model is more appropriate than gene-wise t-tests.
library(DESeq2)

# Load count matrix and metadata.
# Counts should be genes x samples, and metadata should have one row per sample.
counts <- read.csv("counts_final_FINAL_v3.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("meta.csv", row.names = 1)

# Make sure the sample order matches between the count matrix and metadata.
counts <- counts[, rownames(metadata)]

# Convert condition to a factor and set normal as the reference group.
# This makes the tumor-vs-normal contrast interpretable.
metadata$condition <- relevel(factor(metadata$condition), ref = "normal")

# Use a pre-specified low-count filter to remove genes with too little information.
# This avoids testing genes that are mostly zero while preventing post-hoc threshold tuning.
keep <- rowSums(counts >= 10) >= 3
counts_filtered <- counts[keep, ]

# Build a DESeq2 dataset using the experimental design.
# If batch is available in metadata, use design = ~ batch + condition.
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_filtered)),
  colData = metadata,
  design = ~ condition
)

# Fit the negative binomial model and estimate gene-wise dispersion.
dds <- DESeq(dds)

# Test the tumor vs normal contrast.
res <- results(
  dds,
  contrast = c("condition", "tumor", "normal"),
  alpha = 0.05
)

# DESeq2 reports Benjamini-Hochberg adjusted p-values in the padj column.
# This controls the false discovery rate across the tested genes.
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Optional explicit BH correction, included to show the correction step directly.
res_df$padj_BH_explicit <- p.adjust(res_df$pvalue, method = "BH")

# Define significant genes using adjusted p-values, not raw p-values.
sig_genes <- subset(res_df, padj < 0.05)

write.csv(res_df, "deseq2_all_results.csv", row.names = FALSE)
write.csv(sig_genes, "deseq2_significant_results_fdr05.csv", row.names = FALSE)
```
I replaced the gene-wise t-tests with DESeq2 because RNA-seq counts are discrete, non-negative, and overdispersed. DESeq2 fits a negative binomial model, estimates dispersion, and reports Benjamini–Hochberg adjusted p-values through the padj column. I also used a pre-specified low-count filter and defined significance by FDR rather than raw p-values.

> **Also check the repo ([RNA_seq_analysis](https://github.com/mitreacristina/RNAseq_analysis)):** No `renv.lock` (note: `renv` supersedes the older `packrat` for R dependency management), no `Dockerfile` or Apptainer (formerly Singularity) definition, no `README`, no `.gitignore`. The only commit message is "final version". Raw FASTQ files are tracked directly in git. There are no unit tests. The required packages must be installed manually with no version guidance.

---

#### Task 2A — Code audit table

Complete the table for your chosen language.

| Issue Category | What's Wrong Here | How to Fix It |
|---|---|---|
| Environment / dependency management | The script does not record the R version or package versions, and the repo has no renv.lock, so another person may get different results depending on what is installed on their machine.| Use renv::init() and renv::snapshot() to create a project-specific renv.lock file that records the exact R and package versions.|
| Containerization |There is no Dockerfile or Apptainer/Singularity definition, so even with the script, the operating system and system libraries are not reproducible. |Add a Dockerfile or Apptainer definition that installs the required R, Bioconductor, and system dependencies.|
| Hardcoded file paths |The script uses setwd("/Users/mitreacristina/Desktop/rnaseq_project"), which only works on the original author’s computer. | Use an RStudio Project with relative paths, or use here::here() so the script can run from the project directory on any machine. |
| File naming / version control |The file name counts_final_FINAL_v3.csv suggests manual versioning and makes it unclear which file is the true input.|Use clear, stable file names such as counts_raw.csv or counts_gene_by_sample.csv, and track changes through Git commits or DVC rather than file-name edits.|
| Large files in git |Raw FASTQ files are tracked directly in git, which makes the repository too large and mixes source code with raw data storage. |Store FASTQ files in GEO, SRA, Zenodo, OSF, or use DVC/Git LFS, while keeping only lightweight code and metadata in Git. |
| Statistical test choice for count data |The script uses t.test() on RNA-seq counts, but RNA-seq counts are discrete, non-negative, and overdispersed rather than normally distributed. |Use DESeq2, which models RNA-seq counts with a negative binomial distribution and estimates gene-specific dispersion. |
| Multiple testing correction |The script reports genes with raw p < 0.05 after testing thousands of genes, which will produce many false positives.|Use Benjamini–Hochberg FDR correction through DESeq2’s padj column or p.adjust(pvals, method = "BH"), and define significance using adjusted p-values. |
| Post-hoc filtering threshold |The comment says the low-count filter was chosen after looking at the results, which makes the analysis circular and can bias significance.|Pre-specify the filtering rule before analysis, such as keeping genes with at least 10 counts in at least 3 samples, and document the rule in the Methods. |
| Manual candidate selection |The candidate genes were picked by hand after inspecting the results, which increases the risk of p-hacking or HARKing. |Define candidate selection criteria before looking at results, such as padj < 0.05 and abs(log2FoldChange) > 1, or validate candidates in an independent cohort. |
| Random seed / stochastic reproducibility |The DESeq2 model itself is mostly deterministic, but the script does not document this, and any later stochastic steps such as sampling, plotting, or resampling would not be reproducible.|Add set.seed(42) before any stochastic step and state when a seed is not needed because the analysis is deterministic. |
| Code documentation |The script has no author, date, R version, purpose statement, input description, or output description. |Add a script header, comments explaining why each step is performed, and a README that describes how to run the RStudio project from start to finish. |
| Unit testing |There are no tests to make sure the metadata and count matrix match, that required columns exist, or that the output contains expected DESeq2 result columns. |Use testthat to test sample ID matching, metadata columns, count matrix dimensions, and presence of log2FoldChange, pvalue, and padj in the output. |
| Commit message quality |The only commit message is “final version,” which does not explain what changed or why. |Use descriptive commit messages such as “add DESeq2 model,” “apply BH FDR correction,” or “add metadata validation tests.” |

#### Task 2B — Rewrite the stats section

Rewrite the testing section (the loop or `apply` block) using a statistically appropriate method for RNA-seq count data.

RNA-seq counts are discrete, non-negative integers that exhibit **overdispersion** — the variance exceeds the mean, which violates both the normality assumption of the t-test and the equidispersion assumption of a Poisson model. The **negative binomial distribution** explicitly models this overdispersion and is the standard choice for RNA-seq differential expression.

Your rewrite must:

1. Use the correct package and key function calls with their main arguments — pseudocode is acceptable but must show the analysis pipeline steps, not just name the package
   - Python: `pydeseq2` (`DeseqDataSet` → `deseq2()` → `DeseqStats` → `results_df`)
   - R: `DESeq2` (`DESeqDataSetFromMatrix` → `DESeq` → `results`)
2. Apply Benjamini–Hochberg FDR correction and name the function used
3. Set a random seed (Python) or document why one is not needed (R)
4. Include at least one comment per major step explaining *why*, not just *what*

#### Task 2C — Environment reproducibility

Rather than simply writing an environment file, answer the following:

1. Name two specific things that would **silently break** if a collaborator ran the original script two years from now without any pinned versions. Be concrete — name a package and a type of change that could occur.
First, a future version of DESeq2 could change defaults for dispersion estimation, independent filtering, or result formatting, which could change adjusted p-values or output columns without making the script obviously fail. Second, a future version of readr, utils, or tibble/data.frame handling could change how sample IDs, row names, or strings are imported, causing the metadata and count matrix to misalign silently.

2. Write a minimal `environment.yml` (Python/Conda) or `renv` initialization (R) pinning at least 5 relevant packages to specific versions.
# Initialize a project-local R environment
install.packages("renv")
renv::init()

# Install required packages
install.packages(c(
  "BiocManager",
  "here",
  "readr",
  "dplyr",
  "testthat"
))

BiocManager::install(c(
  "DESeq2",
  "RnaSeqSampleSize"
))

# Record exact package versions in renv.lock
renv::snapshot()

3. In one sentence: when is a pinned environment file alone *not* sufficient for full reproducibility, and what additional tool addresses this?
A pinned renv.lock file is not sufficient when system libraries, operating system versions, or external command-line tools affect the analysis; a Docker or Apptainer container addresses this by preserving the full computational environment.

---

### Step 3 — FAIR & Bias Assessment *(~9 min)*

#### Task 3A — FAIR score

Rate each dimension on a 1–5 scale (1 = completely absent, 5 = fully compliant) and justify in one sentence. Refer to specific evidence from the methods excerpt.

| Dimension | Score (1–5) | Justification (cite specific evidence from the methods) |
|---|---|---|
| Findable | 1 | The raw FASTQ files are not deposited in GEO, SRA, Zenodo, or another indexed repository, so the data cannot be found through standard search or accession systems. |
| Accessible | 2| “Available upon reasonable request” is better than no access, but it depends on author response and does not guarantee long-term availability.|
| Interoperable |1 |The samples are labeled only Patient_01 through Patient_20 with no complete metadata file, so the data cannot be easily interpreted or integrated with other datasets. |
| Reusable |1 |The lack of license, provenance, environment file, complete metadata, and documented pipeline makes reuse difficult and potentially ambiguous. |

> Consider: Does "available upon reasonable request" satisfy Accessibility? Does the absence of a `LICENSE` file affect Reusability? Would depositing raw data in NCBI Gene Expression Omnibus (GEO) or Zenodo change the Findability score?

#### Task 3B — Bias identification

List at least **three distinct sources of bias or confounding** present in this study. For each one provide:
1. The **type of bias** (e.g., selection bias, confirmation bias, technical confounding)
2. The **specific claim in the methods excerpt** that is undermined — do not write generic answers like "it threatens validity"
3. One concrete **mitigation strategy** with a named tool or approach

- Samples collected in 2018 were processed together and samples collected in 2022 were processed together. This undermines the claim that observed expression differences are caused by tumor-normal status because sequencing year or library preparation batch could be driving the signal. A concrete mitigation would be to randomize tumor and normal samples across library preparation batches and include batch as a covariate in a DESeq2 or pydeseq2 design formula, such as ~ batch + condition.
- Three tumor samples from Patient_01 were re-sequenced and included as separate entries. This undermines the claim that the study has 20 independent samples because resequenced samples from the same patient are technical replicates, not new biological individuals. A concrete mitigation would be to collapse technical replicates using a tool such as DESeq2’s collapseReplicates() or to model patient identity explicitly.
- All patients were male, aged 55–70, and recruited at a single institution. This undermines any broad claim that the tumor-normal signature generalizes across patient populations. A concrete mitigation would be to recruit a multi-institution cohort with sex, age, ancestry, and clinical covariates recorded in a structured metadata file.
- The methods say that the top 12 candidate genes were selected after seeing the data. This undermines the claim that these genes were hypothesis-driven candidates because the selection was post-hoc. A concrete mitigation would be to pre-register candidate selection rules on OSF or define candidates using pre-specified FDR and effect-size thresholds.

#### Task 3C — Team composition

The methods note that both authors are from the same lab and share the same clinical specialty (oncology). Answer both parts:

1. **Confirmation bias:** How does same-lab authorship increase the risk of confirmation bias, and what is one structural safeguard that can be put in place at the analysis stage?
Same-lab authorship increases the risk of confirmation bias because both authors may share the same assumptions about the disease, the dataset, and what counts as a convincing result. A structural safeguard would be to use a blinded or pre-registered analysis plan, where filtering rules, statistical models, and candidate-selection criteria are finalized before inspecting differential expression results.

2. **Domain blind spots:** What type of expertise is missing from this two-person oncology team that would have caught the core statistical flaw in this script, and how would you recruit that expertise at the study design phase?
The missing expertise is statistical genomics or computational biology, especially someone familiar with RNA-seq count models, multiple testing, and reproducible pipelines. I would recruit a biostatistician or bioinformatics methods collaborator during the study design phase, before sequencing, so they can help define sample size, randomization, metadata, batch structure, and the differential expression model.
---

### Step 4 — Power Calculation & Design Proposal *(~8 min)*

------------------------------
#### Task 4A — Was the study powered?

Cohort A has n = 6 per group. Use the formula below to estimate whether this study was adequately powered to detect a 2-fold change in expression.

> **Important caveat:** The formula below applies to a two-sample t-test on continuous, normally distributed data — the same (incorrect) test used in the script. This gives a lower-bound estimate. For a properly designed RNA-seq study, use purpose-built tools such as `RnaSeqSampleSize` (R/Bioconductor) or `RNASeqPower` (R), which use negative binomial models. The rule of thumb for pilot RNA-seq studies is a minimum of 3 biological replicates per group; for discovery studies, ≥10–20 per group is typical after accounting for multiple testing.

---

##### Part 1 — Worked example: t-test power (microarray context)

The t-test power formula is appropriate for **continuous, normally distributed data** — for example, log₂ microarray intensities. We work through it here to build intuition before moving to the correct model for RNA-seq.

You are designing a microarray study comparing tumor vs. normal tissue. From pilot data:

| Parameter | Symbol | Value |
|---|---|---|
| Significance threshold | α | 0.05 (two-tailed) |
| Desired power | 1 − β | 0.80 |
| Standard deviation of log₂ intensities | σ | 0.8 |
| Minimum difference to detect | δ | 1.0 (a 2-fold change on the log₂ scale) |

**Formula:**

$$n \approx 2 \times \left[\frac{(z_{\alpha/2} + z_\beta) \times \sigma}{\delta}\right]^2$$

where z_α/2 = 1.96 and z_β = 0.84.

**Calculation:**

$$n \approx 2 \times \left[\frac{(1.96 + 0.84) \times 0.8}{1.0}\right]^2 = 2 \times [2.24]^2 \approx \mathbf{11 \text{ per group}}$$

Verify with code — both snippets should return n = 11:

**R**
```r
library(pwr)
result <- pwr.t.test(
  d           = 1.0 / 0.8,   # Cohen's d = delta / sigma
  sig.level   = 0.05,
  power       = 0.80,
  type        = "two.sample",
  alternative = "two.sided"
)
ceiling(result$n)
```

**Python**
```python
import math
from statsmodels.stats.power import TTestIndPower

n = TTestIndPower().solve_power(
    effect_size = 1.0 / 0.8,  # Cohen's d = delta / sigma
    alpha       = 0.05,
    power       = 0.80,
    alternative = "two-sided"
)
print(math.ceil(n))
```

> **Why this formula does not apply to RNA-seq:** RNA-seq counts are discrete, non-negative integers whose variance scales with the mean — they follow a **negative binomial distribution** with Var(X) = μ + μ²·φ, where φ is the gene-specific dispersion (φ = BCV², the squared biological coefficient of variation). The t-test formula assumes constant variance and normality. More importantly, this is a high-throughput experiment testing 18,000 genes simultaneously — the single-gene formula has no concept of multiple testing and will severely underestimate the required sample size.

---

#### Part 2 — Correct power calculation for RNA-seq

Because we are selecting DE genes across the full transcriptome, the power calculation must account for the FDR correction applied across all 18,000 tests. Use `RnaSeqSampleSize` (R/Bioconductor), which models the negative binomial distribution and derives the per-gene significance threshold directly from the genome-wide FDR target (Zhao et al., *BMC Bioinformatics*, 2018).

**Parameters:**

| Parameter | Value | Notes |
|---|---|---|
| Total genes tested | 18,000 | As stated in the methods |
| Expected DE genes | 900 | 5% of 18,000 — a reasonable prior for tumor vs. normal |
| Fold change to detect | 2 | Minimum biologically meaningful change |
| Mean count (control) | 5 | Average read depth for a DE gene of interest |
| Dispersion (φ = BCV²) | 0.16 | BCV = 0.4, typical for human RNA-seq; φ = 0.4² = 0.16 |
| FDR threshold | 0.05 | Benjamini–Hochberg |
| Target power | 0.80 | |

**R**
```r
# install if needed: BiocManager::install("RnaSeqSampleSize")
library(RnaSeqSampleSize)

n <- sample_size(
  power   = 0.80,   # desired power (1 - Type II error rate)
  m       = 18000,  # total number of genes tested
  m1      = 900,    # expected number of truly DE genes
  f       = 0.05,   # FDR threshold (Benjamini-Hochberg)
  rho     = 2,      # fold change to detect (Treatment / Control)
  lambda0 = 5,      # mean read count for a DE gene in the control group
  phi0    = 0.16,   # dispersion = BCV^2 = 0.4^2; typical for human RNA-seq
  k       = 1,      # ratio of group sizes (1 = balanced design)
  w       = 1       # ratio of normalization factors (1 = no systematic bias)
)
print(n)
```

**Python** — call the R function directly via `rpy2`:
```python
import math
from rpy2.robjects.packages import importr
from rpy2 import robjects

sample_size = importr("RnaSeqSampleSize").sample_size

n = sample_size(
    power   = robjects.FloatVector([0.80]),
    m       = robjects.IntVector([18000]),
    m1      = robjects.IntVector([900]),
    f       = robjects.FloatVector([0.05]),
    rho     = robjects.FloatVector([2]),
    lambda0 = robjects.FloatVector([5]),
    phi0    = robjects.FloatVector([0.16]),
    k       = robjects.FloatVector([1]),
    w       = robjects.FloatVector([1])
)
print(math.ceil(float(n[0])))
```

> **Expected result:** Approximately **20–25 samples per group** — consistent with published benchmarks for human RNA-seq studies with BCV = 0.4 and a 2-fold detection threshold. This is three to four times the n = 6 used in Cohort A, and five times the n = 4 in Cohort B.

**Answer these questions:**

1. What minimum n per group does `sample_size()` return with these parameters?
- Using the RNA-seq-specific power calculation with RnaSeqSampleSize::sample_size(), the required sample size is approximately 20–25 samples per group. This is much larger than the sample size used in the scenario. Cohort A has only 6 tumor and 6 normal samples, and Cohort B has only 4 tumor and 4 normal samples. Therefore, the study is underpowered for transcriptome-wide differential expression after accounting for overdispersion and FDR correction.
The t-test calculation gives 11 samples per group, but that is only the simplified worked example. It is not the correct final answer for RNA-seq because it assumes continuous normally distributed data and does not account for testing 18,000 genes.
2. How does the answer change if you assume only 1% of genes are truly DE (m1 = 180)? What does this tell you about how sensitive the power estimate is to your assumptions?
- If only 1% of genes are truly differentially expressed, meaning m1 = 180 instead of m1 = 900, the required sample size increases. This happens because there are fewer true biological signals to detect among the same 18,000 tests, so it becomes harder to find real differentially expressed genes while still controlling the false discovery rate.
This shows that the power estimate is sensitive to the assumptions built into the calculation. If the expected number of truly DE genes is too optimistic, then the study may look more powered than it actually is. In a real RNA-seq study, assumptions about dispersion, mean count, fold change, FDR threshold, and the expected proportion of DE genes all strongly affect the final sample-size estimate
3. The scenario uses n = 6 (Cohort A) and n = 4 (Cohort B). Beyond missing true DE genes, what is one additional scientific consequence of running an underpowered high-throughput study?
- Beyond missing true differentially expressed genes, an underpowered high-throughput study can produce unstable and inflated effect-size estimates. This means that the genes that appear significant may look more strongly changed than they really are, simply because only the most extreme noisy results pass the threshold. As a result, researchers may prioritize the wrong genes for follow-up experiments, build a grant around weak signals, or report findings that do not replicate in a larger cohort.

#### Task 4B — Redesign in 5 bullets

Write exactly **5 bullet points** describing an improved study design, one per course principle listed below. Each bullet must include: (a) one named tool or approach, (b) one sentence explaining which specific flaw it addresses and why, and (c) one limitation — something this fix does *not* solve.
- Statistics: Use RnaSeqSampleSize before sequencing and DESeq2 or pydeseq2 for differential expression, because this addresses the underpowered design and inappropriate t-test by modeling RNA-seq counts with a negative binomial distribution; this does not solve demographic sampling bias if the cohort is still narrow.
- Team diversity: Add collaborators that have a strong background in statistics and computational reproducibility to the study design stage, because their perspectives wouldve likely caught the t-test, FDR, and pipeline issues before publication; this does not eliminate the need for external validation in an independent cohort.
- FAIR: Deposit raw FASTQ files and processed count matrices in GEO or SRA and archive the analysis release on Zenodo, because this makes the data findable, accessible, and citable; this does not by itself ensure that the metadata are complete or biologically balanced.
- Bias/confounding: Randomize tumor and normal samples across sequencing runs and include batch in the design formula, because this prevents year or sequencing run from masquerading as tumor-normal biology; it does not remove all possible confounders, such as institution-level or demographic bias.
- 
Principles to address (one each): **statistics · software engineering · FAIR · bias/confounding · team diversity**

> Python tools to draw from: Snakemake, Zenodo, Docker, GEO, `pydeseq2`, `pytest`, DVC, OSF  
> R tools to draw from: Nextflow, Zenodo, Apptainer, GEO, `DESeq2`, `testthat`, `renv`, `RnaSeqSampleSize`

---

### Step 5 — Synthesis *(~7 min)*

#### Task 5A — Reflection

In 3–5 sentences: *Which single flaw in this study do you think is most common in published bioinformatics literature, and why is it so persistent?* Draw on at least one concept from the course. There is no single correct answer — we are looking for critical reasoning and specificity.
- The most common flaw in published bioinformatics literature is probably incomplete reproducibility rather than total absence of code. 
Many papers provide a script or notebook, but the analysis still depends on undocumented package versions, hardcoded paths, missing metadata, or manual choices that are not captured in the repository. 
I think this happens because bioinformatics projects can develop through exploration, and the final paper presents the cleaned-up story rather than the full computational history. 
From a reproducibility perspective, this is dangerous because a result can appear rigorous while still being difficult or impossible to regenerate exactly.

#### Task 5B — Prioritization *(~7 min)*

You have one week before your collaborator submits the grant. You cannot fix everything. **Choose the single most important flaw to address first.** In 3–5 sentences, defend your choice. Your answer must: (1) name the flaw, (2) explain what harm it causes if left unfixed, and (3) acknowledge the strongest counter-argument for fixing a different flaw instead.
- The single most important flaw to address first is the inappropriate statistical analysis, especially the use of raw t-tests without multiple testing correction.
If left unfixed, the claims made may be based on false positives, which means the follow-up grant could be built around genes that do not actually replicate. 
The strongest counter-argument is that batch effects or pseudoreplication may be even more damaging because they can create systematic false signal before the statistical test is even applied. 
However, fixing the differential expression model, FDR correction, and replicate handling together is the fastest way to determine whether there is any credible signal worth pursuing.

> **Optional extension:** Find one real published paper and identify one reproducibility or FAIR issue in its methods or data availability statement. Write two sentences. 

---

### Deliverables *(~10 min for exercise prep and deliverables check)*

Submit all of the following as a pull request to this repository or single PDF or `.md` file via the course portal.

1. Completed annotation checklist (Task 1A) with one sentence per item
2. Completed code audit table (Task 2A), rewritten stats section (Task 2B), and environment reproducibility answers (Task 2C)
3. FAIR scores with evidence-based justifications (Task 3A), bias analysis (Task 3B), and team composition answers (Task 3C parts 1 and 2)
4. Power calculation with workthrough shown and genome-wide interpretation (Task 4A), plus 5-bullet redesign with tool, rationale, and limitation per bullet (Task 4B)
5. Reflection (Task 5A) and prioritization argument (Task 5B)

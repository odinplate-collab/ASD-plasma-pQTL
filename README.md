# Genome-wide plasma pQTL analysis of ASD — analysis code

Analysis code accompanying *"Genome-wide plasma pQTL analysis of ASD reveals complement pathways"*.
The repository is organized by analysis pipeline.

> **Note on data.** Large input files (WGS/VCF, TMT/DIA proteomics matrices, pQTL summary
> statistics, single-cell atlases) are **not** included here — see the **Data availability**
> statement in the paper.
>
> **Scope.** This repository contains the **analysis** code. Scripts that only render figures
> (and perform no statistical analysis) are not included; plotting code is retained only where
> it is interleaved with the analysis it depends on.

---

## Repository layout

| Folder / script | Analysis |
|-----------------|----------|
| **`01_proteomics_TMT/`** | |
| `differential_abundance.R` | Differential protein abundance (ASD vs TD, limma with age/sex), family-shared expression exclusion, and neutrality-weighted priority score. Self-contained; starts from the normalized expression matrix. |
| `covariate_regression.R` | Sex and age effect regression on protein abundance (limma with family duplicate-correlation; covariates: diagnosis, batch, channel). |
| **`02_pQTL_mapping/`** | |
| `pqtl_mapping.R` | Genome-wide cis/trans pQTL mapping with Matrix eQTL (cis ±1 Mb, BH-FDR < 0.05). |
| **`03_mendelian_randomization/`** | |
| `mr_sensitivity.R` | Two-sample MR with instrument-strength and pleiotropy sensitivity: per-SNP F-statistics, Steiger directionality, MR-Egger intercept, MR-PRESSO, Cochran's Q. |
| **`04_colocalization/`** | |
| `coloc_discovery_pqtl.R` | Bayesian colocalization (`coloc.abf`) of the discovery cis-pQTLs against ASD and neurodevelopmental-trait GWAS. |
| `coloc_ukbppp_pqtl.R` | Colocalization using external UKB-PPP pQTL summary statistics against the same GWAS. |
| **`05_single_cell/`** | |
| `single_cell_analysis.R` | Cortical snRNA-seq re-analysis (Velmeshev 2019; Wamsley 2025): cell-type case–control testing with donor-clustered robust standard errors (CRSE), ligand–receptor coupling, and downstream-program scoring. |
| **`06_validation_DIA/`** | |
| `dia_validation.R` | Independent DIA cohort (n = 459): differential-abundance replication, cross-platform fold-change concordance, and genotype–abundance pQTL replication. |
| **`07_risk_score/`** | |
| `risk_score_auc.py` | Diagnostic panel and multi-omics (peptide / pQTL / WGS) ROC-AUC classification models. |
| **`08_replication_UKBPPP/`** | |
| `effect_size_replication.R` | Effect-size concordance of discovery cis-pQTLs against UKB-PPP. |
| `regional_replication.R` | Regional comparison of cis-pQTL association signals (discovery vs UKB-PPP). |

---

## Software / dependencies

**R (≥ 4.4)** — `MatrixEQTL`, `data.table`, `TwoSampleMR`, `MRPRESSO`, `coloc`,
`Seurat` (v5), `sandwich`, `lmtest`, `limma`, `pROC`, `ggplot2`, `Cairo`, `tidyverse`.
**Python** — `07_risk_score/risk_score_auc.py`: `scikit-learn`, `pandas`, `numpy`.
**PLINK v1.9** — LD clumping of MR instruments with a 1000G reference.

External datasets used (see Data availability): PGC ASD GWAS (2019); UK Biobank Pharma
Proteomics Project (Sun et al. 2023); Niu et al. plasma pQTL (2025); Velmeshev et al. (2019)
and Wamsley et al. (2025) cortical snRNA-seq; and GWAS for colocalization (PGC3 SCZ, ADHD,
IQ, educational attainment).

---

## Running the code / paths

All internal paths are **relative to the project (data) root**. Either run scripts with the
working directory set to that root, or export an environment variable:

```bash
export ASD_ROOT=/path/to/ASD          # data root containing TMT/, pQTL/, GWAS/, ...
```

Arrange the deposited data under `ASD_ROOT` as referenced in each script. The analysis
boundary is the **normalized/batch-corrected protein expression matrix** and the
**genotype / pQTL inputs**; upstream raw-MS processing is described in Methods and is not
part of this runnable code.

## Contact
Questions about the code: Jae Won Oh (odinplate@naver.com).

Correspondence: Min-Sik Kim, Hee Jeong Yoo, Joon-Yong An (see paper).

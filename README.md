# Plasma pQTL analysis of autism spectrum disorder: analysis code

Analysis code for a study of genetic regulation of the plasma proteome in autism spectrum
disorder (ASD), integrating whole-genome sequencing with TMT plasma proteomics.

**Scope.** This repository contains analysis code only: statistical computations and the
writing of result tables. It contains no plotting code and no data.

---

## Data

Input data are not included.

- **Plasma proteomics:** ProteomeXchange/PRIDE PXD082931 (TMT discovery cohort) and PXD082910
  (DIA validation cohort).
- **pQTL summary statistics:** Zenodo, https://doi.org/10.5281/zenodo.22068294.
- **Individual-level genotypes:** not publicly available (identifiable human data).
- **Public resources:**
  - ASD GWAS (Grove et al. 2019) and the other outcome GWAS listed in
    `03_genotype_diagnosis_MR/genotype_diagnosis_mr.py`
  - Niu et al. (2025) plasma pQTL (GWAS Catalog GCST90453052–GCST90453597)
  - UK Biobank Pharma Proteomics Project
  - CSF pQTL of Western et al. (2024; GWAS Catalog)
  - eQTL Catalogue and GTEx v8
  - ASD cortex statistics of Gandal et al. (2018, 2022) and Zhang et al. (2023)
  - cortex proteome estimates of Abraham et al. (2019) and Fatemi et al. (2024)
  - snRNA-seq count matrices of Velmeshev et al. (2019; UCSC Cell Browser, dataset
    `autism`) and Wamsley et al. (2024; UCSC Cell Browser, dataset `asd-psychencode`)
  - CellChatDB and OmniPath

## Paths

Every script resolves paths from one environment variable, `ASD_ROOT` (default: current
directory):

```bash
export ASD_ROOT=/path/to/data_root
```

Expected layout under `ASD_ROOT`:

| Folder | Contents |
|---|---|
| `TMT/` | TMT peptide/protein tables, sample information |
| `Validation/` | DIA matrices and metadata; TMT differential abundance table |
| `pQTL/` | genotype dosages (GRCh37), pQTL inputs and outputs (`pQTL/pQTL_re/`) |
| `reference/` | FragPipe search FASTA (UniProt human, one protein per gene) |
| `GWAS/`, `Co_localization/` | GWAS summary statistics; colocalization outputs |
| `derived/` | intermediate tables: sample list, peptide matrix, protein maps, protein-level cis statistics, outcome GWAS and LD reference cut to protein windows, Niu et al. windows |
| `panel12/` | inputs and outputs for the 12-protein panel |
| `brain/` | cortex and single-nucleus inputs and outputs |
| `single_cell/` | Velmeshev et al. raw count matrix |
| `calibration/`, `output/` | outputs |

Each script header lists its exact inputs and outputs.

---

## Scripts

| Script | Analysis |
|---|---|
| **01_proteomics_TMT/** | |
| `differential_abundance.R` | Family-shared expression exclusion (ASD-only up/down proteins), ASD vs others model with age and sex, neutrality-weighted priority score |
| `covariate_regression.R` | Sex and age effects on protein abundance |
| **02_pQTL_mapping/** | |
| `pqtl_mapping.R` | Peptide matrix preparation (VSN, MinProb imputation, ComBat on TMT batch, unique peptides), peptide principal components, cis/trans Matrix eQTL in ASD (n = 48) and TD (n = 42), ASD-only associations and differential-abundance flags, overlap counts |
| `pqtl_calibration_permutation.R`, `pqtl_calibration_table.R` | Calibration: observed versus permuted association counts by P-value threshold (observed run and 100 permutations of genotype labels in ASD; enrichment and empirical FDR) |
| `protein_level_cis.R` | Protein-level cis association, ASD only / TD only / combined (n = 90, diagnosis-adjusted) |
| **03_genotype_diagnosis_MR/** | |
| `instrument_selection.py` | Lead cis variant, F statistic, MAF and retention for the ASD-elevated cis-pQTL proteins |
| `genotype_diagnosis_mr.py` | Genotype + diagnosis and interaction models, within-group slopes, variance explained; Wald-ratio Mendelian randomization with this study's instruments; concordance with Niu et al. |
| `niu_mr_colocalization.py` | Proteome-wide Mendelian randomization (Niu et al. and UKB-PPP instruments) and colocalization |
| `summary_statistics.py` | DerSimonian–Laird pooling, slope agreement, MBL2 within-genotype tests, MR and concordance summaries |
| **04_colocalization/** | |
| `coloc_discovery_pqtl.R` | coloc.abf of this study's cis-pQTLs with the ASD, schizophrenia, ADHD, intelligence and educational attainment GWAS |
| `coloc_niu_pqtl.R` | coloc.abf of Niu et al. pQTL with the ASD GWAS at the same loci |
| **05_prioritization/** | |
| `protein_fisher_test.R` | Genotype × diagnosis Fisher's exact test for ASD pQTL variants combined with differential abundance; uses the cis and ASD-only trans associations of the April 2022 mapping (input files dated 20220413) |
| `prioritization_sets.R` | Intersection of the three prioritization criteria (12 proteins) and set sizes |
| **06_DIA_validation/** | |
| `dia_validation.R` | DIA preprocessing, pQTL replication, cross-platform concordance, pathway replication, five-group statistics |
| **07_diagnostic_models/** | |
| `extract_12_regions.py` | Lead-variant and regional genotypes for the 12 panel proteins |
| `panel_auc.py` | Single-marker balanced AUC; 12-peptide panel and peptide / pQTL / WGS models with repeated cross-validation |
| `panel_null.py` | Random 12-protein panels and label permutation through the same pipeline |
| **08_cortex/a_tissue_and_colocalization/** | |
| `gtex_expression.py` | GTEx v8 median TPM of the 12 panel genes |
| `fetch_csf_pqtl.R`, `fetch_brain_eqtl.R` | CSF pQTL and brain eQTL summary statistics for the 12 panel genes |
| `coloc_plasma_csf.R`, `coloc_plasma_brain.R`, `coloc_best_per_gene.py` | Plasma–CSF and plasma–brain colocalization |
| **08_cortex/b_cortex_transcriptome/** | |
| `module_gene_sets.py` | Pre-specified gene sets, including seven complement cascade layers |
| `gandal2022_primary.py`, `replication_modules.py` | Complement layers and gene sets in ASD cortex (competitive tests); other disorders and datasets |
| `cortex_result_tables.py` | Result tables combining the cortex, plasma and single-nucleus results (complement cascade, 12 panel genes across disorders, C3 and IGFBP5 across layers) |
| **08_cortex/c_single_nucleus/** | |
| `wamsley_pseudobulk.py`, `wamsley_lr_fetch.py`, `velmeshev_pseudobulk.py`, `velmeshev_export.py` | Raw-count pseudobulk per donor × cell type |
| `donor_qc.py` | Donor sex from expression, Dup15q exclusion, donor overlap |
| `lr_pairs_cellchat.R`, `lr_pairs_omnipath.py` | Curated receptors of the panel proteins |
| `sc_de_both.R` | limma-voom per cell type, both datasets |
| `panel12_meta.py`, `microglia_complement_meta.R` | Fixed-effect meta-analysis across the two datasets |

### Run order

1. `01_proteomics_TMT`
2. `02_pQTL_mapping`: `pqtl_mapping.R` → `protein_level_cis.R`. Calibration (independent of
   the rest): `pqtl_calibration_permutation.R` with `PERM_SEED=0` (observed) and
   `PERM_SEED=1` … `100` (permutations), then `pqtl_calibration_table.R`.
3. `03_genotype_diagnosis_MR`: `instrument_selection.py` → `genotype_diagnosis_mr.py`;
   `niu_mr_colocalization.py` (independent) → `summary_statistics.py`
4. `04_colocalization`
5. `05_prioritization`: `protein_fisher_test.R` → `prioritization_sets.R`
6. `06_DIA_validation`
7. `07_diagnostic_models`: `extract_12_regions.py` → `panel_auc.py`, `panel_null.py`
8. `08_cortex`:
   - `module_gene_sets.py` and the two `lr_pairs_*` scripts first
   - then the pseudobulk scripts → `donor_qc.py` → `sc_de_both.R` → the two meta-analysis
     scripts
   - the cortex-transcriptome scripts, then `cortex_result_tables.py` last
   - the colocalization scripts are independent

## Software

- **R 4.4.3:**
  - MatrixEQTL 2.3, data.table 1.17.0, dplyr 1.1.4, tidyverse 2.0.0, readxl 1.4.5
  - limma 3.62.2, edgeR 4.4.2, sva 3.54.0, vsn 3.74.0, DEP 1.28.0, preprocessCore 1.68.0
  - coloc 5.2.3, ensembldb 2.30.0, EnsDb.Hsapiens.v75, Biostrings 2.74.1, Rsamtools 2.22.0
  - CellChat 2.2.0
- **Python 3.9:** numpy 1.26, pandas 2.1, scipy 1.13, statsmodels 0.14, scikit-learn 1.6,
  openpyxl 3.1, bed-reader

## Reproducibility notes

- **Random seeds:** fixed where randomness is used (diagnostic models: 20260922; calibration
  permutation: `PERM_SEED`).
- **Unseeded imputation:** the MinProb imputation of missing peptide intensities in
  `pqtl_mapping.R` was run without a seed. Re-running it reproduces the original matrices
  only up to that imputation (median absolute difference 3.5e-4).
- **Calibration run time:** each calibration run is a full genome-wide cis and trans scan
  (about 30 min); the 101 runs are independent and can be run in parallel.
- **Tied lead variants:** where several variants in perfect LD share identical statistics, the
  variant used is fixed in `instrument_selection.py`.
- **Remote resources:** GTEx, OmniPath, the eQTL Catalogue, the GWAS Catalog and the UCSC Cell
  Browser are queried remotely; results may change if these resources are updated.

## License and contact

MIT License (see `LICENSE`). Questions about the code: Jae Won Oh (odinplate@naver.com).

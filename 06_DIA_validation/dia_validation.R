# DIA validation cohort: preprocessing, pQTL replication, cross-platform concordance,
# pathway-level replication and five-group summaries for the 12-protein panel.
# Inputs (relative to ASD_ROOT):
#   Validation/ASD_459_report.pg_matrix.tsv        DIA-NN protein-group matrix (459 samples)
#   Validation/ASD_459_report.pr_matrix.tsv        DIA-NN precursor matrix
#   Validation/20260311_ASD_DIA_metadata.xlsx      DIA sample groups (ASD/Father/Mother/Sibling/Control)
#   Validation/peptides.tsv                        TMT peptide matrix (discovery cohort)
#   Validation/tmt_differential_abundance.xlsx     TMT differential abundance, ASD vs TD (sheet 1)
#   Validation/Total_ASD_sample_info.xlsx          TMT / DIA / WGS identifier mapping
#   Validation/ASD_pQTL_imputated_case_rm29903_n_alt.tsv   genotype (alt allele count), cases
#   Validation/ASD_pQTL_imputated_ctrl_rm29903_n_alt.tsv   genotype (alt allele count), controls
#   derived/asd_cis_pqtl_peptide.xlsx              ASD cis-pQTL, peptide level
# Outputs (relative to ASD_ROOT):
#   panel12/dia_validation/stats_pqtl_replication.csv     genotype statistics (TMT and DIA)
#   panel12/dia_validation/stats_dep_concordance.csv      TMT log2FC vs DIA delta z (386 proteins)
#   panel12/dia_validation/stats_pathway_replication.csv  fraction of DEPs per pathway
#   panel12/data/dia_group_summary.tsv            n, mean, SE per protein and DIA group
#   panel12/data/dia_ASD_vs_Control.tsv           ASD vs Control (TD) Welch t-test
# Usage: Rscript dia_validation.R

library(tidyverse)
library(readxl)
library(preprocessCore)

ROOT     <- Sys.getenv("ASD_ROOT", ".")
base_dir <- file.path(ROOT, "Validation")
derived_dir <- file.path(ROOT, "derived")
out_dir  <- file.path(ROOT, "panel12/dia_validation")
dat_dir  <- file.path(ROOT, "panel12/data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dat_dir, showWarnings = FALSE, recursive = TRUE)

dia_pg_file      <- file.path(base_dir, "ASD_459_report.pg_matrix.tsv")
dia_pr_file      <- file.path(base_dir, "ASD_459_report.pr_matrix.tsv")
dia_meta_file    <- file.path(base_dir, "20260311_ASD_DIA_metadata.xlsx")
tmt_pep_file     <- file.path(base_dir, "peptides.tsv")
tmt_dep_file     <- file.path(base_dir, "tmt_differential_abundance.xlsx")
sample_info_file <- file.path(base_dir, "Total_ASD_sample_info.xlsx")
geno_case_file   <- file.path(base_dir, "ASD_pQTL_imputated_case_rm29903_n_alt.tsv")
geno_ctrl_file   <- file.path(base_dir, "ASD_pQTL_imputated_ctrl_rm29903_n_alt.tsv")
pqtl_cis_file    <- file.path(derived_dir, "asd_cis_pqtl_peptide.xlsx")

# Diagnosis labels in the 'subject' column of Total_ASD_sample_info.xlsx (Korean data values,
# written as Unicode escapes): "\uc790\ud3d0" = autism (ASD), "\uc815\uc0c1" = typically developing (TD).
SUBJECT_LABEL_ASD <- "\uc790\ud3d0"
SUBJECT_LABEL_TD  <- "\uc815\uc0c1"

risk_genes   <- c("C3","AHSG","MBL2","C1RL","POSTN","BCHE","BTD","ANPEP",
                  "CLEC3B","IGFBP5","PTGDS","QSOX1")
risk_uniprot <- c(C3="P01024", AHSG="P02765", MBL2="P11226", C1RL="Q9NZP8",
                  POSTN="Q15063", BCHE="P06276", BTD="P43251", ANPEP="P15144",
                  CLEC3B="P05452", IGFBP5="P24593", PTGDS="P41222",
                  QSOX1="O00391")
uniprot_to_gene <- setNames(names(risk_uniprot), risk_uniprot)

# Per-row z-score; rows with zero or undefined SD become NA
zscore_rows <- function(m) {
  t(apply(m, 1, function(x) {
    mu <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(NA, length(x)))
    (x - mu) / s
  }))
}


# ==============================================================================
# 1. DIA protein matrix: shared-protein filter, log2, quantile normalization, z-score
# ==============================================================================

cat("Loading data...\n")

dia_raw <- read_tsv(dia_pg_file, show_col_types = FALSE)
dia_sample_cols <- grep("ASD.*\\.raw$", colnames(dia_raw), value = TRUE)
dia_sample_ids  <- str_extract(dia_sample_cols, "ASD_\\d+(?:_\\d+)?(?=_2ug)")

# Drop columns without a parsable sample ID or with a duplicated ID
ok <- !is.na(dia_sample_ids) & !duplicated(dia_sample_ids)
if (any(!ok)) {
  cat(sprintf("WARNING: dropping %d problematic columns\n", sum(!ok)))
  dia_sample_cols <- dia_sample_cols[ok]
  dia_sample_ids  <- dia_sample_ids[ok]
}

dia_clean <- dia_raw %>% select(Protein.Ids, all_of(dia_sample_cols))
colnames(dia_clean) <- c("Protein.Ids", dia_sample_ids)

# DIA metadata; "Control" is the TD group of the DIA cohort
dia_meta <- read_excel(dia_meta_file) %>%
  mutate(Group = factor(Group, levels = c("ASD","Father","Mother","Sibling","Control")))

ctrl_samples <- dia_meta %>% filter(Group == "Control") %>% pull(Sample_id)
ctrl_idx <- which(dia_sample_ids %in% ctrl_samples)
rest_idx <- which(!dia_sample_ids %in% ctrl_samples)

dia_mat <- dia_clean %>% select(all_of(dia_sample_ids)) %>% as.matrix()
dia_mat[dia_mat == 0] <- NA

# Keep proteins detected in >= 70% of BOTH Control and the remaining samples
detect_ctrl <- apply(dia_mat[, ctrl_idx], 1, function(x) mean(!is.na(x) & x > 0))
detect_rest <- apply(dia_mat[, rest_idx], 1, function(x) mean(!is.na(x) & x > 0))
shared <- detect_ctrl >= 0.7 & detect_rest >= 0.7

cat(sprintf("  Shared protein filtering: %d -> %d proteins (>= 70%% in both groups)\n",
            nrow(dia_mat), sum(shared)))

dia_mat   <- dia_mat[shared, ]
dia_clean <- dia_clean[shared, ]

dia_log2  <- log2(dia_mat)
dia_qnorm <- normalize.quantiles(dia_log2)
colnames(dia_qnorm) <- dia_sample_ids

dia_z <- zscore_rows(dia_qnorm)
colnames(dia_z) <- dia_sample_ids

dia_processed <- tibble(Protein.Ids = dia_clean$Protein.Ids) %>%
  bind_cols(as_tibble(dia_z))

cat(sprintf("  DIA protein: %d proteins x %d samples (shared + quantile norm + z-score)\n",
            nrow(dia_processed), length(dia_sample_ids)))


# ==============================================================================
# 2. DIA precursor matrix -> peptide level (same filter and normalization)
# ==============================================================================

cat("Loading DIA precursor matrix...\n")
dia_pr_raw <- read_tsv(dia_pr_file, show_col_types = FALSE)
dia_pr_sample_cols <- grep("^ASD_", colnames(dia_pr_raw), value = TRUE)

# UniProt-to-gene mapping from the precursor matrix
pr_gene_map <- dia_pr_raw %>%
  select(Protein.Ids, Genes) %>%
  distinct(Protein.Ids, Genes) %>%
  filter(!is.na(Genes), Genes != "")

# Collapse charge states: median per Protein.Ids + Stripped.Sequence
dia_pep <- dia_pr_raw %>%
  select(Protein.Ids, Stripped.Sequence, all_of(dia_pr_sample_cols)) %>%
  group_by(Protein.Ids, Stripped.Sequence) %>%
  summarise(across(all_of(dia_pr_sample_cols), ~median(.x, na.rm = TRUE)), .groups = "drop")

pep_mat <- dia_pep %>% select(all_of(dia_pr_sample_cols)) %>% as.matrix()
pep_mat[pep_mat == 0] <- NA

pep_ctrl_idx <- which(dia_pr_sample_cols %in% ctrl_samples)
pep_rest_idx <- which(!dia_pr_sample_cols %in% ctrl_samples)
pep_detect_ctrl <- apply(pep_mat[, pep_ctrl_idx], 1, function(x) mean(!is.na(x)))
pep_detect_rest <- apply(pep_mat[, pep_rest_idx], 1, function(x) mean(!is.na(x)))
pep_shared <- pep_detect_ctrl >= 0.7 & pep_detect_rest >= 0.7

cat(sprintf("  Peptide shared filtering: %d -> %d peptides (>= 70%% in both groups)\n",
            nrow(pep_mat), sum(pep_shared)))

dia_pep <- dia_pep[pep_shared, ]
pep_mat <- pep_mat[pep_shared, ]

pep_log2  <- log2(pep_mat)
pep_qnorm <- normalize.quantiles(pep_log2)
colnames(pep_qnorm) <- dia_pr_sample_cols

pep_z <- zscore_rows(pep_qnorm)
colnames(pep_z) <- dia_pr_sample_cols

dia_pep_processed <- dia_pep %>%
  select(Protein.Ids, Stripped.Sequence) %>%
  bind_cols(as_tibble(pep_z))


# ==============================================================================
# 3. TMT data, sample ID mapping, peptide-level cis-pQTL leads
# ==============================================================================

# TMT differential abundance, ASD vs TD (protein level)
tmt_dep <- read_excel(tmt_dep_file, sheet = 1)

# TMT peptide matrix, z-scored per peptide across all TMT samples
tmt_pep <- read_tsv(tmt_pep_file, show_col_types = FALSE)
tmt_sample_cols <- setdiff(colnames(tmt_pep), c("Peptide","Protein","Gene"))

tmt_pep_z <- tmt_pep %>%
  rowwise() %>%
  mutate(
    row_mean = mean(c_across(all_of(tmt_sample_cols)), na.rm = TRUE),
    row_sd   = sd(c_across(all_of(tmt_sample_cols)), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(across(all_of(tmt_sample_cols), ~(.x - row_mean) / row_sd)) %>%
  select(-row_mean, -row_sd)

# Sample ID mapping (TMT channel, DIA sample, WGS genotype column)
si <- read_excel(sample_info_file, col_names = FALSE)
id_mapping <- tibble(
  subject  = as.character(si[[6]])[-1],
  family   = as.character(si[[4]])[-1],
  subj_n   = as.character(si[[5]])[-1],
  genomics = as.character(si[[15]])[-1],
  tmt_id   = as.character(si[[19]])[-1]
) %>%
  filter(tmt_id != "NO", !is.na(tmt_id)) %>%
  mutate(family_int = as.integer(as.numeric(family)),
         subj_int   = as.integer(as.numeric(subj_n)),
         dia_id     = paste0("ASD_", family_int, "_", subj_int),
         geno_id    = str_replace_all(genomics, "-", "_"))

# ASD-only cis-pQTL at peptide level (column 'gene' holds the peptide ID <UniProt>_<sequence>)
pqtl_final <- read_excel(pqtl_cis_file) %>%
  transmute(GeneSymbol, snp_id = snps, region = "Cis",
            p_cond = pvalue, beta_cond = beta, peptide_id = gene)

# Lead variant per protein (smallest P over all peptides of the protein)
pqtl_leads <- pqtl_final %>%
  filter(GeneSymbol %in% risk_genes) %>%
  arrange(GeneSymbol,
          factor(region, levels = c("Cis", "Trans")),
          p_cond) %>%
  group_by(GeneSymbol) %>%
  slice(1) %>%
  ungroup() %>%
  rename(snps = snp_id, beta = beta_cond, pvalue = p_cond) %>%
  mutate(pep_sequence = str_extract(peptide_id, "(?<=_)[A-Z]+$"))

cat("Lead pQTL SNPs (cis-prioritized):\n")
print(pqtl_leads %>% select(GeneSymbol, snps, region, beta, pvalue, peptide_id, pep_sequence))

cat("\nRisk proteins after shared filtering:\n")
for (g in risk_genes) {
  uid <- risk_uniprot[g]
  found <- uid %in% dia_processed$Protein.Ids
  cat(sprintf("  %s (%s): %s\n", g, uid, ifelse(found, "PRESENT", "FILTERED OUT")))
}


# ==============================================================================
# 4. pQTL replication: TMT discovery peptide and DIA peptide/protein
# ==============================================================================

cat("\n=== pQTL replication ===\n")

geno_case <- read_tsv(geno_case_file, show_col_types = FALSE)
geno_ctrl <- read_tsv(geno_ctrl_file, show_col_types = FALSE)

# TMT sample grouping by column name
tmt_asd_cols <- grep("_03$", colnames(tmt_pep), value = TRUE)
tmt_td_cols  <- grep("^[A-G]\\d$", colnames(tmt_pep), value = TRUE)

# Genotyped TMT individuals (ASD + TD)
geno_mapping_tmt <- id_mapping %>%
  filter(geno_id %in% c(colnames(geno_case), colnames(geno_ctrl)),
         tmt_id %in% c(tmt_asd_cols, tmt_td_cols)) %>%
  mutate(group = case_when(
    grepl(SUBJECT_LABEL_ASD, subject) ~ "ASD",
    grepl(SUBJECT_LABEL_TD,  subject) ~ "TD",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group)) %>%
  select(tmt_id, geno_id, group)

# Genotyped DIA individuals (ASD, case genotype file)
geno_mapping_dia <- id_mapping %>%
  filter(geno_id %in% colnames(geno_case), dia_id %in% dia_sample_ids) %>%
  select(dia_id, geno_id)

cat(sprintf("  TMT matched: %d (ASD=%d, TD=%d)\n",
            nrow(geno_mapping_tmt),
            sum(geno_mapping_tmt$group == "ASD"),
            sum(geno_mapping_tmt$group == "TD")))
cat(sprintf("  DIA matched: %d\n", nrow(geno_mapping_dia)))

# Variant-peptide pair per protein: the best (smallest P) cis-pQTL pair whose peptide is
# quantified in DIA (DIA tested at peptide level); if no pQTL peptide of the protein is
# quantified in DIA, the protein's lead pair is used and DIA is tested at protein level.
dia_pep_avail <- dia_pep_processed %>%
  distinct(Protein.Ids, Stripped.Sequence)

best_pairs <- map_dfr(risk_genes, function(g) {
  uid <- risk_uniprot[g]
  dia_seqs <- dia_pep_avail %>% filter(Protein.Ids == uid) %>% pull(Stripped.Sequence)

  matched <- pqtl_final %>%
    filter(GeneSymbol == g) %>%
    mutate(pep_seq = str_extract(peptide_id, "(?<=_)[A-Z]+$")) %>%
    filter(pep_seq %in% dia_seqs) %>%
    arrange(p_cond)

  if (nrow(matched) > 0) {
    best <- matched %>% slice(1)
    tibble(gene = g, snp = best$snp_id, region = best$region,
           peptide_id = best$peptide_id, pep_seq = best$pep_seq,
           tmt_beta = best$beta_cond, tmt_p = best$p_cond,
           dia_level = "peptide")
  } else {
    lead <- pqtl_leads %>% filter(GeneSymbol == g)
    tibble(gene = g, snp = lead$snps, region = lead$region,
           peptide_id = lead$peptide_id, pep_seq = lead$pep_sequence,
           tmt_beta = lead$beta, tmt_p = lead$pvalue,
           dia_level = "protein")
  }
})

cat("\nBest pairs per gene:\n")
print(best_pairs %>% select(gene, snp, region, peptide_id, dia_level))

pqtl_rep <- list()

for (i in seq_len(nrow(best_pairs))) {
  gene       <- best_pairs$gene[i]
  snp_id     <- best_pairs$snp[i]
  snp_region <- best_pairs$region[i]
  pep_id     <- best_pairs$peptide_id[i]
  pep_seq    <- best_pairs$pep_seq[i]
  dia_level  <- best_pairs$dia_level[i]
  uniprot    <- risk_uniprot[gene]
  chr_pos    <- str_extract(snp_id, "^\\d+:\\d+")

  geno_row_c <- geno_case %>% filter(variantid == snp_id)
  geno_row_t <- geno_ctrl %>% filter(variantid == snp_id)

  # TMT: discovery peptide (z-scored), ASD + TD
  tmt_row <- tmt_pep_z %>% filter(Peptide == pep_id)

  tmt_df <- NULL
  if (nrow(tmt_row) > 0 & (nrow(geno_row_c) > 0 | nrow(geno_row_t) > 0)) {
    tmt_df <- map_dfr(seq_len(nrow(geno_mapping_tmt)), function(j) {
      gid <- geno_mapping_tmt$geno_id[j]
      tid <- geno_mapping_tmt$tmt_id[j]
      gval <- NA_integer_
      if (gid %in% colnames(geno_row_c)) gval <- as.integer(geno_row_c[[gid]])
      if (is.na(gval) && gid %in% colnames(geno_row_t)) gval <- as.integer(geno_row_t[[gid]])
      pval <- if (tid %in% colnames(tmt_row)) as.numeric(tmt_row[[tid]]) else NA_real_
      tibble(sample = tid, genotype = gval, abundance = pval,
             group = geno_mapping_tmt$group[j])
    }) %>%
      filter(!is.na(genotype), !is.na(abundance))
  }

  # DIA: genotyped ASD individuals, peptide or protein level
  dia_df <- NULL
  dia_row <- if (dia_level == "peptide") {
    dia_pep_processed %>% filter(Protein.Ids == uniprot, Stripped.Sequence == pep_seq)
  } else {
    dia_processed %>% filter(Protein.Ids == uniprot)
  }
  if (nrow(dia_row) > 0 & nrow(geno_row_c) > 0) {
    dia_df <- map_dfr(seq_len(nrow(geno_mapping_dia)), function(j) {
      tibble(sample    = geno_mapping_dia$dia_id[j],
             genotype  = as.integer(geno_row_c[[geno_mapping_dia$geno_id[j]]]),
             abundance = as.numeric(dia_row[[geno_mapping_dia$dia_id[j]]]))
    }) %>%
      filter(!is.na(genotype), !is.na(abundance))
  }

  # abundance ~ genotype, fitted only with >= 5 observations and >= 2 genotype classes;
  # otherwise beta and P are NA (in DIA, C3, POSTN, BTD and IGFBP5 have a single
  # genotype class because no genotyped ASD individual carries the minor allele).
  tmt_b <- NA; tmt_p_val <- NA; dia_b <- NA; dia_p_val <- NA

  if (!is.null(tmt_df) && nrow(tmt_df) >= 5 && length(unique(tmt_df$genotype)) >= 2) {
    fit <- lm(abundance ~ genotype, data = tmt_df)
    tmt_b <- coef(fit)["genotype"]
    tmt_p_val <- summary(fit)$coefficients["genotype","Pr(>|t|)"]
  }

  if (!is.null(dia_df) && nrow(dia_df) >= 5 && length(unique(dia_df$genotype)) >= 2) {
    fit <- lm(abundance ~ genotype, data = dia_df)
    dia_b <- coef(fit)["genotype"]
    dia_p_val <- summary(fit)$coefficients["genotype","Pr(>|t|)"]
  }

  concordant <- if (!is.na(tmt_b) & !is.na(dia_b)) sign(tmt_b) == sign(dia_b) else NA

  pqtl_rep[[gene]] <- tibble(
    gene = gene, snp = snp_id, chr_pos = chr_pos, region = snp_region,
    peptide = pep_id, dia_level = dia_level,
    tmt_beta = tmt_b, tmt_p = tmt_p_val, n_tmt = if (!is.null(tmt_df)) nrow(tmt_df) else 0,
    dia_beta = dia_b, dia_p = dia_p_val, n_dia = if (!is.null(dia_df)) nrow(dia_df) else 0,
    concordant = concordant)

  cat(sprintf("  %s [%s] %s (%s): TMT b=%.3f p=%.2e (N=%d) | DIA[%s] b=%.3f p=%.2e (N=%d) | %s\n",
              gene, snp_region, pep_id, chr_pos,
              tmt_b, tmt_p_val, if (!is.null(tmt_df)) nrow(tmt_df) else 0,
              dia_level, dia_b, dia_p_val,
              if (!is.null(dia_df)) nrow(dia_df) else 0,
              ifelse(is.na(concordant), "N/A", ifelse(concordant, "concordant", "discordant"))))
}

pqtl_rep_df <- bind_rows(pqtl_rep)
write_csv(pqtl_rep_df, file.path(out_dir, "stats_pqtl_replication.csv"))

cat(sprintf("  DIA genotype test run: %d / %d proteins; same direction as TMT: %d / %d\n",
            sum(!is.na(pqtl_rep_df$dia_beta)), nrow(pqtl_rep_df),
            sum(pqtl_rep_df$concordant, na.rm = TRUE), sum(!is.na(pqtl_rep_df$concordant))))


# ==============================================================================
# 5. Cross-platform concordance: TMT log2FC vs DIA difference in mean z
# ==============================================================================

cat("\n=== Cross-platform concordance ===\n")

asd_ids  <- dia_meta %>% filter(Group == "ASD") %>% pull(Sample_id)
ctrl_ids <- dia_meta %>% filter(Group == "Control") %>% pull(Sample_id)

# DIA ASD - Control difference in mean protein z-score
dia_fc_all <- dia_processed %>%
  rowwise() %>%
  mutate(
    asd_m  = mean(c_across(any_of(asd_ids)), na.rm = TRUE),
    ctrl_m = mean(c_across(any_of(ctrl_ids)), na.rm = TRUE),
    dia_fc = asd_m - ctrl_m
  ) %>%
  ungroup() %>%
  select(Protein.Ids, dia_fc)

# Proteins on both platforms, matched by gene symbol (first DIA protein group per gene)
dep_conc <- tmt_dep %>%
  select(Gene, tmt_fc = ASD_nonASD_fc, tmt_p = ASD_nonASD_p_val) %>%
  mutate(tmt_log2fc = log2(tmt_fc)) %>%
  inner_join(pr_gene_map %>% rename(Gene = Genes) %>%
               distinct(Gene, .keep_all = TRUE), by = "Gene") %>%
  inner_join(dia_fc_all, by = "Protein.Ids") %>%
  mutate(is_risk = Gene %in% risk_genes) %>%
  filter(!is.na(dia_fc), !is.na(tmt_log2fc))

fc_r <- cor(dep_conc$tmt_log2fc, dep_conc$dia_fc, use = "complete.obs")
cat(sprintf("N = %d, Pearson r = %.3f\n", nrow(dep_conc), fc_r))

write_csv(dep_conc, file.path(out_dir, "stats_dep_concordance.csv"))


# ==============================================================================
# 6. Pathway-level replication: fraction of DEPs per pathway
# ==============================================================================

cat("\n=== Pathway replication ===\n")

pathway_genes <- list(
  "Complement & Immune Activation" =
    c("C1RL","C1R","C1S","C1QB","C3","C6","C5","C7","C8A","C8B",
      "MBL2","MASP1","MASP2","CFH","CFHR3","CFP","CD46","SERPING1"),
  "ECM Remodeling" =
    c("COL6A1","COL1A1","POSTN","ABI3BP","COMP","LUM","MMP2",
      "TNXB","RECK","VASN","ACAN","FN1"),
  "Integrin Signaling" =
    c("ITGAV","ITGA1","ITGA2","ITGA5","ITGB1","ITGB3","ITGB5",
      "DDR1","DDR2","SDC1","SDC4","LRP1"),
  "IGF1 Signaling"    = c("IGF1","IGFBP5","IGFALS","IGFBP3"),
  "Coagulation"       = c("F2","F5","F11","F13A1","VTN","PLG","SERPIND1"),
  "Synaptic Adhesion" = c("NCAM2","OLFM1","NCAM1","NRCAM")
)

# TMT DEP: P < 0.05 in the TMT differential abundance analysis
tmt_dep_sig   <- tmt_dep %>% filter(ASD_nonASD_p_val < 0.05) %>% pull(Gene)
dia_all_genes <- pr_gene_map$Genes

# DIA DEP: Welch t-test ASD vs Control on protein z-scores, P < 0.05
dia_pvals <- dia_processed %>%
  inner_join(pr_gene_map %>% distinct(Protein.Ids, .keep_all = TRUE), by = "Protein.Ids") %>%
  rowwise() %>%
  mutate(
    p = tryCatch({
      a <- c_across(any_of(asd_ids)); b <- c_across(any_of(ctrl_ids))
      if (sum(!is.na(a)) >= 3 & sum(!is.na(b)) >= 3) t.test(a, b)$p.value else NA_real_
    }, error = function(e) NA_real_)
  ) %>%
  ungroup() %>%
  select(Protein.Ids, Genes, p)

dia_dep_sig <- dia_pvals %>% filter(!is.na(p), p < 0.05) %>% pull(Genes)

pw_res <- map_dfr(names(pathway_genes), function(pw) {
  g <- pathway_genes[[pw]]
  tibble(pathway = pw,
         tmt_frac = sum(g %in% tmt_dep_sig) / max(sum(g %in% tmt_dep$Gene), 1),
         dia_frac = sum(g %in% dia_dep_sig) / max(sum(g %in% dia_all_genes), 1),
         tmt_n_dep = sum(g %in% tmt_dep_sig), dia_n_dep = sum(g %in% dia_dep_sig),
         tmt_n_meas = sum(g %in% tmt_dep$Gene), dia_n_meas = sum(g %in% dia_all_genes))
})

print(pw_res)
write_csv(pw_res, file.path(out_dir, "stats_pathway_replication.csv"))


# ==============================================================================
# 7. DIA five-group summaries for the 12 proteins
# ==============================================================================

cat("\n=== DIA group summaries ===\n")

# Protein-level z-scores of the 12 proteins in all DIA samples with a group label
risk_expr <- dia_processed %>%
  filter(Protein.Ids %in% risk_uniprot) %>%
  mutate(Gene = uniprot_to_gene[Protein.Ids]) %>%
  select(Gene, all_of(dia_sample_ids)) %>%
  pivot_longer(-Gene, names_to = "Sample_id", values_to = "abundance") %>%
  left_join(dia_meta, by = "Sample_id") %>%
  filter(!is.na(Group), !is.na(abundance))

group_summary <- risk_expr %>%
  group_by(protein = Gene, Group) %>%
  summarise(n = n(), mean = mean(abundance), se = sd(abundance) / sqrt(n()),
            .groups = "drop")

# ASD vs Control (TD): Welch t-test when both groups have > 2 values; diff = ASD - Control
asd_vs_ctrl <- risk_expr %>%
  group_by(protein = Gene) %>%
  group_modify(function(d, key) {
    a <- d$abundance[d$Group == "ASD"]
    r <- d$abundance[d$Group == "Control"]
    if (length(a) > 2 && length(r) > 2) {
      tt <- t.test(a, r)
      tibble(diff = unname(diff(rev(tt$estimate))), p = tt$p.value,
             n_ASD = length(a), n_Ctrl = length(r))
    } else {
      tibble(diff = NA_real_, p = NA_real_, n_ASD = length(a), n_Ctrl = length(r))
    }
  }) %>%
  ungroup()

cat(sprintf("DIA ASD vs Control, nominal P < 0.05: %d / %d\n",
            sum(asd_vs_ctrl$p < 0.05, na.rm = TRUE), sum(!is.na(asd_vs_ctrl$p))))

write_tsv(group_summary, file.path(dat_dir, "dia_group_summary.tsv"))
write_tsv(asd_vs_ctrl,   file.path(dat_dir, "dia_ASD_vs_Control.tsv"))

cat("\n=== COMPLETE ===\n")

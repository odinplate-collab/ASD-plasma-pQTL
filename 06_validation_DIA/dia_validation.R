################################################################################
## ASD Plasma pQTL Validation Analysis
## DIA validation cohort: replication analyses
##
## Main Panels (protein level):
##   a. pQTL Replication
##   b. Expression Validation
##   c. Cross-platform Concordance
##   d. Genome-wide DEP Concordance
##   e. Pathway-level Replication
##
## Supplementary (peptide level):
##   s1. pQTL Replication (top 3 SNPs per gene, all peptides)
##   s2. Expression Validation (all pQTL peptides per gene)
################################################################################

# ==============================================================================
# SECTION 0: SETUP
# ==============================================================================

library(tidyverse)
library(readxl)
library(patchwork)
library(ggrepel)
library(preprocessCore)

base_dir   <- file.path(Sys.getenv("ASD_ROOT", "."), "Validation")
output_dir <- file.path(base_dir, "plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

dia_pg_file      <- file.path(base_dir, "ASD_459_report.pg_matrix.tsv")
dia_pr_file      <- file.path(base_dir, "ASD_459_report.pr_matrix.tsv")
dia_meta_file    <- file.path(base_dir, "20260311_ASD_DIA_metadata.xlsx")
tmt_pep_file     <- file.path(base_dir, "peptides.tsv")
tmt_dep_file     <- file.path(base_dir, "Supplementary_table2.xlsx")
pqtl_final_file  <- file.path(base_dir, "ASD_pQTL_final_outputs.xlsx")
sample_info_file <- file.path(base_dir, "Total_ASD_sample_info.xlsx")
geno_case_file   <- file.path(base_dir, "ASD_pQTL_imputated_case_rm29903_n_alt.tsv")
geno_ctrl_file <- file.path(base_dir, "ASD_pQTL_imputated_ctrl_rm29903_n_alt.tsv")

risk_genes   <- c("C3","AHSG","MBL2","VASN","COL6A1","POSTN","VNN1","CFHR3","C1RL")
risk_uniprot <- c(C3="P01024", AHSG="P02765", MBL2="P11226", VASN="Q6EMK4",
                  COL6A1="P12109", POSTN="Q15063", VNN1="O95497",
                  CFHR3="Q02985", C1RL="Q9NZP8")
uniprot_to_gene <- setNames(names(risk_uniprot), risk_uniprot)

theme_ng <- function(base_size = 7) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text             = element_text(family = "Helvetica", colour = "black"),
      axis.text        = element_text(size = rel(0.9), colour = "black"),
      axis.title       = element_text(size = rel(1), face = "plain"),
      axis.line        = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
      strip.background = element_blank(),
      strip.text       = element_text(size = rel(1), face = "bold"),
      legend.title     = element_text(size = rel(0.9)),
      legend.text      = element_text(size = rel(0.85)),
      legend.key.size  = unit(3, "mm"),
      plot.title       = element_text(size = rel(1.1), face = "bold", hjust = 0),
      plot.margin      = margin(4, 4, 4, 4),
      panel.grid       = element_blank()
    )
}

pal_group <- c(ASD = "#C44E52", Father = "#4C72B0", Mother = "#8DA0CB",
               Sibling = "#E6AB02", Control = "#55A868")
pal_geno  <- c("0" = "#4C72B0", "1" = "#CCB974", "2" = "#C44E52")


# ==============================================================================
# SECTION 1: LOAD & PREPROCESS
# ==============================================================================

cat("Loading data...\n")

# --- 1a. DIA protein matrix (pg level) ---
dia_raw <- read_tsv(dia_pg_file, show_col_types = FALSE)
dia_sample_cols <- grep("ASD.*\\.raw$", colnames(dia_raw), value = TRUE)
dia_sample_ids  <- str_extract(dia_sample_cols, "ASD_\\d+(?:_\\d+)?(?=_2ug)")

ok <- !is.na(dia_sample_ids) & !duplicated(dia_sample_ids)
if (any(!ok)) {
  cat(sprintf("WARNING: dropping %d problematic columns\n", sum(!ok)))
  dia_sample_cols <- dia_sample_cols[ok]
  dia_sample_ids  <- dia_sample_ids[ok]
}

dia_clean <- dia_raw %>% select(Protein.Ids, all_of(dia_sample_cols))
colnames(dia_clean) <- c("Protein.Ids", dia_sample_ids)

# --- 1b. DIA metadata ---
dia_meta <- read_excel(dia_meta_file) %>%
  mutate(Family = as.character(as.integer(Family)),
         Group  = factor(Group, levels = c("ASD","Father","Mother","Sibling","Control")))

# --- Identify Control vs Rest indices ---
ctrl_samples <- dia_meta %>% filter(Group == "Control") %>% pull(Sample_id)
ctrl_idx <- which(dia_sample_ids %in% ctrl_samples)
rest_idx <- which(!dia_sample_ids %in% ctrl_samples)

# --- Shared protein filtering + quantile normalization ---
dia_mat <- dia_clean %>% select(all_of(dia_sample_ids)) %>% as.matrix()
dia_mat[dia_mat == 0] <- NA

# Keep proteins detected in >= 70% of BOTH Control and Rest
detect_ctrl <- apply(dia_mat[, ctrl_idx], 1, function(x) mean(!is.na(x) & x > 0))
detect_rest <- apply(dia_mat[, rest_idx], 1, function(x) mean(!is.na(x) & x > 0))
shared <- detect_ctrl >= 0.7 & detect_rest >= 0.7

cat(sprintf("  Shared protein filtering: %d -> %d proteins (>= 70%% in both groups)\n",
            nrow(dia_mat), sum(shared)))

dia_mat   <- dia_mat[shared, ]
dia_clean <- dia_clean[shared, ]

# Log2 transform
dia_log2 <- log2(dia_mat)

# Quantile normalization
dia_qnorm <- normalize.quantiles(dia_log2)
colnames(dia_qnorm) <- dia_sample_ids
rownames(dia_qnorm) <- rownames(dia_log2)

# Z-score per protein
dia_z <- t(apply(dia_qnorm, 1, function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA, length(x)))
  (x - m) / s
}))
colnames(dia_z) <- dia_sample_ids

dia_processed <- tibble(Protein.Ids = dia_clean$Protein.Ids) %>%
  bind_cols(as_tibble(dia_z))

cat(sprintf("  DIA protein: %d proteins x %d samples (shared + quantile norm + z-score)\n",
            nrow(dia_processed), length(dia_sample_ids)))

# --- 1c. TMT DEP results ---
tmt_dep <- read_excel(tmt_dep_file, sheet = 1)

# --- 1d. TMT peptide matrix ---
tmt_pep <- read_tsv(tmt_pep_file, show_col_types = FALSE)

# convert TMT peptide values to z-scores
tmt_sample_cols <- setdiff(colnames(tmt_pep), c("Peptide","Protein","Gene"))

tmt_pep_z <- tmt_pep %>%
  rowwise() %>%
  mutate(
    row_mean = mean(c_across(all_of(tmt_sample_cols)), na.rm = TRUE),
    row_sd   = sd(c_across(all_of(tmt_sample_cols)), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(across(all_of(tmt_sample_cols), ~( .x - row_mean) / row_sd)) %>%
  select(-row_mean, -row_sd)


tmt_gene <- tmt_pep %>%
  filter(!is.na(Gene)) %>%
  group_by(Gene) %>%
  summarise(across(all_of(tmt_sample_cols), ~mean(.x, na.rm = TRUE)), .groups = "drop")

# --- 1e. Sample ID mapping ---
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

xplatform <- id_mapping %>%
  filter(tmt_id %in% tmt_sample_cols, dia_id %in% dia_sample_ids)

# --- 1f. UniProt-to-Gene mapping ---
pr_gene_map <- read_tsv(dia_pr_file, col_types = cols(),
                        col_select = c("Protein.Ids","Genes")) %>%
  distinct(Protein.Ids, Genes) %>%
  filter(!is.na(Genes), Genes != "")

# --- 1g. pQTL lead SNPs (cis + trans, from final_outputs) ---
pqtl_final <- read_excel(pqtl_final_file, sheet = "final_cond")

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

# --- Check risk proteins survived filtering ---
cat("\nRisk proteins after shared filtering:\n")
for (g in risk_genes) {
  uid <- risk_uniprot[g]
  found <- uid %in% dia_processed$Protein.Ids
  cat(sprintf("  %s (%s): %s\n", g, uid, ifelse(found, "PRESENT", "FILTERED OUT")))
}

cat("\nData loading complete.\n\n")

# ==============================================================================
# SECTION 2: ANALYSIS 1 — pQTL REPLICATION (TMT peptide + DIA peptide/protein)
# ==============================================================================

cat("=== Analysis 1: pQTL Replication ===\n")

geno_case <- read_tsv(geno_case_file, show_col_types = FALSE)
geno_ctrl <- read_tsv(geno_ctrl_file, show_col_types = FALSE)

# TMT peptide z-score
tmt_pep_z <- tmt_pep %>%
  rowwise() %>%
  mutate(
    row_mean = mean(c_across(all_of(tmt_sample_cols)), na.rm = TRUE),
    row_sd   = sd(c_across(all_of(tmt_sample_cols)), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(across(all_of(tmt_sample_cols), ~(.x - row_mean) / row_sd)) %>%
  select(-row_mean, -row_sd)

# TMT sample grouping
tmt_asd_cols <- grep("_03$", colnames(tmt_pep), value = TRUE)
tmt_td_cols  <- grep("^[A-G]\\d$", colnames(tmt_pep), value = TRUE)

# Genotype mapping for TMT (ASD + TD)
geno_mapping_tmt <- id_mapping %>%
  filter(geno_id %in% c(colnames(geno_case), colnames(geno_ctrl)),
         tmt_id %in% c(tmt_asd_cols, tmt_td_cols)) %>%
  mutate(group = case_when(
    grepl("자폐", subject) ~ "ASD",   # "자폐" = autism (label used in the source metadata)
    grepl("정상", subject) ~ "TD",    # "정상" = normal/typical (label used in the source metadata)
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group)) %>%
  select(tmt_id, geno_id, group)

# Genotype mapping for DIA
geno_mapping_dia <- id_mapping %>%
  filter(geno_id %in% colnames(geno_case), dia_id %in% dia_sample_ids) %>%
  select(dia_id, geno_id)

cat(sprintf("  TMT matched: %d (ASD=%d, TD=%d)\n",
            nrow(geno_mapping_tmt),
            sum(geno_mapping_tmt$group == "ASD"),
            sum(geno_mapping_tmt$group == "TD")))
cat(sprintf("  DIA matched: %d\n", nrow(geno_mapping_dia)))

# --- Load DIA peptide data if not yet loaded ---
if (!exists("dia_pep_processed")) {
  cat("Loading DIA precursor matrix...\n")
  dia_pr_raw <- read_tsv(dia_pr_file, show_col_types = FALSE)
  dia_pr_sample_cols <- grep("^ASD_", colnames(dia_pr_raw), value = TRUE)

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

  dia_pep <- dia_pep[pep_shared, ]
  pep_mat <- pep_mat[pep_shared, ]

  pep_log2 <- log2(pep_mat)
  pep_qnorm <- normalize.quantiles(pep_log2)
  colnames(pep_qnorm) <- dia_pr_sample_cols

  pep_z <- t(apply(pep_qnorm, 1, function(x) {
    m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(NA, length(x)))
    (x - m) / s
  }))
  colnames(pep_z) <- dia_pr_sample_cols

  dia_pep_processed <- dia_pep %>%
    select(Protein.Ids, Stripped.Sequence) %>%
    bind_cols(as_tibble(pep_z))
}

# --- Build best DIA peptide map per gene ---
# For each gene: find pQTL SNP-peptide pair where peptide exists in DIA
# If none, fall back to protein level

# First: check all pQTL pairs for DIA peptide overlap
dia_pep_avail <- dia_pep_processed %>%
  distinct(Protein.Ids, Stripped.Sequence)

best_pairs <- map_dfr(risk_genes, function(g) {
  uid <- risk_uniprot[g]

  # All DIA peptides for this protein
  dia_seqs <- dia_pep_avail %>% filter(Protein.Ids == uid) %>% pull(Stripped.Sequence)

  # Find pQTL pairs matching DIA peptides
  matched <- pqtl_final %>%
    filter(GeneSymbol == g) %>%
    mutate(pep_seq = str_extract(peptide_id, "(?<=_)[A-Z]+$")) %>%
    filter(pep_seq %in% dia_seqs) %>%
    arrange(p_cond)

  if (nrow(matched) > 0) {
    # Best pQTL pair with DIA peptide available
    best <- matched %>% slice(1)
    tibble(gene = g, snp = best$snp_id, region = best$region,
           peptide_id = best$peptide_id, pep_seq = best$pep_seq,
           tmt_beta = best$beta_cond, tmt_p = best$p_cond,
           dia_level = "peptide")
  } else {
    # Fall back to protein level — use original lead SNP
    lead <- pqtl_leads %>% filter(GeneSymbol == g)
    tibble(gene = g, snp = lead$snps, region = lead$region,
           peptide_id = lead$peptide_id, pep_seq = lead$pep_sequence,
           tmt_beta = lead$beta, tmt_p = lead$pvalue,
           dia_level = "protein")
  }
})

cat("\nBest pairs per gene:\n")
print(best_pairs %>% select(gene, snp, region, peptide_id, dia_level))

# --- Main loop ---
pqtl_rep   <- list()
pqtl_plots <- list()

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

  # --- TMT: always peptide level (z-scored) ---
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
      filter(!is.na(genotype), !is.na(abundance)) %>%
      mutate(genotype_f = factor(genotype), platform = "TMT")
  }

  # --- DIA: peptide or protein level ---
  dia_df <- NULL
  if (dia_level == "peptide") {
    dia_pep_row <- dia_pep_processed %>%
      filter(Protein.Ids == uniprot, Stripped.Sequence == pep_seq)
    if (nrow(dia_pep_row) > 0 & nrow(geno_row_c) > 0) {
      dia_df <- map_dfr(seq_len(nrow(geno_mapping_dia)), function(j) {
        tibble(sample   = geno_mapping_dia$dia_id[j],
               genotype = as.integer(geno_row_c[[geno_mapping_dia$geno_id[j]]]),
               abundance = as.numeric(dia_pep_row[[geno_mapping_dia$dia_id[j]]]))
      }) %>%
        filter(!is.na(genotype), !is.na(abundance)) %>%
        mutate(genotype_f = factor(genotype), platform = "DIA", group = "ASD")
    }
  } else {
    # Protein level
    dia_row <- dia_processed %>% filter(Protein.Ids == uniprot)
    if (nrow(dia_row) > 0 & nrow(geno_row_c) > 0) {
      dia_df <- map_dfr(seq_len(nrow(geno_mapping_dia)), function(j) {
        tibble(sample   = geno_mapping_dia$dia_id[j],
               genotype = as.integer(geno_row_c[[geno_mapping_dia$geno_id[j]]]),
               abundance = as.numeric(dia_row[[geno_mapping_dia$dia_id[j]]]))
      }) %>%
        filter(!is.na(genotype), !is.na(abundance)) %>%
        mutate(genotype_f = factor(genotype), platform = "DIA", group = "ASD")
    }
  }

  # --- Stats ---
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
    tmt_beta = tmt_b, tmt_p = tmt_p_val, n_tmt = if(!is.null(tmt_df)) nrow(tmt_df) else 0,
    dia_beta = dia_b, dia_p = dia_p_val, n_dia = if(!is.null(dia_df)) nrow(dia_df) else 0,
    concordant = concordant)

  combined <- bind_rows(tmt_df, dia_df)
  if (!is.null(combined) && nrow(combined) > 0) {
    pqtl_plots[[gene]] <- combined %>%
      mutate(gene_label = gene, chr_pos = chr_pos, dia_level = dia_level)
  }

  cat(sprintf("  %s [%s] %s (%s): TMT b=%.3f p=%.2e (N=%d) | DIA[%s] b=%.3f p=%.2e (N=%d) | %s\n",
              gene, snp_region, pep_id, chr_pos,
              tmt_b, tmt_p_val, if(!is.null(tmt_df)) nrow(tmt_df) else 0,
              dia_level, dia_b, dia_p_val,
              if(!is.null(dia_df)) nrow(dia_df) else 0,
              ifelse(is.na(concordant), "N/A", ifelse(concordant, "concordant", "discordant"))))
}

pqtl_rep_df  <- bind_rows(pqtl_rep)
pqtl_plot_df <- bind_rows(pqtl_plots)

# --- Panel A ---
plot_order <- intersect(c("MBL2","C1RL","C3","AHSG","CFHR3","POSTN","VASN","COL6A1","VNN1"),
                        unique(pqtl_plot_df$gene_label))

pqtl_plot_df <- pqtl_plot_df %>%
  filter(gene_label %in% plot_order) %>%
  mutate(
    gene_label = factor(gene_label, levels = plot_order),
    platform = factor(platform, levels = c("TMT","DIA")),
    group = factor(group, levels = c("ASD","TD")),
    panel_x = interaction(platform, group, sep = "\n")
  )

x_levels <- c("TMT\nASD", "TMT\nTD", "DIA\nASD")
pqtl_plot_df$panel_x <- factor(pqtl_plot_df$panel_x, levels = x_levels)

# Facet labels: Gene + chr:pos + DIA level
facet_labels <- pqtl_rep_df %>%
  filter(gene %in% plot_order) %>%
  mutate(facet = sprintf("%s\n%s\nDIA: %s", gene, chr_pos, dia_level))
facet_map <- setNames(facet_labels$facet, facet_labels$gene)

pqtl_plot_df$facet <- facet_map[as.character(pqtl_plot_df$gene_label)]
pqtl_plot_df$facet <- factor(pqtl_plot_df$facet, levels = facet_map[plot_order])

stat_labels <- pqtl_rep_df %>%
  filter(gene %in% plot_order) %>%
  mutate(
    label = sprintf("TMT: b=%.2f, p=%.1e\nDIA: b=%.2f, p=%.1e",
                    tmt_beta, tmt_p, dia_beta, dia_p),
    facet = factor(facet_map[gene], levels = facet_map[plot_order])
  )

panel_a <- ggplot(pqtl_plot_df, aes(x = panel_x, y = abundance)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.3,
               fill = "grey95", colour = "black") +
  geom_jitter(aes(colour = genotype_f), width = 0.15, size = 1.2,
              alpha = 0.75, shape = 16) +
  geom_vline(xintercept = 2.5, linetype = "dashed", linewidth = 0.2,
             colour = "grey50") +
  geom_text(data = stat_labels, aes(label = label),
            x = Inf, y = Inf, hjust = 1.02, vjust = 1.1,
            size = 1.2, colour = "grey30", lineheight = 0.8,
            inherit.aes = FALSE) +
  facet_wrap(~facet, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = c("0" = "#4C72B0", "1" = "#E6AB02", "2" = "#C44E52"),
                      name = "Genotype") +
  labs(x = NULL, y = "Abundance (z-score)") +
  theme_ng() +
  theme(legend.position = "right", legend.key.size = unit(2.5, "mm"),
        axis.text.x = element_text(size = 5))

cat("Panel A done.\n")
# ==============================================================================
# SECTION 3: ANALYSIS 2 — EXPRESSION VALIDATION (protein level)
# ==============================================================================

cat("\n=== Analysis 2: Expression Validation (protein level) ===\n")

risk_expr <- dia_processed %>%
  filter(Protein.Ids %in% risk_uniprot) %>%
  mutate(Gene = uniprot_to_gene[Protein.Ids]) %>%
  select(Gene, all_of(dia_sample_ids)) %>%
  pivot_longer(-Gene, names_to = "Sample_id", values_to = "abundance") %>%
  left_join(dia_meta, by = "Sample_id") %>%
  filter(!is.na(Group), !is.na(abundance))

expr_stats <- risk_expr %>%
  group_by(Gene) %>%
  summarise(
    asd_mean     = mean(abundance[Group == "ASD"], na.rm = TRUE),
    non_asd_mean = mean(abundance[Group != "ASD"], na.rm = TRUE),
    p_wilcox     = tryCatch(wilcox.test(abundance[Group == "ASD"],
                                         abundance[Group != "ASD"])$p.value,
                            error = function(e) NA),
    n_asd  = sum(Group == "ASD" & !is.na(abundance)),
    n_rest = sum(Group != "ASD" & !is.na(abundance)),
    .groups = "drop"
  ) %>%
  mutate(dia_direction = ifelse(asd_mean > non_asd_mean, "Up", "Down")) %>%
  left_join(
    tmt_dep %>% filter(Gene %in% risk_genes) %>%
      transmute(Gene, tmt_fc = ASD_nonASD_fc,
                tmt_direction = ifelse(ASD_nonASD_fc > 1, "Up", "Down")),
    by = "Gene"
  ) %>%
  mutate(concordant = dia_direction == tmt_direction)

cat("Expression validation:\n")
print(expr_stats %>% select(Gene, dia_direction, p_wilcox, tmt_direction, concordant))

# colours matched to the discovery-cohort group plots
pal_group <- c(Father = "#4472C4", Mother = "#B07AA1", Sibling = "#E6AB02",
               Control = "#59A14F", ASD = "#E15759")

# group order matched to the discovery-cohort group plots
risk_expr$Group <- factor(risk_expr$Group,
                           levels = c("Father","Mother","Sibling","Control","ASD"))

panel_b <- ggplot(risk_expr, aes(x = Group, y = abundance, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, linewidth = 0.3,
               alpha = 0.3) +
  geom_jitter(width = 0.18, size = 0.6, alpha = 0.6, shape = 16) +
  facet_wrap(~Gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal_group) +
  scale_colour_manual(values = pal_group) +
  labs(x = NULL, y = "Protein abundance (z-score)") +
  theme_ng() +
  theme(legend.position = "bottom",
        legend.margin = margin(-2, 0, 0, 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# ==============================================================================
# SECTION 4: ANALYSIS 3 — CROSS-PLATFORM CONCORDANCE
# ==============================================================================

cat("\n=== Analysis 3: Cross-platform Concordance ===\n")

xplat_data <- list()
for (gene in risk_genes) {
  uniprot <- risk_uniprot[gene]
  dia_row <- dia_processed %>% filter(Protein.Ids == uniprot)
  tmt_row <- tmt_dep %>% filter(Gene == gene)
  if (nrow(dia_row) == 0 | nrow(tmt_row) == 0) next
  for (k in seq_len(nrow(xplatform))) {
    t_id <- xplatform$tmt_id[k]; d_id <- xplatform$dia_id[k]
    if (!(t_id %in% colnames(tmt_row)) | !(d_id %in% colnames(dia_row))) next
    xplat_data[[length(xplat_data)+1]] <- tibble(
      gene = gene, sample = d_id,
      tmt = as.numeric(tmt_row[[t_id]]),
      dia = as.numeric(dia_row[[d_id]]))
  }
}
xplat_df <- bind_rows(xplat_data) %>% filter(!is.na(tmt), !is.na(dia))

xplat_cor <- xplat_df %>%
  group_by(gene) %>%
  summarise(r = cor(tmt, dia, use = "complete.obs"), n = n(), .groups = "drop")

overall_r <- cor(xplat_df$tmt, xplat_df$dia, use = "complete.obs")
cat(sprintf("Overall r = %.3f (N pairs = %d)\n", overall_r, nrow(xplat_df)))
print(xplat_cor)

cor_labels <- xplat_cor %>%
  mutate(label = sprintf("r = %.2f", r), gene = factor(gene, levels = risk_genes))
xplat_df$gene <- factor(xplat_df$gene, levels = risk_genes)

panel_c <- ggplot(xplat_df, aes(x = tmt, y = dia)) +
  geom_point(size = 0.7, alpha = 0.45, colour = "#4C72B0", shape = 16) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.35,
              colour = "#C44E52", linetype = "solid") +
  geom_text(data = cor_labels, aes(label = label),
            x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
            size = 1.8, colour = "grey25", inherit.aes = FALSE) +
  facet_wrap(~gene, scales = "free", nrow = 2) +
  labs(x = "TMT (z-score)", y = "DIA (z-score)") +
  theme_ng()

cat("Panel C done.\n")


# ==============================================================================
# SECTION 5: ANALYSIS 4 — GENOME-WIDE DEP CONCORDANCE
# ==============================================================================

cat("\n=== Analysis 4: DEP Concordance ===\n")

asd_ids  <- dia_meta %>% filter(Group == "ASD") %>% pull(Sample_id)
ctrl_ids <- dia_meta %>% filter(Group == "Control") %>% pull(Sample_id)

dia_fc_all <- dia_processed %>%
  rowwise() %>%
  mutate(
    asd_m  = mean(c_across(any_of(asd_ids)), na.rm = TRUE),
    ctrl_m = mean(c_across(any_of(ctrl_ids)), na.rm = TRUE),
    dia_fc = asd_m - ctrl_m
  ) %>%
  ungroup() %>%
  select(Protein.Ids, dia_fc)

dep_conc <- tmt_dep %>%
  select(Gene, tmt_fc = ASD_nonASD_fc, tmt_p = ASD_nonASD_p_val) %>%
  mutate(tmt_log2fc = log2(tmt_fc)) %>%
  inner_join(pr_gene_map %>% rename(Gene = Genes) %>%
               distinct(Gene, .keep_all = TRUE), by = "Gene") %>%
  inner_join(dia_fc_all, by = "Protein.Ids") %>%
  mutate(is_risk = Gene %in% risk_genes) %>%
  filter(!is.na(dia_fc), !is.na(tmt_log2fc))

fc_r <- cor(dep_conc$tmt_log2fc, dep_conc$dia_fc, use = "complete.obs")
n_conc <- sum((dep_conc$tmt_log2fc > 0 & dep_conc$dia_fc > 0) |
              (dep_conc$tmt_log2fc < 0 & dep_conc$dia_fc < 0))
cat(sprintf("N = %d, r = %.3f, direction concordance = %.1f%%\n",
            nrow(dep_conc), fc_r, n_conc/nrow(dep_conc)*100))

panel_d <- ggplot(dep_conc, aes(x = tmt_log2fc, y = dia_fc)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey70", linetype = "dashed") +
  geom_point(data = dep_conc %>% filter(!is_risk),
             colour = "grey65", size = 0.5, alpha = 0.45, shape = 16) +
  geom_point(data = dep_conc %>% filter(is_risk),
             colour = "#C44E52", size = 1.5, alpha = 0.85, shape = 16) +
  geom_smooth(method = "lm", colour = "black", se = TRUE,
              linewidth = 0.4, fill = "grey85", alpha = 0.25) +
  geom_text_repel(data = dep_conc %>% filter(is_risk),
                  aes(label = Gene), size = 1.9, fontface = "italic",
                  colour = "#C44E52", max.overlaps = 20,
                  segment.size = 0.15, segment.color = "grey50",
                  min.segment.length = 0) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.4,
           label = sprintf("r = %.2f\nn = %d", fc_r, nrow(dep_conc)),
           size = 2.2, colour = "black") +
  labs(x = expression(TMT~log[2]*"(FC)"),
       y = expression(DIA~Delta*"(z-score)")) +
  theme_ng()

cat("Panel D done.\n")


# ==============================================================================
# SECTION 6: ANALYSIS 5 — PATHWAY REPLICATION
# ==============================================================================

cat("\n=== Analysis 5: Pathway Replication ===\n")

pathway_genes <- list(
  "Complement &\nImmune Activation" =
    c("C1RL","C1R","C1S","C1QB","C3","C6","C5","C7","C8A","C8B",
      "MBL2","MASP1","MASP2","CFH","CFHR3","CFP","CD46","SERPING1"),
  "ECM\nRemodeling" =
    c("COL6A1","COL1A1","POSTN","ABI3BP","COMP","LUM","MMP2",
      "TNXB","RECK","VASN","ACAN","FN1"),
  "Integrin\nSignaling" =
    c("ITGAV","ITGA1","ITGA2","ITGA5","ITGB1","ITGB3","ITGB5",
      "DDR1","DDR2","SDC1","SDC4","LRP1"),
  "IGF1\nSignaling"     = c("IGF1","IGFBP5","IGFALS","IGFBP3"),
  "Coagulation"         = c("F2","F5","F11","F13A1","VTN","PLG","SERPIND1"),
  "Synaptic\nAdhesion"  = c("NCAM2","OLFM1","NCAM1","NRCAM")
)

tmt_dep_sig   <- tmt_dep %>% filter(ASD_nonASD_p_val < 0.05) %>% pull(Gene)
dia_all_genes <- pr_gene_map$Genes

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

cat("Pathway replication:\n"); print(pw_res)

pw_long <- pw_res %>%
  select(pathway, TMT = tmt_frac, DIA = dia_frac) %>%
  pivot_longer(-pathway, names_to = "Platform", values_to = "frac") %>%
  mutate(pathway = factor(pathway, levels = rev(unique(pw_res$pathway))))

panel_e <- ggplot(pw_long, aes(x = frac, y = pathway, fill = Platform)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5) +
  scale_fill_manual(values = c(TMT = "#4C72B0", DIA = "#C44E52")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Fraction DEP in pathway", y = NULL) +
  theme_ng() +
  theme(legend.position = c(0.85, 0.15), legend.background = element_blank())

cat("Panel E done.\n")


# ==============================================================================
# SECTION 7: SAVE MAIN PLOTS
# ==============================================================================

cat("\n=== Saving main plots ===\n")

n_pqtl_genes <- length(plot_order)
panel_a_width <- max(88, n_pqtl_genes * 17)

ggsave(file.path(output_dir, "panel_a_pqtl_replication.pdf"),
       panel_a, width = panel_a_width, height = 80, units = "mm", device = cairo_pdf)
ggsave(file.path(output_dir, "panel_b_expression_validation.pdf"),
       panel_b, width = 180, height = 30, units = "mm", device = cairo_pdf)
ggsave(file.path(output_dir, "panel_c_crossplatform.pdf"),
       panel_c, width = 180, height = 85, units = "mm", device = cairo_pdf)
ggsave(file.path(output_dir, "panel_d_dep_concordance.pdf"),
       panel_d, width = 50, height = 45, units = "mm", device = cairo_pdf)
ggsave(file.path(output_dir, "panel_e_pathway_replication.pdf"),
       panel_e, width = 85, height = 65, units = "mm", device = cairo_pdf)

composite <- (panel_a) /
  (panel_b) /
  (panel_c | panel_d) /
  (plot_spacer() | panel_e) +
  plot_layout(heights = c(1, 2.2, 2, 1.5)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(size = 8, face = "bold",
                                                        family = "Helvetica")))

ggsave(file.path(output_dir, "DIA_validation_summary.pdf"),
       composite, width = 180, height = 260, units = "mm", device = cairo_pdf)

write_csv(pqtl_rep_df,    file.path(output_dir, "stats_pqtl_replication.csv"))
write_csv(expr_stats,     file.path(output_dir, "stats_expression_validation.csv"))
write_csv(as.data.frame(xplat_cor), file.path(output_dir, "stats_crossplatform_cor.csv"))
write_csv(dep_conc,       file.path(output_dir, "stats_dep_concordance.csv"))
write_csv(pw_res,         file.path(output_dir, "stats_pathway_replication.csv"))

cat("Main plots saved.\n\n")


# ==============================================================================
# SECTION 8: SUPPLEMENTARY — PEPTIDE-LEVEL VALIDATION
# ==============================================================================

cat("=== Supplementary: Peptide-level validation ===\n")
# --- 8a. Load & preprocess DIA precursor matrix ---
cat("Loading DIA precursor matrix...\n")
dia_pr_raw <- read_tsv(dia_pr_file, show_col_types = FALSE)
dia_pr_sample_cols <- grep("^ASD_", colnames(dia_pr_raw), value = TRUE)

# Collapse charge states: median per Protein.Ids + Stripped.Sequence
dia_pep <- dia_pr_raw %>%
  select(Protein.Ids, Stripped.Sequence, all_of(dia_pr_sample_cols)) %>%
  group_by(Protein.Ids, Stripped.Sequence) %>%
  summarise(across(all_of(dia_pr_sample_cols), ~median(.x, na.rm = TRUE)), .groups = "drop")

# Shared peptide filtering: >= 70% detected in both Control and Rest
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

# Log2 -> quantile norm -> z-score
pep_log2 <- log2(pep_mat)

pep_qnorm <- normalize.quantiles(pep_log2)
colnames(pep_qnorm) <- dia_pr_sample_cols

pep_z <- t(apply(pep_qnorm, 1, function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA, length(x)))
  (x - m) / s
}))
colnames(pep_z) <- dia_pr_sample_cols

dia_pep_processed <- dia_pep %>%
  select(Protein.Ids, Stripped.Sequence) %>%
  bind_cols(as_tibble(pep_z))

cat(sprintf("  DIA peptide: %d peptides x %d samples (shared + quantile norm + z-score)\n",
            nrow(dia_pep_processed), length(dia_pr_sample_cols)))
# --- 8b. Build SNP-peptide pairs: top 3 SNPs per gene, all peptides ---
pqtl_suppl <- pqtl_final %>%
  filter(GeneSymbol %in% risk_genes) %>%
  group_by(GeneSymbol) %>%
  mutate(snp_rank = dense_rank(p_cond)) %>%
  filter(snp_rank <= 3) %>%
  ungroup() %>%
  mutate(pep_sequence = str_extract(peptide_id, "(?<=_)[A-Z]+$"))

# Check DIA availability
pqtl_suppl <- pqtl_suppl %>%
  rowwise() %>%
  mutate(in_dia = nrow(dia_pep_processed %>%
    filter(Protein.Ids == risk_uniprot[GeneSymbol],
           Stripped.Sequence == pep_sequence)) > 0) %>%
  ungroup()

cat(sprintf("  Suppl SNP-peptide pairs: %d total, %d in DIA\n",
            nrow(pqtl_suppl), sum(pqtl_suppl$in_dia)))

# --- 8c. Suppl S1: Peptide-level pQTL replication ---
cat("\n--- Suppl S1: Peptide-level pQTL replication ---\n")

geno_mapping_pep <- geno_mapping_dia %>%
  filter(dia_id %in% dia_pr_sample_cols)

spqtl_rep   <- list()
spqtl_plots <- list()

for (i in seq_len(nrow(pqtl_suppl))) {
  gene       <- pqtl_suppl$GeneSymbol[i]
  snp_id     <- pqtl_suppl$snp_id[i]
  snp_region <- pqtl_suppl$region[i]
  uniprot    <- risk_uniprot[gene]
  pep_seq    <- pqtl_suppl$pep_sequence[i]
  tmt_beta   <- pqtl_suppl$beta_cond[i]
  in_dia     <- pqtl_suppl$in_dia[i]

  if (!in_dia) {
    cat(sprintf("  %s | %s | %s: peptide not in DIA, skip\n", gene, snp_id, pep_seq)); next
  }

  geno_row <- geno_case %>% filter(variantid == snp_id)
  if (nrow(geno_row) == 0) {
    cat(sprintf("  %s | %s: SNP not in genotype, skip\n", gene, snp_id)); next
  }

  dia_pep_row <- dia_pep_processed %>%
    filter(Protein.Ids == uniprot, Stripped.Sequence == pep_seq)

  df <- map_dfr(seq_len(nrow(geno_mapping_pep)), function(j) {
    tibble(sample    = geno_mapping_pep$dia_id[j],
           genotype  = as.integer(geno_row[[geno_mapping_pep$geno_id[j]]]),
           abundance = as.numeric(dia_pep_row[[geno_mapping_pep$dia_id[j]]]))
  }) %>%
    filter(!is.na(genotype), !is.na(abundance)) %>%
    mutate(genotype_f = factor(genotype))

  if (nrow(df) < 5 | length(unique(df$genotype)) < 2) {
    cat(sprintf("  %s | %s | %s: skipped (N=%d, geno=%d)\n",
                gene, snp_id, pep_seq, nrow(df), length(unique(df$genotype)))); next
  }

  fit      <- lm(abundance ~ genotype, data = df)
  dia_beta <- coef(fit)["genotype"]
  dia_p    <- summary(fit)$coefficients["genotype","Pr(>|t|)"]

  pair_id <- sprintf("%s_%s_%s", gene, snp_id, pep_seq)

  spqtl_rep[[pair_id]] <- tibble(
    gene = gene, snp = snp_id, region = snp_region,
    peptide = pep_seq, tmt_beta = tmt_beta,
    dia_beta = dia_beta, dia_p = dia_p, n = nrow(df),
    concordant = sign(tmt_beta) == sign(dia_beta))

  spqtl_plots[[pair_id]] <- df %>%
    mutate(gene_label = sprintf("%s\n%s\n(%s)", gene, pep_seq, snp_region))

  cat(sprintf("  %s [%s] %s: TMT b=%.3f, DIA b=%.3f, p=%.2e, concordant=%s (N=%d)\n",
              gene, snp_region, pep_seq, tmt_beta, dia_beta, dia_p,
              sign(tmt_beta) == sign(dia_beta), nrow(df)))
}

spqtl_rep_df  <- bind_rows(spqtl_rep)
spqtl_plot_df <- bind_rows(spqtl_plots)

if (nrow(spqtl_plot_df) > 0) {
  splot_order <- unique(spqtl_plot_df$gene_label)
  spqtl_plot_df$gene_label <- factor(spqtl_plot_df$gene_label, levels = splot_order)

  sstat_labels <- spqtl_rep_df %>%
    mutate(gene_label = factor(sprintf("%s\n%s\n(%s)", gene, peptide, region),
                               levels = splot_order),
           label = sprintf("b=%.2f\np=%.1e", dia_beta, dia_p))

  n_facets <- length(splot_order)
  n_row <- ceiling(n_facets / 6)

  suppl_s1 <- ggplot(spqtl_plot_df, aes(x = genotype_f, y = abundance)) +
    geom_boxplot(aes(fill = genotype_f), width = 0.55, outlier.shape = NA,
                 linewidth = 0.3, alpha = 0.65) +
    geom_jitter(width = 0.12, size = 0.9, alpha = 0.7, shape = 16) +
    geom_text(data = sstat_labels, aes(label = label),
              x = Inf, y = Inf, hjust = 1.05, vjust = 1.2,
              size = 1.4, colour = "grey30", lineheight = 0.85,
              inherit.aes = FALSE) +
    facet_wrap(~gene_label, scales = "free_y", nrow = n_row) +
    scale_fill_manual(values = pal_geno) +
    labs(x = "Genotype dosage", y = "Peptide abundance (z-score)") +
    theme_ng() +
    theme(legend.position = "none", strip.text = element_text(size = 5))

  cat("Suppl S1 done.\n")
} else {
  suppl_s1 <- NULL
  cat("Suppl S1: no valid pairs.\n")
}

# --- 8d. Suppl S2: Peptide-level expression validation ---
cat("\n--- Suppl S2: Peptide-level expression validation ---\n")

suppl_peptides <- pqtl_suppl %>%
  distinct(GeneSymbol, pep_sequence) %>%
  rowwise() %>%
  mutate(in_dia = nrow(dia_pep_processed %>%
    filter(Protein.Ids == risk_uniprot[GeneSymbol],
           Stripped.Sequence == pep_sequence)) > 0) %>%
  ungroup() %>%
  filter(in_dia)

spep_expr <- map_dfr(seq_len(nrow(suppl_peptides)), function(i) {
  gene    <- suppl_peptides$GeneSymbol[i]
  uniprot <- risk_uniprot[gene]
  pep_seq <- suppl_peptides$pep_sequence[i]

  dia_row <- dia_pep_processed %>%
    filter(Protein.Ids == uniprot, Stripped.Sequence == pep_seq)
  if (nrow(dia_row) == 0) return(NULL)

  dia_row %>%
    select(all_of(dia_pr_sample_cols)) %>%
    pivot_longer(everything(), names_to = "Sample_id", values_to = "abundance") %>%
    mutate(Gene = gene, peptide = pep_seq) %>%
    filter(!is.na(abundance))
})

if (nrow(spep_expr) > 0) {
  spep_expr <- spep_expr %>%
    left_join(dia_meta, by = "Sample_id") %>%
    filter(!is.na(Group)) %>%
    mutate(facet_label = sprintf("%s (%s)", Gene, peptide))

  spep_stats <- spep_expr %>%
    group_by(Gene, peptide) %>%
    summarise(
      asd_mean     = mean(abundance[Group == "ASD"], na.rm = TRUE),
      non_asd_mean = mean(abundance[Group != "ASD"], na.rm = TRUE),
      p_wilcox     = tryCatch(wilcox.test(abundance[Group == "ASD"],
                                           abundance[Group != "ASD"])$p.value,
                              error = function(e) NA),
      n_asd  = sum(Group == "ASD" & !is.na(abundance)),
      n_rest = sum(Group != "ASD" & !is.na(abundance)),
      .groups = "drop"
    ) %>%
    mutate(dia_direction = ifelse(asd_mean > non_asd_mean, "Up", "Down")) %>%
    left_join(
      tmt_dep %>% filter(Gene %in% risk_genes) %>%
        transmute(Gene, tmt_fc = ASD_nonASD_fc,
                  tmt_direction = ifelse(ASD_nonASD_fc > 1, "Up", "Down")),
      by = "Gene"
    ) %>%
    mutate(concordant = dia_direction == tmt_direction)

  cat("Peptide expression validation:\n")
  print(spep_stats %>% select(Gene, peptide, dia_direction, p_wilcox, tmt_direction, concordant))

  facet_order <- spep_expr %>%
    distinct(Gene, facet_label) %>%
    mutate(Gene = factor(Gene, levels = risk_genes)) %>%
    arrange(Gene, facet_label) %>%
    pull(facet_label)
  spep_expr$facet_label <- factor(spep_expr$facet_label, levels = facet_order)

  n_facets_s2 <- length(facet_order)
  n_row_s2 <- ceiling(n_facets_s2 / 5)

  suppl_s2 <- ggplot(spep_expr, aes(x = Group, y = abundance, fill = Group)) +
    geom_violin(linewidth = 0.2, alpha = 0.45, scale = "width", width = 0.75, colour = NA) +
    geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.25,
                 alpha = 0.85, colour = "black") +
    facet_wrap(~facet_label, scales = "free_y", nrow = n_row_s2) +
    scale_fill_manual(values = pal_group) +
    labs(x = NULL, y = "Peptide abundance (z-score)") +
    theme_ng() +
    theme(legend.position = "bottom", legend.margin = margin(-2,0,0,0),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text = element_text(size = 5))

  cat("Suppl S2 done.\n")
} else {
  suppl_s2 <- NULL
  cat("Suppl S2: no peptides matched.\n")
}

# --- 8e. Save peptide-level plots ---
cat("\n=== Saving peptide-level plots ===\n")

if (!is.null(suppl_s1)) {
  s1_h <- max(50, n_row * 45)
  ggsave(file.path(output_dir, "suppl_s1_peptide_pqtl.pdf"),
         suppl_s1, width = 180, height = s1_h, units = "mm", device = cairo_pdf)
}

if (!is.null(suppl_s2)) {
  s2_h <- max(90, n_row_s2 * 45)
  ggsave(file.path(output_dir, "suppl_s2_peptide_expression.pdf"),
         suppl_s2, width = 180, height = s2_h, units = "mm", device = cairo_pdf)
}

if (!is.null(suppl_s1) & !is.null(suppl_s2)) {
  suppl_composite <- suppl_s1 / suppl_s2 +
    plot_layout(heights = c(1, 2)) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold",
                                                          family = "Helvetica")))
  ggsave(file.path(output_dir, "DIA_peptide_validation.pdf"),
         suppl_composite, width = 180, height = s1_h + s2_h,
         units = "mm", device = cairo_pdf)
}

if (exists("spqtl_rep_df") && nrow(spqtl_rep_df) > 0)
  write_csv(spqtl_rep_df, file.path(output_dir, "stats_suppl_peptide_pqtl.csv"))
if (exists("spep_stats") && nrow(spep_stats) > 0)
  write_csv(spep_stats, file.path(output_dir, "stats_suppl_peptide_expression.csv"))

cat(sprintf("\nAll plots saved to: %s\n", output_dir))
cat("\n=== COMPLETE ===\n")

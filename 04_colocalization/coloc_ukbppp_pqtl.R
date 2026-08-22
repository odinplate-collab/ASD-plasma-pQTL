################################################################################
## ASD Colocalization — Strategy 2
## UKB-PPP pQTL (N=33,493) x ASD GWAS (PGC 2019, N=46,351)
##
## Targets: MBL2, C1RL, AHSG, C3, CFH, F11 (C6 not in UKB-PPP)
################################################################################

library(tidyverse)
library(coloc)

setwd(Sys.getenv("ASD_ROOT", "."))

# ==============================================================================
# 1. LOAD ASD GWAS (PGC 2019, GRCh37)
# ==============================================================================

cat("=== Loading ASD GWAS ===\n")
gwas_file <- "GWAS/PGC_GWAS_/iPSYCH-PGC_ASD_Nov2017.gz"
gwas <- read_tsv(gwas_file, show_col_types = FALSE)
gwas <- gwas %>%
  mutate(BETA = log(OR)) %>%
  filter(!is.na(BETA), !is.na(SE), SE > 0, !is.na(P))

n_gwas <- 46351
n_cases <- 18382
s_gwas <- n_cases / n_gwas
cat(sprintf("  ASD GWAS: %d SNPs\n", nrow(gwas)))

# ==============================================================================
# 2. DEFINE TARGETS (GRCh37 coordinates)
# ==============================================================================

targets <- tribble(
  ~gene, ~uniprot, ~olink_id, ~chr, ~gene_start, ~gene_end,
  "MBL2",  "P11226", "OID30759", 10, 54531235, 54535093,
  "C1RL",  "Q9NZP8", "OID30721", 12,  7695532,  7741689,
  "AHSG",  "P02765", "OID30706",  3, 186329052, 186340389,
  "C3",    "P01024", "OID30776", 19,  6677704,  6720650,
  "CFH",   "P08603", "OID30790",  1, 196621008, 196716634,
  "F11",   "P03951", "OID30773",  4, 187187686, 187212699
)

window <- 500000
ukb_dir <- "Co_localization/UKB_PPP"

# UKB-PPP folder name mapping
ukb_folders <- c(
  MBL2 = "MBL2/MBL2_P11226_OID30759_v1_Inflammation_II",
  C1RL = "C1RL/C1RL_Q9NZP8_OID30721_v1_Inflammation_II",
  AHSG = "AHSG/AHSG_P02765_OID30706_v1_Inflammation_II",
  C3   = "C3/C3_P01024_OID30776_v1_Inflammation_II",
  CFH  = "CFH/CFH_P08603_OID30790_v1_Inflammation_II",
  F11  = "F11/F11_P03951_OID30773_v1_Inflammation_II"
)

# ==============================================================================
# 3. RUN COLOC
# ==============================================================================

cat("\n=== Running colocalization (UKB-PPP pQTL x ASD GWAS) ===\n")

coloc_results <- list()

for (i in seq_len(nrow(targets))) {
  gene_name    <- targets$gene[i]
  gene_chr     <- targets$chr[i]
  region_start <- targets$gene_start[i] - window
  region_end   <- targets$gene_end[i] + window

  cat(sprintf("\n--- %s (chr%d:%d-%d) ---\n", gene_name, gene_chr, region_start, region_end))

  # --- Load UKB-PPP pQTL for this chromosome ---
  chr_file <- file.path(ukb_dir, ukb_folders[gene_name],
                        sprintf("discovery_chr%d_%s.gz",
                                gene_chr,
                                str_replace(basename(ukb_folders[gene_name]), "/$", "")))

  if (!file.exists(chr_file)) {
    # Try listing files to find correct name
    folder_path <- file.path(ukb_dir, ukb_folders[gene_name])
    all_files <- list.files(folder_path, pattern = sprintf("chr%d_", gene_chr), full.names = TRUE)
    if (length(all_files) == 0) {
      cat(sprintf("  File not found for %s chr%d, skipping\n", gene_name, gene_chr)); next
    }
    chr_file <- all_files[1]
  }

  cat(sprintf("  Loading: %s\n", basename(chr_file)))
  ukb_pqtl <- read.table(gzfile(chr_file), header = TRUE) %>% as_tibble()

  # UKB-PPP stores GENPOS in GRCh38, whereas the ID field encodes the GRCh37
  # position (chr:pos:allele0:allele1). The ASD GWAS is GRCh37, so the GRCh37
  # position is parsed from ID and used for region filtering and matching.
  ukb_locus <- ukb_pqtl %>%
    mutate(POS37 = suppressWarnings(as.integer(sapply(strsplit(as.character(ID), ":"), `[`, 2)))) %>%
    filter(!is.na(POS37), POS37 >= region_start, POS37 <= region_end,
           !is.na(BETA), !is.na(SE), SE > 0, !is.na(LOG10P)) %>%
    mutate(P = 10^(-LOG10P)) %>%
    distinct(POS37, .keep_all = TRUE)

  cat(sprintf("  UKB-PPP pQTL: %d SNPs in region\n", nrow(ukb_locus)))

  if (nrow(ukb_locus) == 0) {
    cat(sprintf("  No UKB-PPP SNPs for %s, skipping\n", gene_name)); next
  }

  # --- GWAS locus ---
  gwas_locus <- gwas %>%
    filter(CHR == gene_chr, BP >= region_start, BP <= region_end) %>%
    filter(!is.na(BETA), !is.na(SE), SE > 0) %>%
    distinct(BP, .keep_all = TRUE)

  cat(sprintf("  ASD GWAS: %d SNPs in region\n", nrow(gwas_locus)))

  if (nrow(gwas_locus) == 0) {
    cat(sprintf("  No GWAS SNPs for %s, skipping\n", gene_name)); next
  }

  # --- Match by GRCh37 position ---
  merged <- inner_join(
    ukb_locus %>% transmute(pos = POS37, pqtl_beta = BETA, pqtl_se = SE, pqtl_p = P,
                             pqtl_snp = paste0(CHROM, ":", POS37)),
    gwas_locus %>% transmute(pos = BP, gwas_beta = BETA, gwas_se = SE, gwas_p = P),
    by = "pos"
  ) %>% distinct(pos, .keep_all = TRUE)

  cat(sprintf("  Matched SNPs: %d\n", nrow(merged)))

  if (nrow(merged) < 50) {
    cat(sprintf("  Too few matched SNPs (<50) for %s, skipping\n", gene_name)); next
  }

  # --- coloc.abf ---
  n_ukb <- ukb_locus$N[1]  # sample size from data

  result <- tryCatch({
    coloc.abf(
      dataset1 = list(
        beta    = merged$pqtl_beta,
        varbeta = merged$pqtl_se^2,
        N       = n_ukb,
        sdY     = 1,
        type    = "quant",
        snp     = merged$pqtl_snp
      ),
      dataset2 = list(
        beta    = merged$gwas_beta,
        varbeta = merged$gwas_se^2,
        N       = n_gwas,
        s       = s_gwas,
        type    = "cc",
        snp     = merged$pqtl_snp
      )
    )
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message)); return(NULL)
  })

  if (is.null(result)) next

  pp <- result$summary

  # Top SNP
  top_snp_name <- NA_character_; top_snp_pp4 <- NA_real_
  if (!is.null(result$results)) {
    res_df <- as.data.frame(result$results)
    top_idx <- which.max(res_df$SNP.PP.H4)
    top_snp_name <- res_df$snp[top_idx]
    top_snp_pp4  <- res_df$SNP.PP.H4[top_idx]
  }

  min_pqtl_p <- min(merged$pqtl_p, na.rm = TRUE)
  min_gwas_p <- min(merged$gwas_p, na.rm = TRUE)

  cat(sprintf("  PP.H0=%.3f  PP.H1=%.3f  PP.H2=%.3f  PP.H3=%.3f  PP.H4=%.3f\n",
              pp["PP.H0.abf"], pp["PP.H1.abf"], pp["PP.H2.abf"],
              pp["PP.H3.abf"], pp["PP.H4.abf"]))
  cat(sprintf("  min pQTL p=%.2e | min GWAS p=%.2e | N_pQTL=%d\n",
              min_pqtl_p, min_gwas_p, n_ukb))
  cat(sprintf("  Top SNP: %s (PP.H4=%.3f)\n", top_snp_name, top_snp_pp4))

  coloc_results[[gene_name]] <- tibble(
    strategy    = "UKB-PPP x ASD",
    gene        = gene_name,
    n_pqtl      = n_ukb,
    n_snps      = nrow(merged),
    min_pqtl_p  = min_pqtl_p,
    min_gwas_p  = min_gwas_p,
    PP.H0       = pp["PP.H0.abf"],
    PP.H1       = pp["PP.H1.abf"],
    PP.H2       = pp["PP.H2.abf"],
    PP.H3       = pp["PP.H3.abf"],
    PP.H4       = pp["PP.H4.abf"],
    top_snp     = top_snp_name,
    top_snp_pp4 = top_snp_pp4
  )
}

# ==============================================================================
# 4. SUMMARY
# ==============================================================================

coloc_df <- bind_rows(coloc_results)

cat("\n\n========================================\n")
cat("=== STRATEGY 2: UKB-PPP x ASD GWAS ===\n")
cat("========================================\n\n")

coloc_df %>%
  mutate(evidence = case_when(
    PP.H4 > 0.8  ~ "*** STRONG ***",
    PP.H4 > 0.5  ~ "** Moderate **",
    PP.H4 > 0.2  ~ "* Suggestive *",
    PP.H3 > 0.5  ~ "Distinct variants",
    PP.H1 > 0.5  ~ "pQTL only",
    TRUE          ~ "Inconclusive"
  )) %>%
  select(gene, n_pqtl, n_snps, min_pqtl_p, min_gwas_p, PP.H3, PP.H4, evidence) %>%
  arrange(desc(PP.H4)) %>%
  print()

# Compare with Strategy 1
cat("\n=== Comparison: Strategy 1 (our pQTL) vs Strategy 2 (UKB-PPP) ===\n")
coloc_s1 <- read_csv("coloc_results_summary.csv", show_col_types = FALSE) %>%
  select(gene, PP.H4_s1 = PP.H4, PP.H3_s1 = PP.H3)

comparison <- coloc_df %>%
  select(gene, PP.H4_s2 = PP.H4, PP.H3_s2 = PP.H3) %>%
  left_join(coloc_s1, by = "gene") %>%
  mutate(PP.H4_improvement = PP.H4_s2 - PP.H4_s1)

print(comparison)

# Save
write_csv(coloc_df, "coloc_strategy2_UKB_PPP_results.csv")
cat("\nResults saved.\n")


################################################################################
## Colocalization summary plot
## Panel a: PP.H4 heatmap (Strategy 1 + 2 + 3)
## Panel b: Min GWAS p vs PP.H4 scatterplot
################################################################################

library(tidyverse)
library(patchwork)

setwd(Sys.getenv("ASD_ROOT", "."))

theme_ng <- function(base_size = 7) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text             = element_text(family = "Helvetica", colour = "black"),
      axis.text        = element_text(size = rel(0.9), colour = "black"),
      axis.title       = element_text(size = rel(1)),
      axis.line        = element_line(linewidth = 0.3),
      axis.ticks       = element_line(linewidth = 0.3),
      panel.grid       = element_blank()
    )
}

# ==============================================================================
# 1. LOAD ALL RESULTS
# ==============================================================================

# Strategy 1: our pQTL x ASD GWAS
s1 <- read_csv("coloc_results_summary.csv", show_col_types = FALSE) %>%
  mutate(strategy = "Our pQTL\n(N=48)",
         trait = "ASD",
         min_gwas_p = NA_real_)  # need to add if available

# Strategy 2: Niu pQTL x ASD GWAS
s2 <- read_csv("coloc_strategy2_Niu_pQTL_results.csv", show_col_types = FALSE) %>%
  mutate(strategy = "Niu pQTL\n(N=1,909)",
         trait = "ASD")

# Strategy 3: our pQTL x trait GWAS
s3 <- read_csv("coloc_strategy3_trait_GWAS_results.csv", show_col_types = FALSE) %>%
  mutate(strategy = "Our pQTL\n(N=48)")

# ==============================================================================
# 2. COMBINE FOR HEATMAP
# ==============================================================================

# Unify columns
heatmap_data <- bind_rows(
  s1 %>% transmute(gene, trait, strategy, PP.H3, PP.H4, min_gwas_p),
  s2 %>% transmute(gene, trait, strategy, PP.H3, PP.H4, min_gwas_p),
  s3 %>% transmute(gene, trait, strategy, PP.H3, PP.H4, min_gwas_p)
)

# Create combined label: strategy + trait for x-axis
heatmap_data <- heatmap_data %>%
  mutate(
    x_label = case_when(
      strategy == "Our pQTL\n(N=48)" & trait == "ASD"  ~ "ASD\n(Our pQTL)",
      strategy == "Niu pQTL\n(N=1,909)" & trait == "ASD" ~ "ASD\n(Niu pQTL)",
      trait == "SCZ"  ~ "SCZ",
      trait == "ADHD" ~ "ADHD",
      trait == "IQ"   ~ "IQ",
      trait == "EA"   ~ "EA"
    ),
    x_label = factor(x_label, levels = c("ASD\n(Our pQTL)", "ASD\n(Niu pQTL)",
                                          "SCZ", "ADHD", "IQ", "EA")),
    gene = factor(gene, levels = c("MBL2", "AHSG", "C3", "C1RL", "CFH", "C6", "F11"))
  )

# ==============================================================================
# 3. PANEL A: PP.H4 HEATMAP
# ==============================================================================

p_a <- ggplot(heatmap_data, aes(x = x_label, y = gene, fill = PP.H4)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", PP.H4)), size = 1.8, colour = "black") +
  scale_fill_gradient2(low = "white", mid = "#FED976", high = "#E31A1C",
                       midpoint = 0.1, limits = c(0, 0.25),
                       oob = scales::squish,
                       name = "PP.H4") +
  # Add vertical separator between Strategy 2 and 3
  geom_vline(xintercept = 2.5, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  labs(x = NULL, y = NULL, title = "Colocalization posterior probability (PP.H4)") +
  # Annotate strategy groups
  annotate("text", x = 1.5, y = 7.7, label = "pQTL × ASD GWAS",
           size = 2, fontface = "bold", colour = "grey30") +
  annotate("text", x = 4.5, y = 7.7, label = "Our pQTL × Trait GWAS",
           size = 2, fontface = "bold", colour = "grey30") +
  coord_cartesian(clip = "off", ylim = c(0.5, 7.5)) +
  theme_ng(base_size = 7) +
  theme(
    axis.text.x = element_text(size = 5.5, lineheight = 0.9),
    plot.title = element_text(size = 7, face = "bold"),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.25, "cm"),
    legend.title = element_text(size = 5.5),
    legend.text = element_text(size = 5),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5)
  )

# ==============================================================================
# 4. PANEL B: GWAS POWER vs PP.H4 SCATTERPLOT
# ==============================================================================

# Combine all data with min_gwas_p
scatter_data <- bind_rows(
  s2 %>% transmute(gene, trait = "ASD (Niu)", PP.H4, PP.H3, min_gwas_p,
                    strategy = "Niu pQTL × ASD"),
  s3 %>% transmute(gene, trait, PP.H4, PP.H3, min_gwas_p,
                    strategy = "Our pQTL × Trait")
) %>%
  filter(!is.na(min_gwas_p)) %>%
  mutate(
    neg_log_p = -log10(min_gwas_p),
    label = paste0(gene, "×", trait),
    # Highlight key results
    highlight = case_when(
      gene == "CFH" & trait == "SCZ" ~ "CFH×SCZ",
      gene == "MBL2" & trait == "SCZ" ~ "MBL2×SCZ",
      gene == "C1RL" & trait == "ASD (Niu)" ~ "C1RL×ASD",
      gene == "MBL2" & trait == "ADHD" ~ "MBL2×ADHD",
      TRUE ~ NA_character_
    ),
    shape_group = ifelse(strategy == "Niu pQTL × ASD", "pQTL × ASD", "pQTL × Trait")
  )

p_b <- ggplot(scatter_data, aes(x = neg_log_p, y = PP.H4)) +
  geom_point(aes(colour = trait, shape = shape_group), size = 1.5, alpha = 0.7) +
  # Label key points
  ggrepel::geom_text_repel(
    data = scatter_data %>% filter(!is.na(highlight)),
    aes(label = highlight),
    size = 2, fontface = "bold",
    max.overlaps = 20,
    box.padding = 0.4,
    segment.size = 0.2,
    min.segment.length = 0
  ) +
  # Reference lines
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "red", linewidth = 0.3, alpha = 0.5) +
  geom_vline(xintercept = -log10(5e-8), linetype = "dashed", colour = "blue",
             linewidth = 0.3, alpha = 0.5) +
  # Annotations
  annotate("text", x = -log10(5e-8) + 0.3, y = 0.19,
           label = "Genome-wide\nsignificance", size = 1.8, colour = "blue", hjust = 0) +
  annotate("text", x = 0.5, y = 0.8,
           label = "Strong coloc\n(PP.H4 > 0.8)", size = 1.8, colour = "red", hjust = 0, vjust = -0.3) +
  scale_colour_manual(
    values = c("ASD (Niu)" = "#E15759", "SCZ" = "#4E79A7", "ADHD" = "#F28E2B",
               "IQ" = "#76B7B2", "EA" = "#59A14F"),
    name = "GWAS Trait"
  ) +
  scale_shape_manual(values = c("pQTL × ASD" = 17, "pQTL × Trait" = 16),
                     name = "Strategy") +
  scale_x_continuous(breaks = seq(0, 8, 2)) +
  labs(x = expression(-log[10]*"(min GWAS "*italic(p)*" in region)"),
       y = "PP.H4",
       title = "GWAS signal strength determines colocalization power") +
  theme_ng(base_size = 7) +
  theme(
    plot.title = element_text(size = 7, face = "bold"),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 5.5),
    legend.text = element_text(size = 5),
    legend.spacing.y = unit(0.1, "cm")
  )

# ==============================================================================
# 5. COMBINE AND SAVE
# ==============================================================================

p_combined <- p_a + p_b +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(size = 9, face = "bold"))
  )

ggsave("Co_localization/coloc_summary.pdf", p_combined,
       width = 180, height = 80, units = "mm", device = cairo_pdf)

cat("Colocalization summary plot saved: coloc_summary.pdf\n")

# ==============================================================================
# 6. SUPPLEMENTARY TABLE (all results combined)
# ==============================================================================

supp_table <- bind_rows(
  s1 %>% transmute(Strategy = "1: Our pQTL × ASD GWAS",
                    pQTL_source = "Our cohort (N=48)", GWAS = "ASD (PGC 2019)",
                    Gene = gene, n_SNPs = n_snps,
                    PP.H0, PP.H1, PP.H2, PP.H3, PP.H4),
  s2 %>% transmute(Strategy = "2: Niu pQTL × ASD GWAS",
                    pQTL_source = "Niu et al. (N=1,909)", GWAS = "ASD (PGC 2019)",
                    Gene = gene, n_SNPs = n_snps,
                    PP.H0, PP.H1, PP.H2, PP.H3, PP.H4),
  s3 %>% transmute(Strategy = "3: Our pQTL × Trait GWAS",
                    pQTL_source = "Our cohort (N=48)", GWAS = paste0(trait, " (N=", trait_N, ")"),
                    Gene = gene, n_SNPs = n_snps,
                    PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)
)

write_csv(supp_table, "coloc_supplementary_table_all_results.csv")
cat("Supplementary table saved: coloc_supplementary_table_all_results.csv\n")

cat("\n=== DONE ===\n")

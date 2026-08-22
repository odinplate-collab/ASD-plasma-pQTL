################################################################################
## ASD pQTL x GWAS Colocalization Analysis
## coloc.abf for cis-pQTL proteins with MR evidence
##
## Targets: MBL2, C1RL, AHSG, C3, C6, CFH, F11
################################################################################

# ==============================================================================
# SETUP
# ==============================================================================

library(tidyverse)
library(coloc)

setwd(Sys.getenv("ASD_ROOT", "."))

# File paths
pqtl_cis_file <- "pQTL/pQTL_re/Results/cis_pQTL_ASD_result.tsv"
gwas_file     <- "GWAS/PGC_GWAS_/iPSYCH-PGC_ASD_Nov2017.gz"

# # alternative GWAS: SPARK+iPSYCH+PGC
# gwas_spark_full <- read_tsv("GWAS/ASD_SPARK_iPSYCH_PGC/ASD_SPARK_iPSYCH_PGC.tsv.gz",
#                              show_col_types = FALSE)

# # harmonize column names
# gwas <- gwas_spark_full %>%
#   rename(CHR = Chromosome, BP = Position, BETA = Effect, SE = StdErr, P = `P-value`) %>%
#   filter(!is.na(BETA), !is.na(SE), SE > 0)

# n_gwas <- 58794
# n_cases <- 19000   # approximate (see source README for the exact number)
# n_controls <- 39794

# then run the coloc loop as below

# Sample sizes
n_pqtl <- 48   # unrelated ASD for pQTL analysis
 n_gwas <- 46351 # 18382 cases + 27969 controls
 n_cases <- 18382
 n_controls <- 27969

# Target proteins (cis MR significant)
# Gene: lead cis-pQTL peptide, approximate gene position (GRCh37)
targets <- tribble(
  ~gene,  ~chr, ~gene_start, ~gene_end,
  "MBL2",   10,  54531235,   54535093,
  "C1RL",   12,   7695532,    7741689,
  "AHSG",    3, 186329052,  186340389,
  "C3",     19,   6677704,    6720650,
  "C6",      5,  41151418,   41171452,
  "CFH",     1, 196621008,  196716634,
  "F11",     4, 187187686,  187212699
)

# Window: ±500kb around gene
window <- 500000

cat("=== Loading data ===\n")

# ==============================================================================
# LOAD pQTL SUMMARY STATS
# ==============================================================================

cat("Loading pQTL cis results...\n")
pqtl_raw <- read_tsv(pqtl_cis_file, show_col_types = FALSE)
cat(sprintf("  pQTL: %d rows\n", nrow(pqtl_raw)))

# Parse SNP column: chr:pos:ref:alt
pqtl <- pqtl_raw %>%
  mutate(
    chr = as.integer(str_extract(SNP, "^\\d+")),
    pos = as.integer(str_extract(SNP, "(?<=:)\\d+(?=:)")),
    SE  = beta / `t-stat`,
    GeneSymbol = str_extract(gene, "(?<=_)[A-Z0-9]+$") # extract gene from peptide_id?
    # Actually gene column is peptide format: P11226_ALQTEMAR
    # We need to map UniProt to gene symbol
  )

# UniProt to gene mapping for targets
uniprot_map <- c(
  P11226 = "MBL2", Q9NZP8 = "C1RL", P02765 = "AHSG", P01024 = "C3",
  P13671 = "C6", P08603 = "CFH", P03951 = "F11"
)

pqtl <- pqtl %>%
  mutate(
    uniprot = str_extract(gene, "^[A-Z0-9]+"),
    GeneSymbol = uniprot_map[uniprot]
  )

cat(sprintf("  pQTL parsed: %d rows with gene mapping\n", sum(!is.na(pqtl$GeneSymbol))))

# ==============================================================================
# LOAD GWAS SUMMARY STATS
# ==============================================================================

cat("Loading GWAS...\n")
gwas <- read_tsv(gwas_file, show_col_types = FALSE)
cat(sprintf("  GWAS: %d SNPs\n", nrow(gwas)))

# Convert OR to beta for GWAS (log(OR))
gwas <- gwas %>%
  mutate(BETA = log(OR))

# # unify column names after loading the GWAS
# gwas <- gwas_spark_full %>%
#   rename(CHR = Chromosome, BP = Position, SNP = MarkerName,
#          BETA = Effect, SE = StdErr, P = `P-value`)
# ==============================================================================
# RUN COLOC FOR EACH TARGET
# ==============================================================================
cat("\n=== Running colocalization ===\n")

coloc_results <- list()

for (i in seq_len(nrow(targets))) {
  gene_name  <- targets$gene[i]
  gene_chr   <- targets$chr[i]
  region_start <- targets$gene_start[i] - window
  region_end   <- targets$gene_end[i] + window

  cat(sprintf("\n--- %s (chr%d:%d-%d) ---\n", gene_name, gene_chr, region_start, region_end))

  # --- Extract pQTL data for this locus ---
  pqtl_locus <- pqtl %>%
    filter(GeneSymbol == gene_name, chr == gene_chr,
           pos >= region_start, pos <= region_end) %>%
    filter(!is.na(beta), !is.na(SE), SE > 0)

  if (nrow(pqtl_locus) == 0) {
    cat(sprintf("  No pQTL SNPs found for %s, skipping\n", gene_name))
    next
  }

  # If multiple peptides, pick the one with lowest p
  best_peptide <- pqtl_locus %>%
    group_by(gene) %>%
    summarise(n = n(), min_p = min(`p-value`), .groups = "drop") %>%
    arrange(min_p) %>%
    slice(1) %>%
    pull(gene)

  pqtl_locus <- pqtl_locus %>%
    filter(gene == best_peptide) %>%
    distinct(pos, .keep_all = TRUE)

  cat(sprintf("  pQTL: %d SNPs (peptide: %s)\n", nrow(pqtl_locus), best_peptide))

  # --- Extract GWAS data for this locus ---
  gwas_locus <- gwas %>%
    filter(CHR == gene_chr, BP >= region_start, BP <= region_end) %>%
    filter(!is.na(BETA), !is.na(SE), SE > 0, !is.na(P))

  cat(sprintf("  GWAS: %d SNPs\n", nrow(gwas_locus)))

  if (nrow(gwas_locus) == 0) {
    cat(sprintf("  No GWAS SNPs for %s, skipping\n", gene_name))
    next
  }

  # --- Match by position ---
  merged <- inner_join(
    pqtl_locus %>% select(pos, pqtl_beta = beta, pqtl_se = SE, pqtl_p = `p-value`, pqtl_snp = SNP),
    gwas_locus %>% select(pos = BP, gwas_beta = BETA, gwas_se = SE, gwas_p = P, gwas_snp = SNP),
    by = "pos"
  )

  cat(sprintf("  Matched SNPs: %d\n", nrow(merged)))

  if (nrow(merged) < 10) {
    cat(sprintf("  Too few matched SNPs (<10) for %s, skipping\n", gene_name))
    next
  }

  merged <- merged %>% distinct(pos, .keep_all = TRUE)

  # --- Run coloc.abf ---
  result <- tryCatch({
    coloc.abf(
      dataset1 = list(
        beta    = merged$pqtl_beta,
        varbeta = merged$pqtl_se^2,
        N       = n_pqtl,
        sdY     = 1,
        type    = "quant",
        snp     = merged$pqtl_snp
      ),
      dataset2 = list(
        beta    = merged$gwas_beta,
        varbeta = merged$gwas_se^2,
        N       = n_gwas,
        s       = n_cases / n_gwas,
        type    = "cc",
        snp     = merged$pqtl_snp
      )
    )
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    return(NULL)
  })

  if (is.null(result)) next

  # Extract posterior probabilities
  pp <- result$summary
  cat(sprintf("  PP.H0=%.3f  PP.H1=%.3f  PP.H2=%.3f  PP.H3=%.3f  PP.H4=%.3f\n",
              pp["PP.H0.abf"], pp["PP.H1.abf"], pp["PP.H2.abf"],
              pp["PP.H3.abf"], pp["PP.H4.abf"]))

  # Top colocalized SNP
  top_snp_name <- NA_character_
  top_snp_pp4  <- NA_real_
  if (!is.null(result$results)) {
    res_df <- as.data.frame(result$results)
    top_idx <- which.max(res_df$SNP.PP.H4)
    top_snp_name <- res_df$snp[top_idx]
    top_snp_pp4  <- res_df$SNP.PP.H4[top_idx]
    cat(sprintf("  Top SNP: %s (PP.H4=%.3f)\n", top_snp_name, top_snp_pp4))
  }

  coloc_results[[gene_name]] <- tibble(
    gene        = gene_name,
    peptide     = best_peptide,
    chr         = gene_chr,
    n_snps      = nrow(merged),
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
# SUMMARY
# ==============================================================================

coloc_df <- bind_rows(coloc_results)

cat("\n=== COLOCALIZATION SUMMARY ===\n")
print(coloc_df %>% select(gene, n_snps, PP.H3, PP.H4, top_snp))

cat("\nInterpretation:\n")
coloc_df %>%
  mutate(evidence = case_when(
    PP.H4 > 0.8  ~ "Strong colocalization",
    PP.H4 > 0.5  ~ "Moderate colocalization",
    PP.H3 > 0.5  ~ "Distinct causal variants",
    PP.H1 > 0.5  ~ "pQTL only",
    PP.H2 > 0.5  ~ "GWAS only",
    TRUE          ~ "Inconclusive"
  )) %>%
  select(gene, PP.H4, evidence) %>%
  print()

# Save results
write_csv(coloc_df, "coloc_results_summary.csv")

# ==============================================================================
# VISUALIZATION — Regional coloc plots
# ==============================================================================
library(patchwork)
cat("\n=== Generating coloc plots ===\n")

theme_ng <- function(base_size = 7) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text             = element_text(family = "Helvetica", colour = "black"),
      axis.text        = element_text(size = rel(0.9), colour = "black"),
      axis.title       = element_text(size = rel(1)),
      axis.line        = element_line(linewidth = 0.3),
      axis.ticks       = element_line(linewidth = 0.3),
      strip.background = element_blank(),
      strip.text       = element_text(size = rel(1), face = "bold"),
      panel.grid       = element_blank()
    )
}

plot_list <- list()

for (gene_name in names(coloc_results)) {
  gene_chr   <- targets %>% filter(gene == gene_name) %>% pull(chr)
  region_start <- targets %>% filter(gene == gene_name) %>% pull(gene_start) - window
  region_end   <- targets %>% filter(gene == gene_name) %>% pull(gene_end) + window
  best_peptide <- coloc_results[[gene_name]]$peptide

  # Get matched data
  pqtl_locus <- pqtl %>%
    filter(GeneSymbol == gene_name, gene == best_peptide,
           chr == gene_chr, pos >= region_start, pos <= region_end) %>%
    distinct(pos, .keep_all = TRUE)

  gwas_locus <- gwas %>%
    filter(CHR == gene_chr, BP >= region_start, BP <= region_end)

  # pQTL plot
  p1 <- ggplot(pqtl_locus, aes(x = pos / 1e6, y = -log10(`p-value`))) +
    geom_point(size = 0.5, alpha = 0.6, colour = "#C44E52") +
    labs(x = NULL, y = expression(-log[10]*"(p) pQTL"),
         title = sprintf("%s (chr%d)", gene_name, gene_chr)) +
    theme_ng()

  # GWAS plot
  p2 <- ggplot(gwas_locus, aes(x = BP / 1e6, y = -log10(P))) +
    geom_point(size = 0.5, alpha = 0.6, colour = "#4C72B0") +
    labs(x = sprintf("Position (Mb, chr%d)", gene_chr),
         y = expression(-log[10]*"(p) GWAS")) +
    theme_ng()

  # PP.H4 annotation
  pp4 <- coloc_results[[gene_name]]$PP.H4
  p1 <- p1 + annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
                        label = sprintf("PP.H4 = %.3f", pp4),
                        size = 2.5, fontface = "bold")

  combined <- p1 / p2
  plot_list[[gene_name]] <- combined
}

# Save individual plots
for (gene_name in names(plot_list)) {
  ggsave(sprintf("coloc_%s.pdf", gene_name),
         plot_list[[gene_name]],
         width = 88, height = 80, units = "mm", device = cairo_pdf)
}

# Combined summary plot
if (length(plot_list) > 0) {
  library(patchwork)

  all_plots <- wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      tag_levels = "a",
      theme = theme(plot.tag = element_text(size = 8, face = "bold", family = "Helvetica"))
    )

  ggsave("Co_localization/coloc_all_loci.pdf", all_plots,
         width = 180, height = ceiling(length(plot_list) / 2) * 80,
         units = "mm", device = cairo_pdf)
}

cat("\n=== COMPLETE ===\n")
cat("Results saved to: coloc_results_summary.csv\n")
cat("Plots saved to: coloc_*.pdf\n")


################################################################################
## ASD pQTL x Neurodevelopmental Trait GWAS Colocalization
## Strategy 3: Our pQTL (N=48) x large-scale trait GWAS
##
## Traits: SCZ (N=306k), ADHD (N=225k), IQ (N=269k), EA (N=766k)
## Targets: MBL2, C1RL, AHSG, C3, C6, CFH, F11
################################################################################

library(tidyverse)
library(coloc)

setwd(Sys.getenv("ASD_ROOT", "."))

# ==============================================================================
# 1. LOAD pQTL (same as before)
# ==============================================================================

cat("=== Loading pQTL ===\n")
pqtl_cis_file <- "pQTL/pQTL_re/Results/cis_pQTL_ASD_result.tsv"
pqtl_raw <- read_tsv(pqtl_cis_file, show_col_types = FALSE)

pqtl <- pqtl_raw %>%
  mutate(
    chr = as.integer(str_extract(SNP, "^\\d+")),
    pos = as.integer(str_extract(SNP, "(?<=:)\\d+(?=:)")),
    SE  = beta / `t-stat`,
    uniprot = str_extract(gene, "^[A-Z0-9]+")
  )

uniprot_map <- c(P11226 = "MBL2", Q9NZP8 = "C1RL", P02765 = "AHSG", P01024 = "C3",
                 P13671 = "C6", P08603 = "CFH", P03951 = "F11")
pqtl$GeneSymbol <- uniprot_map[pqtl$uniprot]

n_pqtl <- 48

# Targets
targets <- tribble(
  ~gene,  ~chr, ~gene_start, ~gene_end,
  "MBL2",   10,  54531235,   54535093,
  "C1RL",   12,   7695532,    7741689,
  "AHSG",    3, 186329052,  186340389,
  "C3",     19,   6677704,    6720650,
  "C6",      5,  41151418,   41171452,
  "CFH",     1, 196621008,  196716634,
  "F11",     4, 187187686,  187212699
)
window <- 500000

cat(sprintf("  pQTL: %d rows, %d with gene mapping\n", nrow(pqtl), sum(!is.na(pqtl$GeneSymbol))))

# ==============================================================================
# 2. LOAD 4 TRAIT GWAS
# ==============================================================================

cat("\n=== Loading GWAS datasets ===\n")

# --- SCZ ---
cat("Loading SCZ...\n")
scz_raw <- read_tsv("GWAS/SCZ/SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz",
                     show_col_types = FALSE, comment = "##")
scz <- scz_raw %>%
  transmute(CHR = CHROM, BP = POS, BETA = BETA, SE = SE, P = PVAL) %>%
  filter(!is.na(BETA), !is.na(SE), SE > 0, !is.na(P))
n_scz <- 306011  # 69369 cases + 236642 controls
s_scz <- 69369 / n_scz
cat(sprintf("  SCZ: %d SNPs (N=%d)\n", nrow(scz), n_scz))

# --- ADHD ---
cat("Loading ADHD...\n")
adhd_raw <- read.table("GWAS/ADHD/ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta",
                        header = TRUE)
adhd <- adhd_raw %>%
  as_tibble() %>%
  transmute(CHR = CHR, BP = BP, BETA = log(OR), SE = SE, P = P) %>%
  filter(!is.na(BETA), !is.na(SE), SE > 0, !is.na(P), is.finite(BETA))
n_adhd <- 225534  # 38691 cases + 186843 controls
s_adhd <- 38691 / n_adhd
cat(sprintf("  ADHD: %d SNPs (N=%d)\n", nrow(adhd), n_adhd))

# --- IQ ---
cat("Loading IQ...\n")
iq_raw <- read_tsv("GWAS/IQ/Savage_2018/SavageJansen_2018_intelligence_metaanalysis.txt",
                    show_col_types = FALSE)
iq <- iq_raw %>%
  transmute(CHR = CHR, BP = POS, BETA = stdBeta, SE = SE, P = P) %>%
  filter(!is.na(BETA), !is.na(SE), SE > 0, !is.na(P))
n_iq <- 269867
cat(sprintf("  IQ: %d SNPs (N=%d)\n", nrow(iq), n_iq))

# --- EA ---
cat("Loading EA...\n")
ea_raw <- read_tsv("GWAS/education/Okbay_27225129-EduYears_Main/EduYears_Main.txt",
                    show_col_types = FALSE)
ea <- ea_raw %>%
  transmute(CHR = CHR, BP = POS, BETA = Beta, SE = SE, P = Pval) %>%
  filter(!is.na(BETA), !is.na(SE), SE > 0, !is.na(P))
n_ea <- 766345
cat(sprintf("  EA: %d SNPs (N=%d)\n", nrow(ea), n_ea))

# ==============================================================================
# 3. DEFINE GWAS LIST
# ==============================================================================

gwas_list <- list(
  SCZ  = list(data = scz,  N = n_scz,  type = "cc",    s = s_scz,  label = "Schizophrenia"),
  ADHD = list(data = adhd, N = n_adhd, type = "cc",    s = s_adhd, label = "ADHD"),
  IQ   = list(data = iq,   N = n_iq,   type = "quant", s = NA,     label = "Intelligence"),
  EA   = list(data = ea,   N = n_ea,   type = "quant", s = NA,     label = "Educational Attainment")
)

# ==============================================================================
# 4. RUN COLOC: 4 traits x 7 proteins = 28 analyses
# ==============================================================================

cat("\n=== Running colocalization (4 traits x 7 proteins) ===\n")

all_results <- list()

for (trait_name in names(gwas_list)) {
  gwas_info <- gwas_list[[trait_name]]
  gwas_data <- gwas_info$data

  cat(sprintf("\n========== %s (%s, N=%d) ==========\n",
              trait_name, gwas_info$label, gwas_info$N))

  for (i in seq_len(nrow(targets))) {
    gene_name    <- targets$gene[i]
    gene_chr     <- targets$chr[i]
    region_start <- targets$gene_start[i] - window
    region_end   <- targets$gene_end[i] + window

    # --- pQTL locus ---
    pqtl_locus <- pqtl %>%
      filter(GeneSymbol == gene_name, chr == gene_chr,
             pos >= region_start, pos <= region_end) %>%
      filter(!is.na(beta), !is.na(SE), SE > 0)

    if (nrow(pqtl_locus) == 0) {
      cat(sprintf("  %s: no pQTL SNPs, skip\n", gene_name)); next
    }

    # Best peptide
    best_peptide <- pqtl_locus %>%
      group_by(gene) %>%
      summarise(n = n(), min_p = min(`p-value`), .groups = "drop") %>%
      arrange(min_p) %>% slice(1) %>% pull(gene)

    pqtl_locus <- pqtl_locus %>%
      filter(gene == best_peptide) %>%
      distinct(pos, .keep_all = TRUE)

    # --- GWAS locus ---
    gwas_locus <- gwas_data %>%
      filter(CHR == gene_chr, BP >= region_start, BP <= region_end) %>%
      filter(!is.na(BETA), !is.na(SE), SE > 0)

    if (nrow(gwas_locus) == 0) {
      cat(sprintf("  %s: no GWAS SNPs, skip\n", gene_name)); next
    }

    # --- Match by position ---
    merged <- inner_join(
      pqtl_locus %>% select(pos, pqtl_beta = beta, pqtl_se = SE, pqtl_p = `p-value`, pqtl_snp = SNP),
      gwas_locus %>% select(pos = BP, gwas_beta = BETA, gwas_se = SE, gwas_p = P),
      by = "pos"
    ) %>% distinct(pos, .keep_all = TRUE)

    if (nrow(merged) < 10) {
      cat(sprintf("  %s: only %d matched SNPs, skip\n", gene_name, nrow(merged))); next
    }

    # --- coloc.abf ---
    # Dataset 1: pQTL (quantitative)
    d1 <- list(
      beta    = merged$pqtl_beta,
      varbeta = merged$pqtl_se^2,
      N       = n_pqtl,
      sdY     = 1,
      type    = "quant",
      snp     = merged$pqtl_snp
    )

    # Dataset 2: trait GWAS
    if (gwas_info$type == "cc") {
      d2 <- list(
        beta    = merged$gwas_beta,
        varbeta = merged$gwas_se^2,
        N       = gwas_info$N,
        s       = gwas_info$s,
        type    = "cc",
        snp     = merged$pqtl_snp
      )
    } else {
      d2 <- list(
        beta    = merged$gwas_beta,
        varbeta = merged$gwas_se^2,
        N       = gwas_info$N,
        sdY     = 1,
        type    = "quant",
        snp     = merged$pqtl_snp
      )
    }

    result <- tryCatch({
      coloc.abf(dataset1 = d1, dataset2 = d2)
    }, error = function(e) {
      cat(sprintf("  %s: ERROR - %s\n", gene_name, e$message)); return(NULL)
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

    # Min GWAS p in region
    min_gwas_p <- min(merged$gwas_p, na.rm = TRUE)

    cat(sprintf("  %s: %d SNPs | PP.H4=%.3f | min_gwas_p=%.2e | top=%s (%.3f)\n",
                gene_name, nrow(merged), pp["PP.H4.abf"], min_gwas_p,
                top_snp_name, top_snp_pp4))

    pair_id <- paste(trait_name, gene_name, sep = "_")
    all_results[[pair_id]] <- tibble(
      trait       = trait_name,
      trait_label = gwas_info$label,
      trait_N     = gwas_info$N,
      gene        = gene_name,
      peptide     = best_peptide,
      n_snps      = nrow(merged),
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
}

# ==============================================================================
# 5. SUMMARY
# ==============================================================================

coloc_all <- bind_rows(all_results)

cat("\n\n========================================\n")
cat("=== FULL COLOCALIZATION SUMMARY ===\n")
cat("========================================\n\n")

coloc_all %>%
  mutate(evidence = case_when(
    PP.H4 > 0.8  ~ "*** STRONG ***",
    PP.H4 > 0.5  ~ "** Moderate **",
    PP.H4 > 0.2  ~ "* Suggestive *",
    PP.H3 > 0.5  ~ "Distinct variants",
    PP.H1 > 0.5  ~ "pQTL only",
    TRUE          ~ "Inconclusive"
  )) %>%
  select(trait, gene, n_snps, min_gwas_p, PP.H3, PP.H4, evidence) %>%
  arrange(desc(PP.H4)) %>%
  print(n = 28)

# Highlight any hits
cat("\n=== PP.H4 > 0.2 (suggestive or better) ===\n")
hits <- coloc_all %>% filter(PP.H4 > 0.2)
if (nrow(hits) > 0) {
  print(hits %>% select(trait, gene, n_snps, PP.H4, top_snp))
} else {
  cat("  None found.\n")
}

# Save
write_csv(coloc_all, "coloc_strategy3_trait_GWAS_results.csv")
cat("\nResults saved to: coloc_strategy3_trait_GWAS_results.csv\n")

# ==============================================================================
# 6. HEATMAP VISUALIZATION
# ==============================================================================

cat("\n=== Generating heatmap ===\n")

library(patchwork)

theme_ng <- function(base_size = 7) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text             = element_text(family = "Helvetica", colour = "black"),
      axis.text        = element_text(size = rel(0.9), colour = "black"),
      axis.title       = element_text(size = rel(1)),
      axis.line        = element_line(linewidth = 0.3),
      axis.ticks       = element_line(linewidth = 0.3),
      strip.background = element_blank(),
      strip.text       = element_text(size = rel(1), face = "bold"),
      panel.grid       = element_blank()
    )
}

# PP.H4 heatmap
heatmap_df <- coloc_all %>%
  mutate(
    gene = factor(gene, levels = c("MBL2","AHSG","C3","C1RL","CFH","C6","F11")),
    trait = factor(trait, levels = c("SCZ","ADHD","IQ","EA"))
  )

p_heat <- ggplot(heatmap_df, aes(x = trait, y = gene, fill = PP.H4)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", PP.H4)), size = 2.2) +
  scale_fill_gradient2(low = "white", mid = "#FED976", high = "#E31A1C",
                       midpoint = 0.4, limits = c(0, 1),
                       name = "PP.H4") +
  labs(x = "Trait GWAS", y = "Protein (cis-pQTL)",
       title = "Colocalization: pQTL x Neurodevelopmental Trait GWAS") +
  theme_ng(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggsave("Co_localization/coloc_strategy3_heatmap.pdf", p_heat,
       width = 100, height = 80, units = "mm", device = cairo_pdf)

cat("Heatmap saved.\n")
cat("\n=== COMPLETE ===\n")


library(tidyverse)
library(patchwork)

# add trait information to strategy-1 results
coloc_asd <- read_csv("Co_localization/coloc_results_summary.csv",
                       show_col_types = FALSE) %>%
  mutate(trait = "ASD", trait_label = "ASD (PGC 2019)", trait_N = 46351)

# strategy-3 results
coloc_trait <- read_csv("Co_localization/coloc_strategy3_trait_GWAS_results.csv",
                         show_col_types = FALSE)

# combine
coloc_all <- bind_rows(
  coloc_asd %>% select(trait, gene, n_snps, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4),
  coloc_trait %>% select(trait, gene, n_snps, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)
)

# Heatmap
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

heatmap_df <- coloc_all %>%
  mutate(
    gene = factor(gene, levels = c("MBL2","AHSG","C3","C1RL","CFH","C6","F11")),
    trait = factor(trait, levels = c("ASD","SCZ","ADHD","IQ","EA"))
  )

p_heat <- ggplot(heatmap_df, aes(x = trait, y = gene, fill = PP.H4)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", PP.H4)), size = 2.2, colour = "black") +
  scale_fill_gradient2(low = "white", mid = "#FED976", high = "#E31A1C",
                       midpoint = 0.4, limits = c(0, 1),
                       name = "PP.H4") +
  labs(x = NULL, y = NULL,
       title = "Colocalization: cis-pQTL x GWAS") +
  scale_x_discrete(labels = c("ASD\n(N=46k)", "SCZ\n(N=306k)", "ADHD\n(N=225k)",
                               "IQ\n(N=270k)", "EA\n(N=766k)")) +
  theme_ng(base_size = 8) +
  theme(axis.text.x = element_text(size = 6),
        panel.grid = element_blank())

ggsave("Co_localization/coloc_combined_heatmap.pdf",
       p_heat, width = 100, height = 75, units = "mm", device = cairo_pdf)

cat("Combined heatmap saved.\n")

# PP.H3 heatmap (Supplementary)
p_heat_h3 <- ggplot(heatmap_df, aes(x = trait, y = gene, fill = PP.H3)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", PP.H3)), size = 2.2, colour = "black") +
  scale_fill_gradient2(low = "white", mid = "#A6BDDB", high = "#2171B5",
                       midpoint = 0.4, limits = c(0, 1),
                       name = "PP.H3") +
  labs(x = NULL, y = NULL,
       title = "Distinct causal variants (PP.H3): cis-pQTL x GWAS") +
  scale_x_discrete(labels = c("ASD\n(N=46k)", "SCZ\n(N=306k)", "ADHD\n(N=225k)",
                               "IQ\n(N=270k)", "EA\n(N=766k)")) +
  theme_ng(base_size = 8) +
  theme(axis.text.x = element_text(size = 6),
        panel.grid = element_blank())

# PP.H4 and PP.H3 side by side (composite panel)
p_combined <- p_heat + p_heat_h3 +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(size = 8, face = "bold")))

ggsave("Co_localization/coloc_combined_heatmap_H4_H3.pdf",
       p_combined, width = 200, height = 75, units = "mm", device = cairo_pdf)

cat("PP.H4 + PP.H3 heatmap saved.\n")

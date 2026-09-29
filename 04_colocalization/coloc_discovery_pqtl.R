# Bayesian colocalization (coloc.abf, gene span +/- 500 kb) of this study's cis-pQTLs (ASD, n = 48, peptide level) with ASD GWAS and with SCZ, ADHD, IQ and EA GWAS
# Inputs (relative to ASD_ROOT):
#   pQTL/pQTL_re/Results/cis_pQTL_ASD_result.tsv   (Matrix eQTL cis output, ASD n = 48; 02_pQTL_mapping/pqtl_mapping.R)
#   GWAS/PGC_GWAS_/iPSYCH-PGC_ASD_Nov2017.gz       (ASD, PGC 2019)
#   GWAS/SCZ/SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz   (schizophrenia, PGC3)
#   GWAS/ADHD/ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta          (ADHD)
#   GWAS/IQ/Savage_2018/SavageJansen_2018_intelligence_metaanalysis.txt  (intelligence)
#   GWAS/education/Okbay_27225129-EduYears_Main/EduYears_Main.txt       (educational attainment)
# Outputs (relative to ASD_ROOT):
#   Co_localization/coloc_results_summary.csv             (this study's pQTL x ASD GWAS)
#   Co_localization/coloc_strategy3_trait_GWAS_results.csv (this study's pQTL x SCZ, ADHD, IQ, EA GWAS)
# Usage: Rscript 04_colocalization/coloc_discovery_pqtl.R
#
# Colocalization of Niu et al. pQTL (n = 1,909) with the ASD GWAS is in coloc_niu_pqtl.R.

ROOT <- Sys.getenv("ASD_ROOT", ".")

library(tidyverse)
library(coloc)

out_dir <- file.path(ROOT, "Co_localization")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pqtl_cis_file <- file.path(ROOT, "pQTL/pQTL_re/Results/cis_pQTL_ASD_result.tsv")
gwas_file     <- file.path(ROOT, "GWAS/PGC_GWAS_/iPSYCH-PGC_ASD_Nov2017.gz")

n_pqtl <- 48   # unrelated ASD individuals used for pQTL mapping

# Target proteins and gene positions (GRCh37)
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

# Window: +/- 500 kb around the gene
window <- 500000

# ==============================================================================
# pQTL summary statistics (SE = beta / t)
# ==============================================================================
pqtl_raw <- read_tsv(pqtl_cis_file, show_col_types = FALSE)

pqtl <- pqtl_raw %>%
  mutate(
    chr = as.integer(str_extract(SNP, "^\\d+")),
    pos = as.integer(str_extract(SNP, "(?<=:)\\d+(?=:)")),
    SE  = beta / `t-stat`,
    uniprot = str_extract(gene, "^[A-Z0-9]+")   # feature ids are "<UniProt>_<PEPTIDE>"
  )

uniprot_map <- c(P11226 = "MBL2", Q9NZP8 = "C1RL", P02765 = "AHSG", P01024 = "C3",
                 P13671 = "C6", P08603 = "CFH", P03951 = "F11")
pqtl$GeneSymbol <- uniprot_map[pqtl$uniprot]

# ==============================================================================
# Strategy 1: this study's pQTL x ASD GWAS (PGC 2019)
# ==============================================================================
n_gwas  <- 46351   # 18382 cases + 27969 controls
n_cases <- 18382

gwas <- read_tsv(gwas_file, show_col_types = FALSE)
gwas <- gwas %>%
  mutate(BETA = log(OR))

coloc_results <- list()

for (i in seq_len(nrow(targets))) {
  gene_name    <- targets$gene[i]
  gene_chr     <- targets$chr[i]
  region_start <- targets$gene_start[i] - window
  region_end   <- targets$gene_end[i] + window

  pqtl_locus <- pqtl %>%
    filter(GeneSymbol == gene_name, chr == gene_chr,
           pos >= region_start, pos <= region_end) %>%
    filter(!is.na(beta), !is.na(SE), SE > 0)

  if (nrow(pqtl_locus) == 0) next

  # Peptide with the smallest P in the region
  best_peptide <- pqtl_locus %>%
    group_by(gene) %>%
    summarise(n = n(), min_p = min(`p-value`), .groups = "drop") %>%
    arrange(min_p) %>%
    slice(1) %>%
    pull(gene)

  pqtl_locus <- pqtl_locus %>%
    filter(gene == best_peptide) %>%
    distinct(pos, .keep_all = TRUE)

  gwas_locus <- gwas %>%
    filter(CHR == gene_chr, BP >= region_start, BP <= region_end) %>%
    filter(!is.na(BETA), !is.na(SE), SE > 0, !is.na(P))

  if (nrow(gwas_locus) == 0) next

  # Match by GRCh37 position
  merged <- inner_join(
    pqtl_locus %>% select(pos, pqtl_beta = beta, pqtl_se = SE, pqtl_p = `p-value`, pqtl_snp = SNP),
    gwas_locus %>% select(pos = BP, gwas_beta = BETA, gwas_se = SE, gwas_p = P, gwas_snp = SNP),
    by = "pos"
  )

  if (nrow(merged) < 10) next

  merged <- merged %>% distinct(pos, .keep_all = TRUE)

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
    cat(sprintf("  %s: ERROR %s\n", gene_name, e$message))
    return(NULL)
  })

  if (is.null(result)) next

  pp <- result$summary

  # SNP with the highest per-SNP PP.H4
  top_snp_name <- NA_character_
  top_snp_pp4  <- NA_real_
  if (!is.null(result$results)) {
    res_df <- as.data.frame(result$results)
    top_idx <- which.max(res_df$SNP.PP.H4)
    top_snp_name <- res_df$snp[top_idx]
    top_snp_pp4  <- res_df$SNP.PP.H4[top_idx]
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

coloc_df <- bind_rows(coloc_results)
write_csv(coloc_df, file.path(out_dir, "coloc_results_summary.csv"))
print(coloc_df %>% select(gene, n_snps, PP.H3, PP.H4, top_snp))

# ==============================================================================
# Strategy 3: this study's pQTL x neurodevelopmental trait GWAS (SCZ, ADHD, IQ, EA)
# ==============================================================================

# SCZ (PGC3 wave 3, European, autosomes)
scz_raw <- read_tsv(file.path(ROOT, "GWAS/SCZ/SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"),
                    show_col_types = FALSE, comment = "##")
scz <- scz_raw %>%
  transmute(CHR = CHROM, BP = POS, BETA = BETA, SE = SE, P = PVAL) %>%
  filter(!is.na(BETA), !is.na(SE), SE > 0, !is.na(P))
n_scz <- 306011  # 69369 cases + 236642 controls
s_scz <- 69369 / n_scz

# ADHD (iPSYCH + deCODE + PGC, 2022); beta = log(OR)
adhd_raw <- read.table(file.path(ROOT, "GWAS/ADHD/ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta"),
                       header = TRUE)
adhd <- adhd_raw %>%
  as_tibble() %>%
  transmute(CHR = CHR, BP = BP, BETA = log(OR), SE = SE, P = P) %>%
  filter(!is.na(BETA), !is.na(SE), SE > 0, !is.na(P), is.finite(BETA))
n_adhd <- 225534  # 38691 cases + 186843 controls
s_adhd <- 38691 / n_adhd

# Intelligence (Savage et al. 2018); standardized beta
iq_raw <- read_tsv(file.path(ROOT, "GWAS/IQ/Savage_2018/SavageJansen_2018_intelligence_metaanalysis.txt"),
                   show_col_types = FALSE)
iq <- iq_raw %>%
  transmute(CHR = CHR, BP = POS, BETA = stdBeta, SE = SE, P = P) %>%
  filter(!is.na(BETA), !is.na(SE), SE > 0, !is.na(P))
n_iq <- 269867

# Educational attainment
ea_raw <- read_tsv(file.path(ROOT, "GWAS/education/Okbay_27225129-EduYears_Main/EduYears_Main.txt"),
                   show_col_types = FALSE)
ea <- ea_raw %>%
  transmute(CHR = CHR, BP = POS, BETA = Beta, SE = SE, P = Pval) %>%
  filter(!is.na(BETA), !is.na(SE), SE > 0, !is.na(P))
n_ea <- 766345

gwas_list <- list(
  SCZ  = list(data = scz,  N = n_scz,  type = "cc",    s = s_scz,  label = "Schizophrenia"),
  ADHD = list(data = adhd, N = n_adhd, type = "cc",    s = s_adhd, label = "ADHD"),
  IQ   = list(data = iq,   N = n_iq,   type = "quant", s = NA,     label = "Intelligence"),
  EA   = list(data = ea,   N = n_ea,   type = "quant", s = NA,     label = "Educational Attainment")
)

all_results <- list()

for (trait_name in names(gwas_list)) {
  gwas_info <- gwas_list[[trait_name]]
  gwas_data <- gwas_info$data

  for (i in seq_len(nrow(targets))) {
    gene_name    <- targets$gene[i]
    gene_chr     <- targets$chr[i]
    region_start <- targets$gene_start[i] - window
    region_end   <- targets$gene_end[i] + window

    pqtl_locus <- pqtl %>%
      filter(GeneSymbol == gene_name, chr == gene_chr,
             pos >= region_start, pos <= region_end) %>%
      filter(!is.na(beta), !is.na(SE), SE > 0)

    if (nrow(pqtl_locus) == 0) next

    best_peptide <- pqtl_locus %>%
      group_by(gene) %>%
      summarise(n = n(), min_p = min(`p-value`), .groups = "drop") %>%
      arrange(min_p) %>% slice(1) %>% pull(gene)

    pqtl_locus <- pqtl_locus %>%
      filter(gene == best_peptide) %>%
      distinct(pos, .keep_all = TRUE)

    gwas_locus <- gwas_data %>%
      filter(CHR == gene_chr, BP >= region_start, BP <= region_end) %>%
      filter(!is.na(BETA), !is.na(SE), SE > 0)

    if (nrow(gwas_locus) == 0) next

    merged <- inner_join(
      pqtl_locus %>% select(pos, pqtl_beta = beta, pqtl_se = SE, pqtl_p = `p-value`, pqtl_snp = SNP),
      gwas_locus %>% select(pos = BP, gwas_beta = BETA, gwas_se = SE, gwas_p = P),
      by = "pos"
    ) %>% distinct(pos, .keep_all = TRUE)

    if (nrow(merged) < 10) next

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
      cat(sprintf("  %s x %s: ERROR %s\n", trait_name, gene_name, e$message)); return(NULL)
    })

    if (is.null(result)) next

    pp <- result$summary

    top_snp_name <- NA_character_; top_snp_pp4 <- NA_real_
    if (!is.null(result$results)) {
      res_df <- as.data.frame(result$results)
      top_idx <- which.max(res_df$SNP.PP.H4)
      top_snp_name <- res_df$snp[top_idx]
      top_snp_pp4  <- res_df$SNP.PP.H4[top_idx]
    }

    # Smallest GWAS P among matched SNPs
    min_gwas_p <- min(merged$gwas_p, na.rm = TRUE)

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

coloc_all <- bind_rows(all_results)
write_csv(coloc_all, file.path(out_dir, "coloc_strategy3_trait_GWAS_results.csv"))
print(coloc_all %>% select(trait, gene, n_snps, min_gwas_p, PP.H3, PP.H4) %>% arrange(desc(PP.H4)), n = 28)

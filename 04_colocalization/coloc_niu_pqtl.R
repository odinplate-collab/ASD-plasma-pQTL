# Bayesian colocalization (coloc.abf, gene span +/- 500 kb) of Niu et al. plasma pQTL with the
# ASD GWAS for the seven complement and coagulation cis-pQTL proteins. Same loci, window,
# matching and settings as strategy 1 in coloc_discovery_pqtl.R; only the pQTL source differs.
# Inputs (relative to ASD_ROOT):
#   GWAS/plasmapQTL/<accession>.tsv.gz        Niu et al. summary statistics (GWAS Catalog; GRCh37)
#   GWAS/PGC_GWAS_/iPSYCH-PGC_ASD_Nov2017.gz   ASD GWAS (18,382 cases, 27,969 controls)
# Outputs (relative to ASD_ROOT): Co_localization/coloc_strategy2_Niu_pQTL_results.csv
# Usage: Rscript coloc_niu_pqtl.R
suppressPackageStartupMessages({library(data.table); library(coloc)})

ROOT <- Sys.getenv("ASD_ROOT", ".")
out_dir <- file.path(ROOT, "Co_localization")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
gwas_file <- file.path(ROOT, "GWAS/PGC_GWAS_/iPSYCH-PGC_ASD_Nov2017.gz")

# Target proteins, gene positions (GRCh37) and Niu et al. accessions. The sample size of each
# accession is read from its file (column n).
targets <- data.table(
  gene       = c("MBL2", "CFH", "C3", "F11", "AHSG", "C1RL", "C6"),
  chr        = c(10L, 1L, 19L, 4L, 3L, 12L, 5L),
  gene_start = c(54531235L, 196621008L, 6677704L, 187187686L, 186329052L, 7695532L, 41151418L),
  gene_end   = c(54535093L, 196716634L, 6720650L, 187212699L, 186340389L, 7741689L, 41171452L),
  accession  = c("GCST90453167", "GCST90453149", "GCST90453052", "GCST90453099",
                 "GCST90453520", "GCST90453393", "GCST90453597"))
window <- 500000

n_gwas  <- 46351
n_cases <- 18382
gwas <- fread(gwas_file, select = c("CHR", "SNP", "BP", "OR", "SE", "P"))
gwas[, BETA := log(OR)]

res <- list()
for (i in seq_len(nrow(targets))) {
  tg <- targets[i]
  lo <- tg$gene_start - window; hi <- tg$gene_end + window
  niu <- fread(file.path(ROOT, "GWAS/plasmapQTL", paste0(tg$accession, ".tsv.gz")),
               select = c("chromosome", "base_pair_location", "beta", "standard_error", "p_value", "variant_id", "n"))
  niu <- niu[chromosome == tg$chr & base_pair_location >= lo & base_pair_location <= hi &
             !is.na(beta) & !is.na(standard_error) & standard_error > 0]
  niu <- unique(niu, by = "base_pair_location")
  n_pqtl <- max(niu$n)
  g <- gwas[CHR == tg$chr & BP >= lo & BP <= hi & !is.na(BETA) & !is.na(SE) & SE > 0 & !is.na(P)]
  # match by GRCh37 position
  m <- merge(niu[, .(pos = base_pair_location, pqtl_beta = beta, pqtl_se = standard_error,
                     pqtl_p = p_value, pqtl_snp = variant_id)],
             g[, .(pos = BP, gwas_beta = BETA, gwas_se = SE, gwas_p = P)], by = "pos")
  m <- unique(m, by = "pos")
  if (nrow(m) < 10) next
  r <- coloc.abf(
    dataset1 = list(beta = m$pqtl_beta, varbeta = m$pqtl_se^2, N = n_pqtl, sdY = 1,
                    type = "quant", snp = m$pqtl_snp),
    dataset2 = list(beta = m$gwas_beta, varbeta = m$gwas_se^2, N = n_gwas,
                    s = n_cases / n_gwas, type = "cc", snp = m$pqtl_snp))
  pp <- r$summary
  top <- which.max(r$results$SNP.PP.H4)
  res[[tg$gene]] <- data.table(
    strategy = "Niu_pQTL x ASD", gene = tg$gene, n_pqtl = n_pqtl, n_snps = nrow(m),
    min_pqtl_p = min(m$pqtl_p), min_gwas_p = min(m$gwas_p),
    PP.H0 = pp[["PP.H0.abf"]], PP.H1 = pp[["PP.H1.abf"]], PP.H2 = pp[["PP.H2.abf"]],
    PP.H3 = pp[["PP.H3.abf"]], PP.H4 = pp[["PP.H4.abf"]],
    top_snp = as.character(r$results$snp[top]), top_snp_pp4 = r$results$SNP.PP.H4[top])
}
out <- rbindlist(res)
fwrite(out, file.path(out_dir, "coloc_strategy2_Niu_pQTL_results.csv"))
print(out[, .(gene, n_pqtl, n_snps, PP.H3 = round(PP.H3, 4), PP.H4 = round(PP.H4, 4), top_snp)])

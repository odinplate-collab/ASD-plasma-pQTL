# pQTL calibration: genome-wide Matrix eQTL run in ASD (n = 48) with the settings of the
# main mapping (pqtl_mapping.R: peptide matrix pep_unique_ASD.txt; covariates = peptide
# principal components 1-5, age and sex; linear model; cis window 1 Mb), on observed data
# (PERM_SEED = 0) or with the genotype sample labels permuted (PERM_SEED > 0). Genotypes are
# permuted, so that each individual's peptide abundances stay with their covariates. The run
# records the number of cis and trans associations below each P-value threshold; the
# thresholds include the FDR 5% cut-offs of the main mapping (cis P <= 1.024e-5, trans P <= 4.445e-7).
# The calibration uses the observed run and 100 permuted runs (seeds 1-100); see pqtl_calibration_table.R.
# Inputs (relative to ASD_ROOT): pQTL/pQTL_re/ASD_WGS_imputated_variantloc_rm29903.tsv,
#   pqtl_peptide_locus.txt, ASD_pQTL_imputated_case_rm29903_n_alt.tsv, pep_unique_ASD.txt, pc_score_ASD.txt
# Outputs (relative to ASD_ROOT): calibration/outputs/fdrcurve_seed<SEED>.csv
# Usage: PERM_SEED=0 Rscript pqtl_calibration_permutation.R   (then PERM_SEED=1 ... 100)
suppressMessages({library(MatrixEQTL); library(data.table)})
ROOT <- Sys.getenv("ASD_ROOT", ".")
RE   <- file.path(ROOT, "pQTL/pQTL_re")
OUT  <- file.path(ROOT, "calibration/outputs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
log  <- function(...) { cat(format(Sys.time()), "-", ..., "\n"); flush.console() }
SEED <- as.integer(Sys.getenv("PERM_SEED", "0"))

BREAKS <- c(0, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 4.445e-7, 1e-6,
            1.024e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 1)

snpspos <- as.data.frame(fread(file.path(RE, "ASD_WGS_imputated_variantloc_rm29903.tsv")))[, c("variantid", "chr", "position")]
snpspos$chr <- paste0("chr", snpspos$chr)
genepos <- as.data.frame(fread(file.path(RE, "pqtl_peptide_locus.txt")))[, c("id", "chr", "s1", "s2")]

f_snp <- file.path(RE, "ASD_pQTL_imputated_case_rm29903_n_alt.tsv")
f_exp <- file.path(RE, "pep_unique_ASD.txt")
f_cov <- file.path(RE, "pc_score_ASD.txt")
# the three matrices must list the same individuals in the same order
hdr <- function(f) gsub('"', "", gsub("-", "_", strsplit(readLines(f, n = 1), "\t")[[1]][-1]))
stopifnot(identical(hdr(f_snp), hdr(f_exp)), identical(hdr(f_snp), hdr(f_cov)))

mk <- function(fn) { s <- SlicedData$new(); s$fileDelimiter <- "\t"; s$fileOmitCharacters <- "NA"
                     s$fileSkipRows <- 1; s$fileSkipColumns <- 1; s$fileSliceSize <- 5000; s$LoadFile(fn); s }
log("loading ...")
snps <- mk(f_snp); gene <- mk(f_exp); cvrt <- mk(f_cov)
log(sprintf("dims: snps=%d gene=%d covariates=%d samples=%d", nrow(snps), nrow(gene), nrow(cvrt), ncol(snps)))

if (SEED == 0) log("observed data") else {
  set.seed(SEED)
  log("genotype labels permuted, seed =", SEED)
  snps$ColumnSubsample(sample(ncol(snps)))
}

TMP <- file.path(OUT, "_tmp"); dir.create(TMP, showWarnings = FALSE)
tc <- file.path(TMP, sprintf("c%d.txt", SEED)); tt <- file.path(TMP, sprintf("t%d.txt", SEED))
t0 <- Sys.time()
me <- Matrix_eQTL_main(
  snps = snps, gene = gene, cvrt = cvrt,
  output_file_name     = tt, pvOutputThreshold     = 1e-30,
  output_file_name.cis = tc, pvOutputThreshold.cis = 1e-30,
  useModel = modelLINEAR, errorCovariance = numeric(),
  snpspos = snpspos, genepos = genepos, cisDist = 1e6,
  pvalue.hist = BREAKS, min.pv.by.genesnp = FALSE, noFDRsaveMemory = TRUE, verbose = FALSE)
log(sprintf("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

grab <- function(part, lbl) data.table(
  seed = SEED, permuted = (SEED != 0), class = lbl,
  p_upper = BREAKS[-1], count = as.numeric(part$hist.counts),
  cum = cumsum(as.numeric(part$hist.counts)), total = sum(as.numeric(part$hist.counts)))
res <- rbind(grab(me$cis, "cis"), grab(me$trans, "trans"))
fwrite(res, file.path(OUT, sprintf("fdrcurve_seed%d.csv", SEED)))
suppressWarnings(file.remove(tc, tt))

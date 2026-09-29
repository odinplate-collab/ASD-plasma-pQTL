# Download CSF pQTL summary statistics (Western et al. 2024) for the 12 panel proteins and cut each to the gene +/- 250 kb.
# Inputs (relative to ASD_ROOT):
#   brain/data/csf_pqtl_download_plan.tsv  one row per CSF pQTL study to fetch (gene, GWAS Catalog
#     accession, trait, sample_size, GRCh38 chr/start/end, ensembl). Curated from the GWAS Catalog
#     entries of Western et al. 2024 (17 studies for the 12 proteins, after dropping aptamers that
#     target another gene: ACHE for BCHE, NPEPL1 for ANPEP, MASP1 complexes for MBL2).
#   GWAS Catalog FTP (https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics), GRCh38 files
# Outputs (relative to ASD_ROOT):
#   brain/data/csf_pqtl/CSF_<gene>_<accession>_cis250kb.tsv.gz
#   brain/data/csf_pqtl_summary.tsv  (lead variant per study)
# Usage: Rscript fetch_csf_pqtl.R
# Each full file (~200 MB) is downloaded with curl, checked against the FTP md5, cut to the window, then deleted.
suppressMessages({library(data.table); library(tools)})
ROOT <- Sys.getenv("ASD_ROOT", ".")
D <- file.path(ROOT, "brain", "data"); OUT <- file.path(D, "csf_pqtl")
TMPD <- file.path(OUT, "_tmp")
dir.create(TMPD, showWarnings = FALSE, recursive = TRUE)
say <- function(...) { cat(format(Sys.time()), "-", ..., "\n"); flush.console() }
WIN <- 250000L

PLAN <- fread(file.path(D, "csf_pqtl_download_plan.tsv"))
PLAN[, `:=`(lo = pmax(1, start - WIN), hi = end + WIN)]
say("studies to fetch:", nrow(PLAN), "for", uniqueN(PLAN$gene), "proteins")

base <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
blockdir <- function(acc) {
  n <- as.numeric(sub("GCST", "", acc))
  lo <- floor((n - 1)/1000)*1000 + 1
  sprintf("GCST%d-GCST%d", lo, lo + 999)
}
CURL <- Sys.which("curl"); stopifnot(nzchar(CURL))
curl_get <- function(url, dest, resume = TRUE) {
  a <- c("-sS", "-L", "--fail", "--retry", "8", "--retry-delay", "15",
         "--retry-all-errors", "--connect-timeout", "60", "--max-time", "1800")
  if (resume) a <- c(a, "--continue-at", "-")
  system2(CURL, c(a, "-o", shQuote(dest), shQuote(url)), stdout = FALSE, stderr = FALSE)
}
expected_md5 <- function(acc) {
  f <- file.path(TMPD, paste0(acc, ".md5"))
  curl_get(sprintf("%s/%s/%s/md5sum.txt", base, blockdir(acc), acc), f, resume = FALSE)
  if (!file.exists(f)) return(NA_character_)
  l <- readLines(f, warn = FALSE); unlink(f)
  # md5sum.txt also lists the -meta.yaml file, possibly first, so match the filename field exactly
  parts <- strsplit(trimws(l), "[[:space:]]+")
  want <- paste0(acc, ".tsv.gz")
  h <- Filter(function(x) length(x) >= 2 && x[2] == want, parts)
  if (!length(h)) return(NA_character_)
  h[[1]][1]
}
fetch <- function(acc, dest) {
  want <- expected_md5(acc)
  url <- sprintf("%s/%s/%s/%s.tsv.gz", base, blockdir(acc), acc, acc)
  for (att in 1:3) {
    curl_get(url, dest)
    if (file.exists(dest) && file.size(dest) > 1e6 &&
        (is.na(want) || identical(unname(md5sum(dest)), want))) return(TRUE)
    say("    attempt", att, "failed")
    if (!is.na(want) && file.exists(dest)) unlink(dest)
    if (att < 3) Sys.sleep(30)
  }
  FALSE
}

done <- 0L; failed <- character()
for (i in seq_len(nrow(PLAN))) {
  g <- PLAN$gene[i]; acc <- PLAN$accession[i]
  out <- file.path(OUT, sprintf("CSF_%s_%s_cis250kb.tsv.gz", g, acc))
  if (file.exists(out)) { done <- done + 1L; next }
  tmp <- file.path(TMPD, paste0(acc, ".tsv.gz"))
  say(sprintf("[%2d/%d] %-8s %s", i, nrow(PLAN), g, acc))
  if (!fetch(acc, tmp)) { failed <- c(failed, paste(g, acc)); unlink(tmp); next }
  d <- try(fread(tmp, showProgress = FALSE), silent = TRUE)
  if (inherits(d, "try-error")) { failed <- c(failed, paste(g, "read")); unlink(tmp); next }
  cc <- names(d)[names(d) %in% c("chromosome", "chr", "CHR")][1]
  pp <- names(d)[names(d) %in% c("base_pair_location", "pos", "BP")][1]
  sub <- d[as.character(get(cc)) == as.character(PLAN$chr[i]) &
           get(pp) >= PLAN$lo[i] & get(pp) <= PLAN$hi[i]]
  if (!nrow(sub)) {
    say("         zero rows in window -- recorded as unavailable")
    cat(sprintf("%s (%s)\n", g, acc), file = file.path(OUT, "NOT_AVAILABLE.txt"),
        append = TRUE)
  } else {
    sub[, `:=`(gene = g, accession = acc)]
    fwrite(sub, out, sep = "\t")
    say(sprintf("         kept %d rows  (max -log10P %.2f)", nrow(sub),
                max(sub$neg_log_10_p_value, na.rm = TRUE)))
    done <- done + 1L
  }
  unlink(tmp); rm(d, sub); gc(verbose = FALSE)
}
unlink(TMPD, recursive = TRUE)
say(sprintf("finished | written %d | failed %d", done, length(failed)))
if (length(failed)) { print(failed); writeLines(failed, file.path(OUT, "FAILED.txt")) }

FL <- list.files(OUT, pattern = "cis250kb[.]tsv[.]gz$", full.names = TRUE)
if (length(FL)) {
  S <- rbindlist(lapply(FL, function(p) {
    d <- fread(p, showProgress = FALSE)
    k <- which.max(d$neg_log_10_p_value)
    d[, .(gene = gene[1], accession = accession[1], n_variants = .N,
          max_neglog10p = round(max(neg_log_10_p_value, na.rm = TRUE), 3),
          lead = variant_id[k], lead_rsid = rs_id[k],
          lead_beta = beta[k], lead_pos = base_pair_location[k], n = n[k])]
  }))
  setorder(S, -max_neglog10p)
  fwrite(S, file.path(D, "csf_pqtl_summary.tsv"), sep = "\t")
  print(S)
  say("CSF cis signal at P < 5e-8:", S[max_neglog10p > 7.301, uniqueN(gene)],
      "of", S[, uniqueN(gene)], "proteins")
}
say("done")

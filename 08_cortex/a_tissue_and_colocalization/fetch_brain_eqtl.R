# Fetch brain cis-eQTL summary statistics (eQTL Catalogue, gene-level, brain/cortex datasets) for the 12 panel genes, gene +/- 250 kb.
# Inputs (relative to ASD_ROOT):
#   brain/data/panel12_grch38.tsv  gene, ensembl, chr, start, end, strand (GRCh38) of the 12 panel genes
#   brain/data/eqtl_paths.tsv      eQTL Catalogue dataset table (tabix_ftp_paths.tsv from the
#                                  eQTL-Catalogue-resources repository)
#   eQTL Catalogue tabix-indexed files on the EBI FTP (queried remotely; nothing large is downloaded)
# Outputs (relative to ASD_ROOT):
#   brain/data/brain_eqtl/<dataset>__<gene>.tsv.gz
#   brain/data/brain_eqtl_summary.tsv
# Usage: Rscript fetch_brain_eqtl.R
# eQTL Catalogue is GRCh38 and the plasma pQTL GRCh37; windows are taken in GRCh38 here and the two
# sides are matched on rsID in coloc_plasma_brain.R.
suppressPackageStartupMessages({
  library(data.table); library(Rsamtools); library(GenomicRanges); library(IRanges)
})
ROOT <- Sys.getenv("ASD_ROOT", ".")
D <- file.path(ROOT, "brain", "data"); OUT <- file.path(D, "brain_eqtl")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
say <- function(...) { cat(format(Sys.time()), "-", ..., "\n"); flush.console() }
WIN <- 250000L

G <- fread(file.path(D, "panel12_grch38.tsv"))
PATHS <- fread(file.path(D, "eqtl_paths.tsv"))
# gene-level ("ge") datasets whose tissue label mentions brain or cortex
BR <- PATHS[quant_method == "ge" &
            (grepl("brain", tissue_label, ignore.case = TRUE) |
             grepl("cortex", tissue_label, ignore.case = TRUE))]
BR[, url := sub("^ftp://", "https://", ftp_path)]
say("brain gene-level datasets:", nrow(BR))
print(BR[, .(dataset_id, study_label, tissue_label, sample_size)])

COLS <- c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant",
          "ma_samples", "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2",
          "molecular_trait_object_id", "gene_id", "median_tpm", "rsid")

done <- 0L; empty <- 0L; failed <- character()
for (i in seq_len(nrow(BR))) {
  ds <- BR$dataset_id[i]
  tf <- try(TabixFile(BR$url[i]), silent = TRUE)
  if (inherits(tf, "try-error")) { failed <- c(failed, ds); next }
  for (j in seq_len(nrow(G))) {
    g <- G$gene[j]
    f <- file.path(OUT, sprintf("%s__%s.tsv.gz", ds, g))
    if (file.exists(f)) { done <- done + 1L; next }
    gr <- GRanges(as.character(G$chr[j]),
                  IRanges(max(1L, G$start[j] - WIN), G$end[j] + WIN))
    x <- try(scanTabix(tf, param = gr), silent = TRUE)
    if (inherits(x, "try-error")) { failed <- c(failed, paste(ds, g)); next }
    ln <- x[[1]]
    if (!length(ln)) { empty <- empty + 1L; next }
    d <- fread(text = paste(ln, collapse = "\n"), header = FALSE,
               col.names = COLS, showProgress = FALSE)
    d <- d[molecular_trait_id == G$ensembl[j]]
    if (!nrow(d)) { empty <- empty + 1L; next }
    d[, `:=`(dataset_id = ds, study = BR$study_label[i],
             tissue = BR$tissue_label[i], n = BR$sample_size[i], gene_symbol = g)]
    fwrite(d, f, sep = "\t")
    done <- done + 1L
    say(sprintf("%-10s %-32s %-8s %6d variants", ds, BR$tissue_label[i], g, nrow(d)))
  }
}
say(sprintf("finished | written/present %d | empty %d | failed %d",
            done, empty, length(failed)))
if (length(failed)) writeLines(failed, file.path(OUT, "FAILED.txt"))

FL <- list.files(OUT, pattern = "[.]tsv[.]gz$", full.names = TRUE)
SUM <- rbindlist(lapply(FL, function(p) {
  d <- fread(p, showProgress = FALSE)
  d[, .(dataset_id = dataset_id[1], study = study[1], tissue = tissue[1],
        n = n[1], gene = gene_symbol[1], n_variants = .N,
        min_p = min(pvalue), lead_rsid = rsid[which.min(pvalue)],
        lead_beta = beta[which.min(pvalue)], median_tpm = median_tpm[1])]
}))
fwrite(SUM, file.path(D, "brain_eqtl_summary.tsv"), sep = "\t")
say("summary rows:", nrow(SUM))
print(SUM[, .(datasets = uniqueN(dataset_id), min_p = min(min_p),
              n_sig_5e8 = sum(min_p < 5e-8), n_sig_1e5 = sum(min_p < 1e-5),
              median_tpm = round(median(median_tpm), 2)), by = gene][order(min_p)])
say("done")

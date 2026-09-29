# Colocalization of the plasma protein-level cis-pQTL (this study) with brain cis-eQTL (eQTL Catalogue) for the 12 panel genes.
# Inputs (relative to ASD_ROOT):
#   derived/14_ld_reference_479/1000G_EUR_479regions.bim   (GRCh37 position -> rsID)
#   derived/12_cis_summary_all_proteins/protein_level_cis_combined_n90_dxadjusted_ALL.csv.gz
#     (plasma protein-level cis statistics, n = 90, diagnosis-adjusted, GRCh37)
#   brain/data/brain_eqtl/<dataset>__<gene>.tsv.gz   (fetch_brain_eqtl.R; GRCh38)
# Outputs (relative to ASD_ROOT): brain/data/coloc_plasma_vs_brain_eqtl.tsv
# Usage: Rscript coloc_plasma_brain.R
# One coloc.abf test per gene and brain dataset over the shared variants; matched on rsID.
suppressPackageStartupMessages({library(data.table); library(coloc)})
ROOT <- Sys.getenv("ASD_ROOT", ".")
PK <- file.path(ROOT, "derived")
D <- file.path(ROOT, "brain", "data"); EQ <- file.path(D, "brain_eqtl")
say <- function(...) { cat(format(Sys.time()), "-", ..., "\n"); flush.console() }
P12 <- c("AHSG","ANPEP","BCHE","BTD","C1RL","C3","CLEC3B","IGFBP5","MBL2",
         "POSTN","PTGDS","QSOX1")

## ---- rsID map for the GRCh37 plasma variants ---------------------------------
say("reading 1000G EUR bim for the rsID map")
BIM <- fread(file.path(PK, "14_ld_reference_479/1000G_EUR_479regions.bim"),
             header = FALSE, select = c(1, 2, 4, 5, 6),
             col.names = c("chr", "rsid", "pos", "a1", "a2"))
BIM[, key := paste(chr, pos, sep = ":")]
setkey(BIM, key)
say("bim variants:", nrow(BIM))

## ---- plasma cis statistics ------------------------------------------------------
say("reading the protein-level cis summary")
CIS <- fread(cmd = paste("gzip -dc", shQuote(file.path(PK,
  "12_cis_summary_all_proteins/protein_level_cis_combined_n90_dxadjusted_ALL.csv.gz"))),
  showProgress = FALSE)
CIS <- CIS[protein %in% P12]
CIS[, key := paste(chr, pos, sep = ":")]
CIS <- merge(CIS, BIM[, .(key, rsid)], by = "key")
CIS <- CIS[!duplicated(paste(protein, rsid))]
say("cis rows with an rsID:", nrow(CIS),
    "| proteins:", paste(sort(unique(CIS$protein)), collapse = ","))

FL <- list.files(EQ, pattern = "[.]tsv[.]gz$", full.names = TRUE)
say("brain eQTL files:", length(FL))

res <- list()
for (f in FL) {
  b <- fread(f, showProgress = FALSE)
  g <- b$gene_symbol[1]
  a <- CIS[protein == g]
  if (!nrow(a)) next
  b <- b[!is.na(rsid) & rsid != "" & !duplicated(rsid)]
  m <- merge(a[, .(rsid, b_p = beta, se_p = se, p_p = p, maf_p = maf, n_p = n)],
             b[, .(rsid, b_e = beta, se_e = se, p_e = pvalue, maf_e = maf,
                   n_e = an/2, dataset_id, study, tissue)], by = "rsid")
  m <- m[is.finite(b_p) & is.finite(se_p) & se_p > 0 &
         is.finite(b_e) & is.finite(se_e) & se_e > 0 &
         maf_p > 0 & maf_p < 1 & maf_e > 0 & maf_e < 1]
  # windows with fewer than 50 shared variants are not tested
  if (nrow(m) < 50) next
  # coloc.abf default priors (p1 = p2 = 1e-4, p12 = 1e-5)
  r <- try(coloc.abf(
    dataset1 = list(beta = m$b_p, varbeta = m$se_p^2, snp = m$rsid,
                    type = "quant", N = max(m$n_p), MAF = m$maf_p),
    dataset2 = list(beta = m$b_e, varbeta = m$se_e^2, snp = m$rsid,
                    type = "quant", N = max(m$n_e), MAF = m$maf_e)),
    silent = TRUE)
  if (inherits(r, "try-error")) next
  s <- as.list(r$summary)
  res[[length(res) + 1]] <- data.table(
    gene = g, dataset_id = m$dataset_id[1], study = m$study[1], tissue = m$tissue[1],
    n_shared = nrow(m), min_p_plasma = min(m$p_p), min_p_brain = min(m$p_e),
    PP.H0 = s$PP.H0.abf, PP.H1 = s$PP.H1.abf, PP.H2 = s$PP.H2.abf,
    PP.H3 = s$PP.H3.abf, PP.H4 = s$PP.H4.abf)
}
R <- rbindlist(res)
setorder(R, -PP.H4)
fwrite(R, file.path(D, "coloc_plasma_vs_brain_eqtl.tsv"), sep = "\t")
say("coloc tests:", nrow(R))
say("PP.H4 > 0.5:", R[PP.H4 > 0.5, .N], "| PP.H4 > 0.8:", R[PP.H4 > 0.8, .N])
print(R[1:15, .(gene, study, tissue, n_shared,
                p_plasma = signif(min_p_plasma, 2), p_brain = signif(min_p_brain, 2),
                H3 = round(PP.H3, 3), H4 = round(PP.H4, 3))])
say("")
say("best PP.H4 per gene:")
print(R[, .(datasets = .N, best_H4 = round(max(PP.H4), 3),
            best_tissue = tissue[which.max(PP.H4)],
            min_p_brain = signif(min(min_p_brain), 2)), by = gene][order(-best_H4)])
say("done")

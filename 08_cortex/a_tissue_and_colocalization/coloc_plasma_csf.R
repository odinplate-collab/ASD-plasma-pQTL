# Colocalization of the plasma protein-level cis-pQTL (this study) with CSF pQTL (Western et al. 2024) for the 12 panel proteins.
# Inputs (relative to ASD_ROOT):
#   derived/14_ld_reference_479/1000G_EUR_479regions.bim   (GRCh37 position -> rsID)
#   derived/12_cis_summary_all_proteins/protein_level_cis_combined_n90_dxadjusted_ALL.csv.gz
#     (plasma protein-level cis statistics, n = 90, diagnosis-adjusted, GRCh37)
#   brain/data/csf_pqtl/CSF_<gene>_<accession>_cis250kb.tsv.gz   (fetch_csf_pqtl.R; GRCh38)
# Outputs (relative to ASD_ROOT): brain/data/coloc_plasma_vs_csf.tsv
# Usage: Rscript coloc_plasma_csf.R
# The two sides are on different genome builds and are matched on rsID.
suppressPackageStartupMessages({library(data.table); library(coloc)})
ROOT <- Sys.getenv("ASD_ROOT", ".")
PK <- file.path(ROOT, "derived")
D <- file.path(ROOT, "brain", "data"); CS <- file.path(D, "csf_pqtl")
say <- function(...) { cat(format(Sys.time()), "-", ..., "\n"); flush.console() }
P12 <- c("AHSG","ANPEP","BCHE","BTD","C1RL","C3","CLEC3B","IGFBP5","MBL2",
         "POSTN","PTGDS","QSOX1")

BIM <- fread(file.path(PK, "14_ld_reference_479/1000G_EUR_479regions.bim"),
             header = FALSE, select = c(1, 2, 4), col.names = c("chr", "rsid", "pos"))
BIM[, key := paste(chr, pos, sep = ":")]
CIS <- fread(cmd = paste("gzip -dc", shQuote(file.path(PK,
  "12_cis_summary_all_proteins/protein_level_cis_combined_n90_dxadjusted_ALL.csv.gz"))),
  showProgress = FALSE)[protein %in% P12]
CIS[, key := paste(chr, pos, sep = ":")]
CIS <- merge(CIS, BIM[, .(key, rsid)], by = "key")
CIS <- CIS[!duplicated(paste(protein, rsid))]
say("plasma cis rows with rsID:", nrow(CIS))

FL <- list.files(CS, pattern = "cis250kb[.]tsv[.]gz$", full.names = TRUE)
say("CSF files:", length(FL))

res <- list()
for (f in FL) {
  b <- fread(f, showProgress = FALSE)
  g <- b$gene[1]; acc <- b$accession[1]
  a <- CIS[protein == g]
  if (!nrow(a)) { say("no plasma cis rows for", g); next }
  b <- b[!is.na(rs_id) & rs_id != "" & !duplicated(rs_id)]
  b[, maf := pmin(effect_allele_frequency, 1 - effect_allele_frequency)]
  m <- merge(a[, .(rsid, b_p = beta, se_p = se, maf_p = maf, n_p = n,
                   nlp_p = -log10(p))],
             b[, .(rsid = rs_id, b_c = beta, se_c = standard_error, maf_c = maf,
                   n_c = n, nlp_c = neg_log_10_p_value)], by = "rsid")
  m <- m[is.finite(b_p) & is.finite(se_p) & se_p > 0 &
         is.finite(b_c) & is.finite(se_c) & se_c > 0 &
         maf_p > 0 & maf_p < 0.5 & maf_c > 0 & maf_c < 0.5]
  # windows with fewer than 30 shared variants are not tested
  if (nrow(m) < 30) { say(sprintf("%-8s %s: only %d shared variants", g, acc, nrow(m))); next }
  # coloc.abf default priors (p1 = p2 = 1e-4, p12 = 1e-5)
  r <- try(coloc.abf(
    dataset1 = list(beta = m$b_p, varbeta = m$se_p^2, snp = m$rsid,
                    type = "quant", N = max(m$n_p), MAF = m$maf_p),
    dataset2 = list(beta = m$b_c, varbeta = m$se_c^2, snp = m$rsid,
                    type = "quant", N = max(m$n_c), MAF = m$maf_c)), silent = TRUE)
  if (inherits(r, "try-error")) { say("coloc failed for", g); next }
  s <- as.list(r$summary)
  k <- which.max(m$nlp_c)
  res[[length(res) + 1]] <- data.table(
    gene = g, accession = acc, n_shared = nrow(m),
    max_nlp_plasma = round(max(m$nlp_p), 2), max_nlp_csf = round(max(m$nlp_c), 2),
    csf_lead = m$rsid[k], csf_lead_nlp = round(m$nlp_c[k], 2),
    plasma_nlp_at_csf_lead = round(m$nlp_p[k], 2),
    PP.H0 = s$PP.H0.abf, PP.H1 = s$PP.H1.abf, PP.H2 = s$PP.H2.abf,
    PP.H3 = s$PP.H3.abf, PP.H4 = s$PP.H4.abf)
  say(sprintf("%-8s %s  shared %4d  plasma -log10P %5.1f  CSF %5.1f  H4 %.3f",
              g, acc, nrow(m), max(m$nlp_p), max(m$nlp_c), s$PP.H4.abf))
}
R <- rbindlist(res)
setorder(R, -PP.H4)
fwrite(R, file.path(D, "coloc_plasma_vs_csf.tsv"), sep = "\t")
say("")
say("coloc tests:", nrow(R), "| PP.H4 > 0.5:", R[PP.H4 > 0.5, .N],
    "| PP.H4 > 0.8:", R[PP.H4 > 0.8, .N])
print(R[, .(gene, n_shared, plasma = max_nlp_plasma, CSF = max_nlp_csf,
            csf_lead, H3 = round(PP.H3, 3), H4 = round(PP.H4, 3))])
say("done")

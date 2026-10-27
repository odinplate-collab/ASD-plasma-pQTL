# Genotype effect by group for one protein and its external lead instrument (by default the protein with
# the strongest genotype-by-diagnosis interaction in genotype_by_diagnosis_interaction.R).
# Individuals: all with whole-genome sequencing and TMT data; abundance from the protein abundance matrix
# (z units); genotype = number of protein-raising alleles of the Niu et al. lead instrument.
#   slopes      abundance ~ genotype within ASD, TD, fathers and mothers (OLS; 95% CI from the t distribution)
#   genotype    ASD versus TD within each genotype class (Welch t-test)
#   alleles     raising-allele frequency in ASD versus TD (Fisher's exact test on allele counts; all
#               sequenced ASD and TD individuals)
#   homozygotes among individuals homozygous for the raising allele: abundance ~ age in those without ASD (OLS),
#               and abundance ~ ASD + male + age + (1 | family) (lme4/lmerTest)
# Inputs (relative to ASD_ROOT):
#   derived/genotype_x_diagnosis/niu_lead_instruments.tsv, interaction_by_specification.csv
#   TMT/250429/ASD_250701_norm_averaged_mat.r, TMT/Peptide/ASD_250822_pqtl_pep_unique_filtered.r,
#   TMT/Sample_info_no_replicate.txt
#   argument 1: per-individual genotype table written by extract_hardcall_genotypes.py
#   argument 2 (optional): protein symbol
# Outputs (relative to ASD_ROOT), summary statistics only: derived/genotype_x_diagnosis/lead_variant_by_group_<protein>.csv
# Usage: Rscript lead_variant_by_group.R <genotype_table.tsv> [protein]
suppressPackageStartupMessages({library(data.table); library(lme4); library(lmerTest)})
ROOT <- Sys.getenv("ASD_ROOT", ".")
OUT  <- file.path(ROOT, "derived/genotype_x_diagnosis"); ar <- commandArgs(TRUE)
prot <- if (length(ar) >= 2) ar[2] else fread(file.path(OUT, "interaction_by_specification.csv"))[spec == "base"][order(gx_p)][1, protein]

G <- fread(ar[1]); G[, wgs := gsub("-", "_", wgs_id)]
niu <- fread(file.path(OUT, "niu_lead_instruments.tsv"))[Protein == prot]; stopifnot(nrow(niu) == 1)
niu[, c("chr", "pos") := tstrsplit(chrpos, ":")]
snp <- grep(paste0("^", niu$chr, ":", niu$pos, ":(", niu$ea, ":", niu$oa, "|", niu$oa, ":", niu$ea, ")$"), names(G), value = TRUE); stopifnot(length(snp) == 1)
alt <- strsplit(snp, ":")[[1]][4]
flip <- (niu$ea == alt & niu$beta_niu < 0) | (niu$ea != alt & niu$beta_niu > 0)     # TRUE: the reference allele raises abundance
cat("protein:", prot, "| instrument:", niu$rsID, snp, "| raising allele:", if (flip) strsplit(snp, ":")[[1]][3] else alt, "\n")

e <- new.env(); load(file.path(ROOT, "TMT/Peptide/ASD_250822_pqtl_pep_unique_filtered.r"), envir = e); pep_ids <- colnames(e$averaged_mat); rm(e)
e2 <- new.env(); load(file.path(ROOT, "TMT/250429/ASD_250701_norm_averaged_mat.r"), envir = e2); PM <- as.matrix(e2$averaged_mat)
si <- fread(file.path(ROOT, "TMT/Sample_info_no_replicate.txt"), select = c("Sample", "Sex", "Age", "Family", "Family_n", "WGS_ID"))
si[, role := fifelse(Family == "sibiling", "Sibling", fifelse(Family == "non-ASD", "TD", Family))]
si[, fam := fifelse(is.na(Family_n), paste0("TD_", Sample), paste0("F", Family_n))]
si[, `:=`(dx = as.integer(role == "ASD"), male = as.integer(Sex == "M"), age10 = Age / 10, wgs = gsub("-", "_", WGS_ID))]
d0 <- merge(si[Sample %in% pep_ids], G[, .(wgs, g = get(snp))], by = "wgs")[!is.na(g)]
if (flip) d0[, g := 2 - g]
d0[, y := PM[prot, match(Sample, colnames(PM))]]; d <- d0[!is.na(y)]
cat(sprintf("individuals with genotype: %d | with abundance: %d (%s)\n", nrow(d0), nrow(d), paste(names(table(d$role)), table(d$role), collapse = ", ")))

row <- function(analysis, group, n, est, se = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, note = "")
  data.table(protein = prot, instrument = niu$rsID, analysis, group, n, estimate = est, se, lo, hi, p, note)
res <- list()
for (r in c("ASD", "TD", "Father", "Mother")) {
  x <- d[role == r]; s <- coef(summary(lm(y ~ g, x)))["g", ]; q <- qt(0.975, nrow(x) - 2)
  res[[length(res) + 1]] <- row("slope per raising allele", r, nrow(x), s[[1]], s[[2]], s[[1]] - q * s[[2]], s[[1]] + q * s[[2]], s[[4]]) }
for (k in 0:2) {
  a <- d[role == "ASD" & g == k, y]; b <- d[role == "TD" & g == k, y]
  res[[length(res) + 1]] <- row("ASD minus TD within genotype", paste(k, "raising alleles"), length(a) + length(b), mean(a) - mean(b),
                                p = if (length(a) > 1 && length(b) > 1) t.test(a, b)$p.value else NA_real_, note = sprintf("%d ASD, %d TD", length(a), length(b))) }
ac <- d0[role %in% c("ASD", "TD"), .(raise = sum(g), other = sum(2 - g), n = .N), by = role]; setkey(ac, role)
ft <- fisher.test(matrix(c(ac["ASD", raise], ac["ASD", other], ac["TD", raise], ac["TD", other]), 2, byrow = TRUE))
res[[length(res) + 1]] <- row("raising-allele frequency", "ASD", ac["ASD", n], ac["ASD", raise / (2 * n)], p = ft$p.value, note = "Fisher's exact test, ASD versus TD")
res[[length(res) + 1]] <- row("raising-allele frequency", "TD", ac["TD", n], ac["TD", raise / (2 * n)])
hz <- d[g == 2]; na <- hz[dx == 0]; s <- coef(summary(lm(y ~ Age, na)))["Age", ]
res[[length(res) + 1]] <- row("homozygotes: abundance ~ age (per year)", "without ASD", nrow(na), s[[1]], s[[2]], p = s[[4]])
m <- suppressMessages(lmer(y ~ dx + male + age10 + (1 | fam), data = hz)); s <- coef(summary(m))["dx", ]
res[[length(res) + 1]] <- row("homozygotes: ASD effect (sex, age, family)", "all homozygotes", nrow(hz), s[[1]], s[[2]], s[[1]] - 1.96 * s[[2]], s[[1]] + 1.96 * s[[2]], s[[ncol(coef(summary(m)))]])
res <- rbindlist(res)
fwrite(res, file.path(OUT, sprintf("lead_variant_by_group_%s.csv", prot)))
options(width = 200); print(res[, .(analysis, group, n, est = round(estimate, 3), lo = round(lo, 2), hi = round(hi, 2), p = signif(p, 2), note)])

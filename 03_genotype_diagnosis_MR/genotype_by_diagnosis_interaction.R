# Genotype-by-diagnosis interaction screen with external instruments, and robustness of the strongest
# interaction.
# Instruments: the lead cis instrument of Niu et al. for every measured protein (extract_hardcall_genotypes.py),
#   kept when its minor allele frequency in this cohort is >= 0.05 and both ASD and the other individuals
#   carry more than one genotype class. Genotype = number of protein-raising alleles in Niu et al.
# Individuals: all with whole-genome sequencing and TMT data (ASD, TD, fathers, mothers, siblings);
#   ASD = 1, everyone else = 0. Abundance: protein abundance matrix (z units).
# Specifications (estimate = genotype x ASD, SD of abundance per raising allele), for every protein:
#   base   y ~ genotype * ASD + male + age + (1 | family)       (screen model; lme4/lmerTest, Satterthwaite df)
#   td89   ASD and TD only, OLS, y ~ genotype * ASD, y standardized across them
#   g_age  base + genotype x age
#   g_sex  base + genotype x sex
#   int    base on rank-based inverse-normal transformed abundance
# For the protein with the smallest base-model interaction P: the genotype x age and genotype x sex terms,
#   and the base model leaving out one family at a time (each unrelated TD individual is its own family).
# Inputs (relative to ASD_ROOT):
#   derived/genotype_x_diagnosis/niu_lead_instruments.tsv  (extract_hardcall_genotypes.py)
#   TMT/250429/ASD_250701_norm_averaged_mat.r              protein abundance matrix (object averaged_mat)
#   TMT/Peptide/ASD_250822_pqtl_pep_unique_filtered.r      peptide matrix (individuals with peptide data)
#   TMT/Sample_info_no_replicate.txt                       sample information
#   argument 1: per-individual genotype table written by extract_hardcall_genotypes.py
# Outputs (relative to ASD_ROOT), summary statistics only: derived/genotype_x_diagnosis/
#   interaction_by_specification.csv, top_protein_robustness.csv, top_protein_competing_terms.csv,
#   top_protein_leave_one_family_out.csv
# Usage: Rscript genotype_by_diagnosis_interaction.R <genotype_table.tsv>
suppressPackageStartupMessages({library(data.table); library(lme4); library(lmerTest)})
ROOT <- Sys.getenv("ASD_ROOT", ".")
OUT  <- file.path(ROOT, "derived/genotype_x_diagnosis"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

G <- fread(commandArgs(TRUE)[1]); G[, wgs := gsub("-", "_", wgs_id)]
niu <- fread(file.path(OUT, "niu_lead_instruments.tsv")); niu[, c("chr", "pos") := tstrsplit(chrpos, ":")]
cols <- data.table(snp = setdiff(names(G), c("wgs_id", "wgs"))); cols[, c("chr", "pos", "ref", "alt") := tstrsplit(snp, ":")]
nv <- merge(niu, cols, by = c("chr", "pos"))[(ea == alt & oa == ref) | (ea == ref & oa == alt)]
nv[, flip := (ea == alt & beta_niu < 0) | (ea == ref & beta_niu > 0)]          # TRUE: the reference allele raises abundance

e <- new.env(); load(file.path(ROOT, "TMT/Peptide/ASD_250822_pqtl_pep_unique_filtered.r"), envir = e); pep_ids <- colnames(e$averaged_mat); rm(e)
e2 <- new.env(); load(file.path(ROOT, "TMT/250429/ASD_250701_norm_averaged_mat.r"), envir = e2); PM <- as.matrix(e2$averaged_mat)
si <- fread(file.path(ROOT, "TMT/Sample_info_no_replicate.txt"), select = c("Sample", "Sex", "Age", "Family", "Family_n", "WGS_ID"))
si[, role := fifelse(Family == "sibiling", "Sibling", fifelse(Family == "non-ASD", "TD", Family))]
si[, fam := fifelse(is.na(Family_n), paste0("TD_", Sample), paste0("F", Family_n))]
si[, `:=`(dx = as.integer(role == "ASD"), male = as.integer(Sex == "M"), age10 = Age / 10, wgs = gsub("-", "_", WGS_ID))]
si <- si[Sample %in% pep_ids & wgs %in% G$wgs]
cat(sprintf("individuals: %d | %s | family units: %d\n", nrow(si), paste(names(table(si$role)), table(si$role), collapse = ", "), uniqueN(si$fam)))

build <- function(r) {
  d <- merge(si, G[, .(wgs, g = get(r$snp))], by = "wgs")[!is.na(g)]
  d[, y := PM[r$Protein, match(Sample, colnames(PM))]]; d <- d[!is.na(y)]
  if (r$flip) d[, g := 2 - g]
  if (min(mean(d$g), 2 - mean(d$g)) / 2 < 0.05 || d[dx == 1, uniqueN(g)] < 2 || d[dx == 0, uniqueN(g)] < 2) return(NULL)
  d }
nv <- nv[Protein %in% rownames(PM)]
D <- setNames(lapply(seq_len(nrow(nv)), function(i) build(nv[i])), nv$Protein); D <- D[!vapply(D, is.null, TRUE)]
cat("proteins in the screen:", length(D), "\n")

fit <- function(f, d) suppressMessages(suppressWarnings(lmer(f, data = d)))
co  <- function(m, k) { s <- coef(summary(m)); if (k %in% rownames(s)) unname(s[k, c(1, 2, ncol(s))]) else rep(NA_real_, 3) }
inorm <- function(x) qnorm((rank(x) - 0.5) / length(x))
SPEC <- list(
  base  = function(d) co(fit(y ~ g * dx + male + age10 + (1 | fam), d), "g:dx"),
  td89  = function(d) { x <- d[role %in% c("ASD", "TD")]; x[, y := (y - mean(y)) / sd(y)]; unname(coef(summary(lm(y ~ g * dx, x)))["g:dx", c(1, 2, 4)]) },
  g_age = function(d) co(fit(y ~ g * dx + g * age10 + male + (1 | fam), d), "g:dx"),
  g_sex = function(d) co(fit(y ~ g * dx + g * male + age10 + (1 | fam), d), "g:dx"),
  int   = function(d) { d <- copy(d)[, y := inorm(y)]; co(fit(y ~ g * dx + male + age10 + (1 | fam), d), "g:dx") })

res <- rbindlist(lapply(names(SPEC), function(s) rbindlist(lapply(names(D), function(p) {
  v <- SPEC[[s]](D[[p]]); data.table(spec = s, protein = p, n = if (s == "td89") D[[p]][role %in% c("ASD", "TD"), .N] else nrow(D[[p]]), gx = v[1], gx_se = v[2], gx_p = v[3]) }))))
res[, rank := frank(gx_p, ties.method = "min"), by = spec]; res[, n_prot := .N, by = spec]
res[spec == "base", fdr := p.adjust(gx_p, "BH")]
fwrite(res, file.path(OUT, "interaction_by_specification.csv"))
top <- res[spec == "base"][order(gx_p)][1, protein]; dm <- D[[top]]
cat("strongest base-model interaction:", top, nv[Protein == top, paste(rsID, snp)], "\n")

# genotype x age and genotype x sex terms for the top protein
mA <- fit(y ~ g * dx + g * age10 + male + (1 | fam), dm); mS <- fit(y ~ g * dx + g * male + age10 + (1 | fam), dm)
comp <- data.table(protein = top, term = c("genotype x age (per 10 years)", "genotype x male sex"),
                   est = c(co(mA, "g:age10")[1], co(mS, "g:male")[1]), se = c(co(mA, "g:age10")[2], co(mS, "g:male")[2]), p = c(co(mA, "g:age10")[3], co(mS, "g:male")[3]))
fwrite(comp, file.path(OUT, "top_protein_competing_terms.csv"))

# leave one family out (all proteins refitted, so that the rank of the top protein can be tracked)
fams <- unique(dm$fam)
lofo <- rbindlist(lapply(fams, function(f) {
  r <- rbindlist(lapply(names(D), function(p) { v <- SPEC$base(D[[p]][fam != f]); data.table(protein = p, gx = v[1], gx_se = v[2], gx_p = v[3]) }))
  r[, rank := frank(gx_p, ties.method = "min")]
  cbind(dropped = f, dropped_role = paste(sort(unique(dm[fam == f, role])), collapse = "+"), r[protein == top], second = r[order(gx_p)][protein != top][1, protein]) }))
fwrite(lofo, file.path(OUT, "top_protein_leave_one_family_out.csv"))

sm <- res[protein == top, .(protein, spec, n, n_prot, est = gx, lo = gx - 1.96 * gx_se, hi = gx + 1.96 * gx_se, p = gx_p, rank_best = rank, rank_worst = rank)]
sm <- rbind(sm, lofo[, .(protein = top, spec = "leave_one_family_out", n = nrow(dm) - 1L, n_prot = length(D), est = NA_real_, lo = min(gx), hi = max(gx), p = max(gx_p), rank_best = min(rank), rank_worst = max(rank))])
fwrite(sm, file.path(OUT, "top_protein_robustness.csv"))
options(width = 200); print(sm, digits = 3); print(comp, digits = 3)

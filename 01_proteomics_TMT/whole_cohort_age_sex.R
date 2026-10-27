# Age and sex sensitivity of differential abundance in the whole cohort.
# Protein abundance matrix (512 proteins x 220 individuals: ASD, TD, fathers, mothers, siblings), z units.
# Per protein, linear mixed models with a random intercept for family (each unrelated TD individual
# forms its own family; lme4/lmerTest, REML), ASD individuals compared with all other individuals:
#   adjusted (linear age):  abundance ~ ASD + male + age (per 10 years)
#   adjusted (spline age):  abundance ~ ASD + male + natural spline of age (3 df)
#   without ASD (n = 167):  abundance ~ male + age (per 10 years)
# plus the unadjusted ASD versus TD difference (two-sided Welch t-test) and, within the families
# with an unaffected sibling, the mean ASD minus sibling difference (one-sample t-test across families).
# FDR: Benjamini-Hochberg across proteins, per coefficient.
# Inputs (relative to ASD_ROOT):
#   TMT/250429/ASD_250701_norm_averaged_mat.r   protein abundance matrix (object averaged_mat)
#   TMT/Sample_info_no_replicate.txt             sample information (group, family, sex, age)
# Outputs (relative to ASD_ROOT): derived/whole_cohort_age_sex/per_protein_whole_cohort.csv
# Usage: Rscript whole_cohort_age_sex.R
suppressPackageStartupMessages({library(data.table); library(lme4); library(lmerTest); library(splines)})
ROOT <- Sys.getenv("ASD_ROOT", ".")
OUT  <- file.path(ROOT, "derived/whole_cohort_age_sex"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

load(file.path(ROOT, "TMT/250429/ASD_250701_norm_averaged_mat.r"))     # averaged_mat (proteins x individuals)
si <- fread(file.path(ROOT, "TMT/Sample_info_no_replicate.txt"))
si <- si[Sample %in% colnames(averaged_mat)]
si[, role := fifelse(Family == "sibiling", "Sibling", fifelse(Family == "non-ASD", "TD", Family))]
si[, fam := fifelse(is.na(Family_n), paste0("TD_", Sample), paste0("F", Family_n))]
si[, `:=`(asd = as.integer(role == "ASD"), male = as.integer(Sex == "M"), age10 = Age / 10)]
X <- as.matrix(averaged_mat[, si$Sample]); stopifnot(identical(colnames(X), si$Sample), !anyNA(X))
cat(sprintf("proteins %d | individuals %d | %s\n", nrow(X), ncol(X), paste(names(table(si$role)), table(si$role), collapse = ", ")))

co  <- function(m, k) { s <- coef(summary(m)); if (k %in% rownames(s)) unname(s[k, c("Estimate", "Std. Error", "Pr(>|t|)")]) else rep(NA_real_, 3) }
fit <- function(f, d) suppressMessages(suppressWarnings(lmer(f, data = d, REML = TRUE)))
non  <- si$asd == 0
sibf <- si[role == "Sibling", unique(fam)]

res <- rbindlist(lapply(seq_len(nrow(X)), function(i) {
  d <- copy(si); d[, y := X[i, ]]
  w   <- t.test(d[role == "ASD", y], d[role == "TD", y])
  mB  <- fit(y ~ male + age10 + (1 | fam), d[non])
  mC1 <- fit(y ~ asd + male + age10 + (1 | fam), d)
  mC2 <- fit(y ~ asd + male + ns(Age, df = 3) + (1 | fam), d)
  pr  <- d[fam %in% sibf, .(asd = mean(y[role == "ASD"]), sib = mean(y[role == "Sibling"])), by = fam][!is.na(asd) & !is.na(sib)]
  tp  <- t.test(pr$asd - pr$sib)
  mu  <- d[, .(m = mean(y)), by = role]; g <- setNames(mu$m, mu$role)
  a1 <- co(mC1, "asd"); s1 <- co(mC1, "male"); g1 <- co(mC1, "age10"); a2 <- co(mC2, "asd"); sB <- co(mB, "male"); gB <- co(mB, "age10")
  data.table(Gene = rownames(X)[i],
    mean_ASD = g[["ASD"]], mean_Sibling = g[["Sibling"]], mean_TD = g[["TD"]], mean_Father = g[["Father"]], mean_Mother = g[["Mother"]],
    unadj_asd_td = unname(w$estimate[1] - w$estimate[2]), unadj_p = w$p.value,
    asd_linear = a1[1], asd_linear_se = a1[2], asd_linear_p = a1[3],
    asd_spline = a2[1], asd_spline_se = a2[2], asd_spline_p = a2[3],
    male_same_model = s1[1], male_same_model_se = s1[2], male_same_model_p = s1[3],
    age10_same_model = g1[1], age10_same_model_se = g1[2], age10_same_model_p = g1[3],
    male_without_asd = sB[1], male_without_asd_p = sB[3], age10_without_asd = gB[1], age10_without_asd_p = gB[3],
    asd_minus_sibling = mean(pr$asd - pr$sib), asd_minus_sibling_p = tp$p.value, sibling_families = nrow(pr))
}))
for (k in grep("_p$", names(res), value = TRUE)) res[, (sub("_p$", "_fdr", k)) := p.adjust(get(k), "BH")]
fwrite(res, file.path(OUT, "per_protein_whole_cohort.csv"))
cat(sprintf("ASD effect (linear age) P < 0.05: %d, FDR < 0.05: %d of %d proteins\n", res[asd_linear_p < 0.05, .N], res[asd_linear_fdr < 0.05, .N], nrow(res)))

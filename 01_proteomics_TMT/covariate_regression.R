# Sex and age effects on protein abundance (limma; covariates: diagnosis, TMT batch, TMT channel; duplicate-correlation block)
# Inputs (relative to ASD_ROOT):
#   TMT/250429/ASD_250701_norm_averaged_mat.r   (R workspace: averaged_mat = protein z-scores, proteins x samples;
#                                                sample.info = TMT sample sheet with Sample, Sex, Age, Group, Batch, Channel, Family)
# Outputs (relative to ASD_ROOT):
#   output/covariate_regression_sex.txt         (Gene, beta, pval, FDR, t, direction)
#   output/covariate_regression_age.txt         (Gene, beta, pval, FDR, t, direction)
#   output/covariate_regression_thresholds.txt  (|beta| threshold per covariate)
# Usage: Rscript 01_proteomics_TMT/covariate_regression.R

ROOT <- Sys.getenv("ASD_ROOT", ".")
OUT  <- file.path(ROOT, "output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

library(limma)
library(tidyverse)

load(file.path(ROOT, "TMT/250429/ASD_250701_norm_averaged_mat.r"))   # -> averaged_mat, sample.info

# ==== 0) Match expression matrix to sample metadata ====
stopifnot(is.matrix(averaged_mat) || is.data.frame(averaged_mat))
expr <- as.matrix(averaged_mat)

stopifnot("Sample" %in% colnames(sample.info))
meta <- sample.info %>% mutate(Sample = as.character(Sample))
common <- intersect(colnames(expr), meta$Sample)
if (length(common) < 3) stop("Sample matching failed: check that colnames(averaged_mat) match sample.info$Sample.")

expr <- expr[, common, drop = FALSE]
meta <- meta %>% filter(Sample %in% common) %>% arrange(match(Sample, colnames(expr)))

if (is.null(rownames(expr))) stop("averaged_mat must have gene names as rownames.")
genes <- rownames(expr)

# ==== 1) Covariates: diagnosis, batch, channel, block ====
diag_col <- if ("Group" %in% names(meta)) "Group" else if ("ASD" %in% names(meta)) "ASD" else NA_character_
if (is.na(diag_col)) stop("Diagnosis column not found: sample.info must contain Group or ASD.")

diag_factor <- factor(
  dplyr::recode(as.character(meta[[diag_col]]),
                "ASD" = "ASD", "asd" = "ASD",
                "non-ASD" = "non-ASD", "control" = "non-ASD",
                "typical development" = "non-ASD", .default = as.character(meta[[diag_col]])),
  levels = c("non-ASD", "ASD")
)

sex_factor <- factor(meta$Sex, levels = c("F", "M"))  # reference F; coefficient sex_factorM
if (any(is.na(sex_factor))) stop("Sex must be coded as 'F' / 'M'.")

age_num <- as.numeric(meta$Age)
if (any(is.na(age_num))) stop("Age must be numeric.")
age_sc  <- scale(age_num)  # standardized age

batch_factor   <- factor(meta$Batch)
channel_factor <- factor(meta$Channel)

# Block for duplicateCorrelation: the "Family" column (family role in the sample sheet); "Family_n" if absent
fam_col <- if ("Family" %in% names(meta)) "Family" else if ("Family_n" %in% names(meta)) "Family_n" else NA_character_
if (is.na(fam_col)) stop("A Family or Family_n column is required.")
family_factor <- factor(meta[[fam_col]])

# ==== 2) Design matrix ====
design <- model.matrix(~ sex_factor + age_sc + diag_factor + batch_factor + channel_factor)
colnames(design) <- make.names(colnames(design))
coef_sex <- "sex_factorM"
coef_age <- "age_sc"

# ==== 3) Duplicate correlation and model fit ====
dupcor <- duplicateCorrelation(expr, design, block = family_factor)
fit <- lmFit(expr, design, block = family_factor, correlation = dupcor$consensus)
fit <- eBayes(fit)

# ==== 4) Results ====
extract_coef <- function(fit, coef_name) {
  tt <- topTable(fit, coef = coef_name, n = Inf, sort.by = "none")
  tt %>%
    as_tibble(rownames = "Gene") %>%
    transmute(Gene,
              beta = logFC,        # for a continuous covariate the limma coefficient is the regression slope
              pval = P.Value,
              FDR  = adj.P.Val,
              t    = t)
}

res_sex <- extract_coef(fit, coef_sex)
res_age <- extract_coef(fit, coef_age)

# ==== 5) Significance calls: FDR < 0.05 and |beta| >= threshold ====
# Threshold = 95th percentile of |beta| divided by 3 (fallback: 75th percentile / 2, then 0.1)
call_direction <- function(df, alpha = 0.05) {
  dat <- df %>%
    filter(is.finite(beta), is.finite(FDR))
  beta_cut <- suppressWarnings(quantile(abs(dat$beta), 0.95, na.rm = TRUE) / 3)
  if (!is.finite(beta_cut) || beta_cut <= 0)
    beta_cut <- suppressWarnings(quantile(abs(dat$beta), 0.75, na.rm = TRUE) / 2)
  if (!is.finite(beta_cut) || beta_cut <= 0) beta_cut <- 0.1
  p_use <- pmin(pmax(df$FDR, 1e-300), 1)
  df$direction <- case_when(
    p_use < alpha & df$beta >=  beta_cut ~ "Up",
    p_use < alpha & df$beta <= -beta_cut ~ "Down",
    TRUE ~ "NS")
  list(res = df, beta_cut = unname(beta_cut))
}

sex_call <- call_direction(res_sex)
age_call <- call_direction(res_age)

write.table(sex_call$res, file.path(OUT, "covariate_regression_sex.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(age_call$res, file.path(OUT, "covariate_regression_age.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
thresholds <- data.frame(covariate = c("sex", "age"),
                         beta_cut = c(sex_call$beta_cut, age_call$beta_cut),
                         n_up   = c(sum(sex_call$res$direction == "Up"),   sum(age_call$res$direction == "Up")),
                         n_down = c(sum(sex_call$res$direction == "Down"), sum(age_call$res$direction == "Down")),
                         duplicate_correlation = dupcor$consensus)
write.table(thresholds, file.path(OUT, "covariate_regression_thresholds.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
print(thresholds)

# ============================
# Sex & Age effect regression on protein abundance (limma)
#   - Family correlation handled with duplicateCorrelation
#   - Covariates: diagnosis group, TMT batch, channel
#   - Outputs: res_sex, res_age (Gene, beta, pval, FDR, t)
# ============================

# ---- pkgs ----
# install.packages(c("limma","edgeR","tidyverse","ggrepel","patchwork"))
# Paths are relative to the project data root.
# Set the ASD_ROOT environment variable, or run R with that directory as the working directory.
setwd(Sys.getenv("ASD_ROOT", "."))

library(limma)
library(tidyverse)
library(ggrepel)
library(patchwork)

# ==== 0) Match expression matrix to sample metadata ====
# averaged_mat: rows = genes, cols = samples
# sample.info:  sample-level metadata (must include Sample, Sex, Age, Group, Batch, Channel, Family)
stopifnot(is.matrix(averaged_mat) || is.data.frame(averaged_mat))
expr <- as.matrix(averaged_mat)

# match column names to sample IDs
stopifnot("Sample" %in% colnames(sample.info))
meta <- sample.info %>% mutate(Sample = as.character(Sample))
common <- intersect(colnames(expr), meta$Sample)
if (length(common) < 3) stop("Sample matching failed: check that colnames(averaged_mat) match sample.info$Sample.")

# order samples
expr <- expr[, common, drop = FALSE]
meta <- meta %>% filter(Sample %in% common) %>% arrange(match(Sample, colnames(expr)))

# gene names
if (is.null(rownames(expr))) stop("averaged_mat must have gene names as rownames.")
genes <- rownames(expr)

# ==== 1) Covariates: diagnosis, batch, channel, family ====
# detect diagnosis column: prefer 'Group', otherwise 'ASD'
diag_col <- if ("Group" %in% names(meta)) "Group" else if ("ASD" %in% names(meta)) "ASD" else NA_character_
if (is.na(diag_col)) stop("Diagnosis column not found: sample.info must contain Group or ASD.")

# normalize labels to "ASD" / "non-ASD"
diag_factor <- factor(
  dplyr::recode(as.character(meta[[diag_col]]),
                "ASD" = "ASD", "asd" = "ASD",
                "non-ASD" = "non-ASD", "control" = "non-ASD",
                "typical development" = "non-ASD", .default = as.character(meta[[diag_col]])),
  levels = c("non-ASD","ASD")
)

sex_factor <- factor(meta$Sex, levels = c("F","M"))  # reference level F (coefficient is named SexM)
if (any(is.na(sex_factor))) stop("Sex must be coded as 'F' / 'M'.")

age_num <- as.numeric(meta$Age)
if (any(is.na(age_num))) stop("Age must be numeric.")
age_sc  <- scale(age_num)  # standardize

batch_factor   <- factor(meta$Batch)
channel_factor <- factor(meta$Channel)

# family ID (block for duplicate correlation)
fam_col <- if ("Family" %in% names(meta)) "Family" else if ("Family_n" %in% names(meta)) "Family_n" else NA_character_
if (is.na(fam_col)) stop("A Family or Family_n column is required.")
family_factor <- factor(meta[[fam_col]])

# ==== 2) Design matrix ====
design <- model.matrix(~ sex_factor + age_sc + diag_factor + batch_factor + channel_factor)
colnames(design) <- make.names(colnames(design))
# record coefficient names of interest
coef_sex <- "sex_factorM"
coef_age <- "age_sc"

# ==== 3) Estimate family duplicate-correlation and fit model ====
dupcor <- duplicateCorrelation(expr, design, block = family_factor)
fit <- lmFit(expr, design, block = family_factor, correlation = dupcor$consensus)
fit <- eBayes(fit)

# ==== 4) Result extraction ====
extract_coef <- function(fit, coef_name) {
  tt <- topTable(fit, coef = coef_name, n = Inf, sort.by = "none")
  out <- tt %>%
    as_tibble(rownames = "Gene") %>%
    transmute(Gene,
              beta = logFC,                   # for continuous covariates the limma coef equals the regression slope
              pval = P.Value,
              FDR  = adj.P.Val,
              t    = t)
  out
}

res_sex <- extract_coef(fit, coef_sex)
res_age <- extract_coef(fit, coef_age)

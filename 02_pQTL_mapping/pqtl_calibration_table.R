# Observed versus permuted association counts at each P-value threshold, from the calibration
# runs (pqtl_calibration_permutation.R: observed = seed 0, permuted = seeds 1-100).
# For each class (cis, trans) and threshold: observed count, mean and 2.5th-97.5th percentiles
# of the permuted counts, enrichment = observed / mean permuted, and empirical FDR (%) =
# 100 x mean permuted / observed.
# FDR 5% cut-offs of the main mapping: cis P <= 1.024e-5, trans P <= 4.445e-7.
# Inputs (relative to ASD_ROOT): calibration/outputs/fdrcurve_seed<0..100>.csv
# Outputs (relative to ASD_ROOT): calibration/outputs/calibration_permutation_summary.csv
# Usage: Rscript pqtl_calibration_table.R
suppressMessages(library(data.table))
ROOT <- Sys.getenv("ASD_ROOT", ".")
OUT  <- file.path(ROOT, "calibration/outputs")
files <- list.files(OUT, pattern = "^fdrcurve_seed[0-9]+\\.csv$", full.names = TRUE)
d <- rbindlist(lapply(files, fread))
d[, cum := as.numeric(cum)]
obs  <- d[seed == 0, .(class, p_threshold = p_upper, observed = cum)]
perm <- d[seed > 0, .(n_permutations = .N, null_mean = mean(cum),
                      null_lo = quantile(cum, 0.025), null_hi = quantile(cum, 0.975)),
          by = .(class, p_threshold = p_upper)]
m <- merge(obs, perm, by = c("class", "p_threshold"))[p_threshold <= 1e-4]
m[, enrichment := fifelse(null_mean > 0, observed / null_mean, NA_real_)]
m[, empirical_FDR_pct := fifelse(observed > 0, 100 * null_mean / observed, NA_real_)]
setorder(m, class, p_threshold)
fwrite(m, file.path(OUT, "calibration_permutation_summary.csv"))
print(m[p_threshold %in% c(1.024e-5, 4.445e-7)])

# Export of the protein abundance matrix (proteins x individuals, z units) to a tab-separated
# table, for the Python genotype and diagnosis models.
# Inputs (relative to ASD_ROOT): TMT/250429/ASD_250701_norm_averaged_mat.r (object averaged_mat)
# Outputs (relative to ASD_ROOT): derived/protein_matrix/protein_abundance_matrix.tsv.gz (rows: gene; columns: TMT sample)
# Usage: Rscript export_protein_matrix.R
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("ASD_ROOT", ".")
OUT  <- file.path(ROOT, "derived/protein_matrix"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
e <- new.env(); load(file.path(ROOT, "TMT/250429/ASD_250701_norm_averaged_mat.r"), envir = e)
PM <- as.matrix(e$averaged_mat)
fwrite(data.table(gene = rownames(PM), PM), file.path(OUT, "protein_abundance_matrix.tsv.gz"), sep = "\t")
cat("proteins:", nrow(PM), "| individuals:", ncol(PM), "\n")

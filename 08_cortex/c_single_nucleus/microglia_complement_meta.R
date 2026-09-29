# Microglial C1q, C3 and complement receptors: fixed-effect meta-analysis of the two snRNA-seq datasets (Wamsley, Velmeshev).
# Inputs (relative to ASD_ROOT): brain/data/lr_singlecell/de_genes_both_datasets.tsv  (sc_de_both.R)
# Outputs (relative to ASD_ROOT): brain/data/lr_singlecell/microglia_complement_meta.tsv
#   per gene: logFC / P / FDR in each dataset, meta log2FC, meta SE, meta P, heterogeneity P
# Usage: Rscript microglia_complement_meta.R
# Inverse-variance fixed effect with SE = |logFC / t| from limma-voom; heterogeneity by Cochran's Q (1 df).
# All rows of the microglia (MG) cell type are used for each gene.
suppressPackageStartupMessages({library(data.table)})
ROOT <- Sys.getenv("ASD_ROOT", ".")
SC <- file.path(ROOT, "brain", "data", "lr_singlecell")
GENES <- c("C1QA", "C1QB", "C1QC", "C3", "C3AR1", "VSIG4", "C5AR1", "ITGB2", "ITGAM", "ITGAX")
RECEPTORS <- c("C3AR1", "VSIG4", "C5AR1", "ITGB2", "ITGAM", "ITGAX")   # the six microglial complement receptors

de <- fread(file.path(SC, "de_genes_both_datasets.tsv"))[celltype == "MG" & gene %in% GENES]
de[, se := abs(logFC / t)]
meta <- de[, {w <- 1 / se^2; b <- sum(logFC * w) / sum(w); s <- sqrt(1 / sum(w))
              q <- sum(w * (logFC - b)^2)
              .(n_datasets = .N, meta = b, meta_se = s, meta_p = 2 * pnorm(-abs(b / s)),
                het_p = pchisq(q, 1, lower.tail = FALSE))}, by = gene]
out <- merge(dcast(de, gene ~ dataset, value.var = c("logFC", "P.Value", "adj.P.Val")), meta, by = "gene")
out[, receptor := gene %in% RECEPTORS]
out <- out[match(intersect(GENES, gene), gene)]
fwrite(out, file.path(SC, "microglia_complement_meta.tsv"), sep = "\t")
print(out[, .(gene, receptor, meta = round(meta, 3), meta_p = signif(meta_p, 3), het_p = signif(het_p, 3))])
cat("receptors with meta P < 0.05:", out[receptor & meta_p < 0.05, paste(gene, collapse = " ")], "\n")

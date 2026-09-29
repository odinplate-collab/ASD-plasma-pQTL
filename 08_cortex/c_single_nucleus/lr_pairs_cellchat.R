# Ligand-receptor pairs from CellChatDB.human in which one of the 12 panel proteins is the ligand (complex receptors expanded to subunit genes).
# Inputs (relative to ASD_ROOT): none (CellChatDB.human from the CellChat R package)
# Outputs (relative to ASD_ROOT): brain/data/lr_pairs_panel12_cellchat.tsv
# Usage: Rscript lr_pairs_cellchat.R
suppressPackageStartupMessages({library(CellChat); library(data.table)})
ROOT <- Sys.getenv("ASD_ROOT", ".")
db <- CellChatDB.human
I <- as.data.table(db$interaction)
P12 <- c("AHSG","ANPEP","BCHE","BTD","C1RL","C3","CLEC3B","IGFBP5","MBL2","POSTN","PTGDS","QSOX1")
cx <- as.data.table(db$complex, keep.rownames = "complex")
# a complex name is replaced by its subunit genes
expand <- function(x) {
  if (x %in% cx$complex) { r <- unlist(cx[complex == x, -1]); r[nzchar(r)] } else x
}
hit <- I[sapply(ligand, function(l) any(expand(l) %in% P12))]
out <- hit[, .(interaction_name, ligand, receptor, pathway_name, annotation,
               receptor_genes = sapply(receptor, function(r) paste(expand(r), collapse = "+")))]
fwrite(out, file.path(ROOT, "brain", "data", "lr_pairs_panel12_cellchat.tsv"), sep = "\t")

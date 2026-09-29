# Pseudobulk differential expression (limma-voom) per cortical cell type in two snRNA-seq datasets (Wamsley et al. 2024, Velmeshev et al. 2019) with one identical pipeline, plus expression calls for the 12 panel genes and their receptors.
# Inputs (relative to ASD_ROOT):
#   brain/data/lr_pairs_panel12_cellchat.tsv, lr_pairs_panel12_omnipath.tsv   (lr_pairs_cellchat.R, lr_pairs_omnipath.py)
#   brain/data/wamsley/pb_groups.tsv, pb_counts.tsv.gz                        (wamsley_pseudobulk.py)
#   brain/data/lr_singlecell/wamsley_lr_counts.tsv, wamsley_lr_detect.tsv     (wamsley_lr_fetch.py)
#   brain/data/lr_singlecell/donor_qc_wamsley.tsv                             (donor_qc.py)
#   brain/data/velmeshev/pb_groups.tsv, pb_counts_symbol.tsv.gz               (velmeshev_pseudobulk.py, velmeshev_export.py)
#   brain/data/lr_singlecell/velmeshev_lr_detect.tsv                          (velmeshev_export.py)
# Outputs (relative to ASD_ROOT), all in brain/data/lr_singlecell/:
#   de_genes_both_datasets.tsv               limma-voom topTable (coefficient ASD vs control) per dataset and cell type
#   celltype_info_both_datasets.tsv          donors (ASD / control), genes tested and nuclei per dataset and cell type
#   panel12_by_celltype_both_datasets.tsv    12 panel genes: mean CPM, % nuclei detecting, expression call, DE statistics
#   receptors_by_celltype_both_datasets.tsv  the same for the curated receptor genes
# Usage: Rscript sc_de_both.R
# Pipeline (identical for both datasets): raw-count pseudobulk per donor x broad cell type, library size =
# total UMI of the contributing nuclei, profiles with < 20 nuclei dropped, genes with >= 1 CPM in >= 25% of
# profiles tested with limma-voom, ~ dx + region + sex + age (region dropped when only one is present),
# on a gene universe shared by both datasets. Wamsley donor sex is taken from XIST / Y-gene expression
# (donor_qc.py), falling back to Gandal 2022 and then the browser metadata for ambiguous donors; donors
# listed as Dup15q in Gandal 2022 are excluded.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
ROOT <- Sys.getenv("ASD_ROOT", ".")
D <- file.path(ROOT, "brain", "data"); OUT <- file.path(D, "lr_singlecell")
say <- function(...) { cat(format(Sys.time()), "-", ..., "\n"); flush.console() }
MIN_CELLS <- 20L
DET_MIN <- 0.02          # a gene counts as expressed in a cell type only if
AMB_MIN <- 0.10          # >= 2% of nuclei detect it AND its CPM is >= 10% of the
                         # top cell type (below that, ambient RNA is the likelier source)
P12 <- c("AHSG","ANPEP","BCHE","BTD","C1RL","C3","CLEC3B","IGFBP5","MBL2",
         "POSTN","PTGDS","QSOX1")
CT7 <- c("EXN","INN","AST","ODC","OPC","MG","END")
outf <- function(n) file.path(OUT, paste0(n, ".tsv"))

## curated receptors of the panel proteins: CellChatDB.human + OmniPath, plus IGF1R / IGF2R
## (receptors of the IGFs bound by IGFBP5) and VSIG4 / C5AR1 (complement receptors not in the databases)
LRc <- fread(file.path(D, "lr_pairs_panel12_cellchat.tsv"))
LRo <- fread(file.path(D, "lr_pairs_panel12_omnipath.tsv"))
PAIRS <- unique(rbind(
  LRc[, .(ligand = sub("^PGD2-", "", ligand), receptor = receptor,
          subunits = receptor_genes, source = "CellChatDB")],
  LRo[, .(ligand, receptor, subunits = receptor, source = "OmniPath")]))
PAIRS <- PAIRS[, .(source = paste(sort(unique(source)), collapse = "+")),
               by = .(ligand, receptor, subunits)]
PAIRS <- rbind(PAIRS, data.table(
  ligand = c("IGFBP5", "IGFBP5", "C3", "C3"), receptor = c("IGF1R", "IGF2R", "VSIG4", "C5AR1"),
  subunits = c("IGF1R", "IGF2R", "VSIG4", "C5AR1"),
  source = c("indirect (IGF receptor)", "indirect (IGF receptor)",
             "complement receptor, not in DBs", "C5a receptor, not in DBs")))
REC <- unique(unlist(strsplit(PAIRS$subunits, "+", fixed = TRUE)))

rd <- function(f) {
  x <- if (grepl("gz$", f)) fread(cmd = paste("gzip -dc", shQuote(f)), showProgress = FALSE)
       else fread(f, showProgress = FALSE)
  m <- as.matrix(x[, -1]); rownames(m) <- x$gene; m
}

## ---- Wamsley -------------------------------------------------------------
GW <- fread(file.path(D, "wamsley/pb_groups.tsv"))
MW <- rd(file.path(D, "wamsley/pb_counts.tsv.gz"))
LW <- rd(file.path(OUT, "wamsley_lr_counts.tsv"))
ov <- intersect(rownames(MW), rownames(LW))
stopifnot(max(abs(MW[ov, ] - LW[ov, ])) == 0)
MW <- rbind(MW, LW[setdiff(rownames(LW), rownames(MW)), , drop = FALSE])
DW <- rd(file.path(OUT, "wamsley_lr_detect.tsv"))
GW[, `:=`(donor = as.character(individual_ID),
          dx = fifelse(grepl("^ASD", diagnosis, ignore.case = TRUE), "ASD",
               fifelse(grepl("control|ctl|healthy", diagnosis, ignore.case = TRUE),
                       "Control", NA_character_)))]
## donor QC: sex from expression (Gandal 2022 record, then metadata, when ambiguous); Dup15q excluded
QW <- fread(file.path(OUT, "donor_qc_wamsley.tsv"), colClasses = list(character = "donor"))
QW[, sex_final := fifelse(sex_expr != "ambiguous", sex_expr,
                  fifelse(gandal_sex == "M", "XY", fifelse(gandal_sex == "F", "XX", sex_meta)))]
GW[, sex := QW$sex_final[match(donor, QW$donor)]]
stopifnot(!anyNA(GW$sex))
GW[donor %in% QW[dup15q == 1, donor], dx := NA_character_]
say("sex corrected for", QW[sex_final != sub("XYY", "XY", sex_meta), .N],
    "donors | Dup15q excluded:", paste(QW[dup15q == 1, donor], collapse = " "))

## ---- Velmeshev -----------------------------------------------------------
GV <- fread(file.path(D, "velmeshev/pb_groups.tsv"))
MV <- rd(file.path(D, "velmeshev/pb_counts_symbol.tsv.gz"))
DV <- rd(file.path(OUT, "velmeshev_lr_detect.tsv"))
GV[, `:=`(donor = as.character(individual_ID), dx = diagnosis)]
stopifnot(ncol(MW) == nrow(GW), ncol(MV) == nrow(GV))

## ---- collapse to donor x cell type (identical for both) -------------------
collapse <- function(G, M, Dt) {
  G <- copy(G); G[, key := paste(donor, celltype, sep = "|")]
  ix <- split(seq_len(nrow(G)), G$key)
  sm <- function(X) sapply(ix, function(i) rowSums(X[, i, drop = FALSE]))
  M2 <- sm(M); D2 <- sm(Dt)
  G2 <- G[, .(n_cells = sum(n_cells), lib_size = sum(lib_size),
              region = names(which.max(tapply(n_cells, region, sum))),
              dx = dx[1], age = age[1], sex = sex[1]), by = .(key, donor, celltype)]
  G2 <- G2[match(colnames(M2), key)]
  stopifnot(identical(G2$key, colnames(M2)), identical(colnames(D2), colnames(M2)))
  list(G = G2, M = M2, D = D2)
}
W <- collapse(GW, MW, DW)
V <- collapse(GV, MV, DV)

## shared gene universe (Wamsley was fetched gene by gene)
UNIV <- intersect(rownames(W$M), rownames(V$M))
say("Wamsley genes", nrow(W$M), "| Velmeshev genes", nrow(V$M), "| shared universe", length(UNIV))
say("12 panel genes in universe:", paste(intersect(P12, UNIV), collapse = " "),
    "| missing:", paste(setdiff(P12, UNIV), collapse = " "))

## ---- one cell type, one dataset ------------------------------------------
run_one <- function(ds, X, ct, genes) {
  idx <- which(X$G$celltype == ct)
  g <- X$G[idx]; m <- X$M[genes, idx, drop = FALSE]; dt <- X$D[, idx, drop = FALSE]
  keep_s <- !is.na(g$dx) & g$n_cells >= MIN_CELLS & g$lib_size > 0
  g <- g[keep_s]; m <- m[, keep_s, drop = FALSE]; dt <- dt[, keep_s, drop = FALSE]
  nA <- uniqueN(g[dx == "ASD", donor]); nC <- uniqueN(g[dx == "Control", donor])
  if (nA < 4 || nC < 4) return(NULL)
  g[, `:=`(dx = factor(dx, levels = c("Control", "ASD")), region = factor(region),
           sex = factor(sex), age = suppressWarnings(as.numeric(age)))]
  g[is.na(age), age := median(g$age, na.rm = TRUE)]
  y <- DGEList(counts = m, lib.size = g$lib_size)
  cp <- cpm(y)
  # expression summary for every gene before filtering; detect_frac = nuclei with >= 1 UMI / all nuclei
  lr_g <- intersect(rownames(dt), rownames(m))
  expr <- data.table(dataset = ds, celltype = ct, gene = rownames(m),
                     mean_cpm = rowMeans(cp), frac_donors_cpm1 = rowMeans(cp >= 1))
  expr[, detect_frac := NA_real_]
  expr[match(lr_g, gene), detect_frac := rowSums(dt[lr_g, , drop = FALSE]) / sum(g$n_cells)]
  keep_g <- rowMeans(cp >= 1) >= 0.25
  expr[, passes_filter := gene %in% rownames(m)[keep_g]]
  y <- y[keep_g, , keep.lib.sizes = TRUE]
  form <- if (nlevels(droplevels(g$region)) > 1) ~ dx + region + sex + age else ~ dx + sex + age
  des <- model.matrix(form, data = g)
  v <- voom(y, des)
  fit <- eBayes(lmFit(v, des))
  tt <- as.data.table(topTable(fit, coef = "dxASD", n = Inf), keep.rownames = "gene")
  tt[, `:=`(dataset = ds, celltype = ct)]
  info <- data.table(dataset = ds, celltype = ct, n_ASD = nA, n_CTL = nC,
                     genes_tested = nrow(y), n_cells = sum(g$n_cells))
  list(tt = tt, expr = expr, info = info)
}

res <- list()
for (ds in c("Wamsley", "Velmeshev")) {
  X <- if (ds == "Wamsley") W else V
  cts <- if (ds == "Wamsley") CT7 else c(CT7, "NRGN", "NEUMAT")
  for (ct in cts) {
    r <- run_one(ds, X, ct, union(UNIV, intersect(c(P12, REC), rownames(X$M))))
    if (is.null(r)) { say(ds, ct, "skipped (too few donors)"); next }
    res[[paste(ds, ct)]] <- r
    say(sprintf("%-9s %-6s ASD %2d / CTL %2d  genes %4d  nuclei %6d", ds, ct,
                r$info$n_ASD, r$info$n_CTL, r$info$genes_tested, r$info$n_cells))
  }
}

bind <- function(L, k) rbindlist(lapply(L, `[[`, k), fill = TRUE)
TT <- bind(res, "tt"); EXPR <- bind(res, "expr"); INFO <- bind(res, "info")
fwrite(TT, outf("de_genes_both_datasets"), sep = "\t")
fwrite(INFO, outf("celltype_info_both_datasets"), sep = "\t")

## Expression call: detected in >= 2% of nuclei (DET_MIN)
## and mean CPM >= 10% (AMB_MIN) of the highest-expressing of the seven broad cell types (CT7), among
## genes passing the CPM filter. `expressed` takes the highest cell type over all annotated cell types
## (for Velmeshev also NRGN and NEUMAT); `expressed_2pct` restricts it to CT7 (NA for other cell types).
expr_call <- function(X) {
  X[, top_cpm := max(mean_cpm), by = .(dataset, gene)]
  X[, rel_to_top := mean_cpm / top_cpm]
  X[, expressed := passes_filter & !is.na(detect_frac) & detect_frac >= DET_MIN & rel_to_top >= AMB_MIN]
  X[, top_cpm_ct7 := max(mean_cpm[celltype %in% CT7]), by = .(dataset, gene)]
  X[, expressed_2pct := fifelse(celltype %in% CT7,
                                passes_filter & !is.na(detect_frac) & detect_frac >= DET_MIN &
                                  mean_cpm / top_cpm_ct7 >= AMB_MIN, NA)]
  X
}

## ---- the 12 panel genes, every cell type, both datasets ---------------------
P12T <- merge(EXPR[gene %in% P12],
              TT[gene %in% P12, .(dataset, celltype, gene, logFC, t, P.Value, adj.P.Val)],
              by = c("dataset", "celltype", "gene"), all.x = TRUE)
P12T <- expr_call(P12T)
fwrite(P12T, outf("panel12_by_celltype_both_datasets"), sep = "\t")
say(""); say("=== the 12 genes: logFC (P) per cell type; '.' = not expressed enough to test ===")
P12T[, v := fifelse(expressed, sprintf("%+.2f(%.0e)", logFC, P.Value),
            fifelse(passes_filter, sprintf("amb %+.2f", logFC), "."))]
for (ds in c("Wamsley", "Velmeshev")) {
  say(ds); print(dcast(P12T[dataset == ds], gene ~ celltype, value.var = "v"), width = 250)
}

## ---- receptors of the 12 ----------------------------------------------------
RT <- merge(EXPR[gene %in% REC],
            TT[gene %in% REC, .(dataset, celltype, gene, logFC, t, P.Value, adj.P.Val)],
            by = c("dataset", "celltype", "gene"), all.x = TRUE)
RT <- expr_call(RT)
fwrite(RT, outf("receptors_by_celltype_both_datasets"), sep = "\t")
say(""); say("=== receptors of the 12: expressed where, and ASD change (logFC, P) ===")
RT[, v := fifelse(expressed, sprintf("%+.2f(%.0e)", logFC, P.Value),
          fifelse(passes_filter, sprintf("amb %+.2f", logFC), "."))]
for (ds in c("Wamsley", "Velmeshev")) {
  say(ds); print(dcast(RT[dataset == ds], gene ~ celltype, value.var = "v"), width = 250)
}
say("done")

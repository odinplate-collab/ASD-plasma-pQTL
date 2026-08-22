suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(scales); library(ggnewscale)
})

#setwd(file.path(Sys.getenv("ASD_ROOT", "."), "single_cell/"))

load("asd_single_cell.r")

## =========================
## Single-cell analysis: risk proteins -> receptor readiness -> ligand-receptor network
## =========================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix)
  library(dplyr);  library(tidyr);  library(ggplot2)
  library(reshape2)
  library(sandwich); library(lmtest)
  library(scales)
  library(CellChat) # required for circular network plots
})

## ---- 0) Parameters and gene sets ----
expr_layer   <- "data"      # Seurat v5: layer="data" (log-normalized)
cutoff_expr  <- 0.1         # expression cutoff per cell
min_cells    <- 6L          # minimum cells per Region x cluster x gene
min_donors   <- 2L          # minimum number of donors
alpha_p      <- 0.05        # p-value filter for display
topN_edges   <- 25L         # top edges in the circular diagram

# nine risk-protein ligands (no TGF-beta ligand; TGF-beta excluded from L-R pairs)
risk9 <- c("C3","AHSG","VASN","COL6A1","POSTN","VNN1","CFHR3","C1RL")

# receptor panel
receptors <- list(
  Integrin_ECM = c("ITGAV","ITGA1","ITGA2","ITGA5","ITGB1","ITGB3","ITGB5","DDR1","DDR2","LRP1","SDC1","SDC4"),
  Complement  = c("C3AR1","C5AR1","LRP1","CD46")
)

# map risk proteins to axes for L-R pairing (AHSG, VNN1 unassigned and excluded)
lig_axis <- c(
  POSTN="Integrin_ECM", COL6A1="Integrin_ECM", VASN="Integrin_ECM",
  C3="Complement", CFHR3="Complement", C1RL="Complement",
  AHSG=NA, VNN1=NA
)

## ---- 1) Prepare metadata and expression data ----
stopifnot(inherits(merged, "Seurat"))
DefaultAssay(merged) <- "RNA"

md <- merged@meta.data
# check and normalize required columns (names may vary)
md$Region <- if ("Region" %in% names(md)) md$Region else md$region
md$cluster <- md$cluster
md$donor   <- if ("donor_id" %in% names(md)) md$donor_id else if ("individual" %in% names(md)) md$individual else md$donor
md$group   <- factor(if ("Dx" %in% names(md)) md$Dx else if ("diagnosis" %in% names(md)) md$diagnosis else md$group,
                     levels = c("Control","ASD"))
stopifnot(all(c("Region","cluster","donor","group") %in% names(md)))

# expression matrix (sparse)
mat <- GetAssayData(merged[["RNA"]], layer = "data")

## ---- 2) Re-scale to z after cutoff (sparse-friendly, gene-wise) ----
getz <- function(v, cutoff=0.1, ids=NULL){
  keep <- is.finite(v) & v > cutoff
  z <- rep(NA_real_, length(v))
  if (sum(keep) >= 3L) z[keep] <- as.numeric(scale(v[keep]))
  if (!is.null(ids)) names(z) <- ids
  z
}

cells   <- colnames(mat)

# (a) per-gene z for the nine risk proteins
genes6A <- intersect(risk9, rownames(mat))

Z6A <- do.call(rbind, lapply(genes6A, function(g)
  getz(as.numeric(mat[g, ]), cutoff = 0.1, ids = cells)
))
rownames(Z6A) <- genes6A
colnames(Z6A) <- cells

# (b) per-gene z for receptors
gR <- unique(unlist(receptors))
genes6B <- intersect(gR, rownames(mat))

cells <- colnames(mat)
ZR <- do.call(rbind, lapply(genes6B, function(g) getz(as.numeric(mat[g,]), 0.1, ids=cells)))
rownames(ZR) <- genes6B; colnames(ZR) <- cells
## ---- 3) Long-format data (merged with cell metadata) ----
anno <- data.frame(
  cell    = colnames(mat),
  Region  = md$Region,
  cluster = md$cluster,
  donor   = md$donor,
  group   = md$group,
  stringsAsFactors = FALSE
)

df6A <- melt(Z6A, varnames=c("gene","cell"), value.name="z") %>%
  left_join(anno, by="cell") %>%
  filter(!is.na(z))

df6B <- melt(ZR, varnames=c("gene","cell"), value.name="z") %>%
  left_join(anno, by="cell") %>%
  filter(!is.na(z))

# drop combinations with too few cells or no variance
enforce_min <- function(df){
  df %>%
    group_by(Region, cluster, gene) %>%
    filter(n() >= min_cells,
           n_distinct(donor) >= min_donors,
           is.finite(var(z)) && var(z) > 0) %>%
    ungroup()
}
df6A <- enforce_min(df6A)
df6B <- enforce_min(df6B)

## ---- 4) Cell-level tests with donor-clustered robust SE (CRSE) ----
fit_cell_crse_df <- function(dat){
  tb <- table(dat$group)
  if (length(tb)<2 || any(tb==0)) return(tibble::tibble(beta=NA_real_, t=NA_real_, p=NA_real_, n=nrow(dat)))
  m <- lm(z ~ group, data=dat)
  ct <- try(lmtest::coeftest(m, vcov.=sandwich::vcovCL, cluster=~donor), silent=TRUE)
  if (inherits(ct,"try-error") || !"groupASD"%in% rownames(ct))
    return(tibble::tibble(beta=NA_real_, t=NA_real_, p=NA_real_, n=nrow(dat)))
  tibble::tibble(
    beta = unname(coef(m)["groupASD"]),
    t    = unname(ct["groupASD","t value"]),
    p    = unname(ct["groupASD","Pr(>|t|)"]),
    n    = nrow(dat)
  )
}

res6A <- df6A %>%
  dplyr::group_by(Region, cluster, gene) %>%
  dplyr::reframe(fit_cell_crse_df(dplyr::pick(everything()))) %>%  # pick() replaces the deprecated cur_data()
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(beta), !is.na(p))

## ---- 6) Risk-protein bubble plot (p<alpha; colour = effect, size = -log10 p) ----
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(scales)})

# --- map clusters to broad cell types (extend as needed) ---
broad_map <- function(cl){
  if (grepl("^AST", cl))               "Astrocyte"
  else if (grepl("^Endothelial", cl))  "Endothelial"
  else if (grepl("^(L[2-6]|Neu)", cl)) "Excitatory neuron"
  else if (grepl("^IN", cl))           "Interneuron"
  else if (grepl("^Oligo", cl))        "Oligodendrocyte"
  else if (grepl("^OPC", cl))          "OPC"
  else                                  "Other"
}

# --- shared bubble-plot function (left annotation tiles + size/colour legends) ---
bubble_with_sidebars <- function(df, title_txt, file, alpha_p=0.05, beta_cap=0.5,
                                 gene_levels=NULL, region_order=c("ACC","PFC")) {
  dd <- df %>%
    filter(!is.na(beta), !is.na(p)) %>%
    mutate(
      rowlab = paste(Region, cluster, sep=" — "),
      Region = factor(Region, levels = region_order),
      mlog10p = -log10(p),
      beta_c  = pmax(pmin(beta, beta_cap), -beta_cap),
      broad   = vapply(as.character(cluster), broad_map, character(1))
    ) %>%
    filter(p < alpha_p)

  if (nrow(dd) == 0) stop("No rows pass p < alpha_p.")

  if (is.null(gene_levels)) gene_levels <- sort(unique(dd$gene))
  dd$gene  <- factor(dd$gene,  levels = gene_levels)
  rows_ord <- dd %>% distinct(rowlab, Region) %>% arrange(Region, rowlab) %>% pull(rowlab)
  dd$rowlab <- factor(dd$rowlab, levels = rows_ord)

  # data for left annotation tiles
  annot_region <- dd %>% distinct(rowlab, Region) %>% mutate(x=0.6, key=as.character(Region))
  annot_broad  <- dd %>% distinct(rowlab, broad ) %>% mutate(x=0.2, key=as.character(broad))

  # palette
  pal_region <- scale_fill_manual("Region", values=c("ACC"="#4C78A8","PFC"="#F58518"), drop=FALSE)
  pal_broad  <- scale_fill_manual("Broad cell type", values=c(
    "Astrocyte"="#E8B9C5","Endothelial"="#8BC34A","Excitatory neuron"="#A3D0C2",
    "Interneuron"="#F5B14C","Oligodendrocyte"="#4C78A8","OPC"="#00B5AD","Other"="#B0B0B0"
  ), drop=FALSE)

  # main plot
  p <- ggplot() +
    # left annotation tiles
    geom_tile(data=annot_region, aes(x=x, y=rowlab, fill=key), width=0.35, height=0.9) +
    pal_region + ggnewscale::new_scale_fill() +
    geom_tile(data=annot_broad,  aes(x=x, y=rowlab, fill=key), width=0.35, height=0.9) +
    pal_broad  + ggnewscale::new_scale_fill() +
    # main points
    geom_point(data=dd,
               aes(x=gene, y=rowlab, size=mlog10p, fill=beta_c),
               shape=21, color="grey25", stroke=0.25) +
    scale_fill_gradient2(low="#2B6CB0", mid="white", high="#C53030",
                         limits=c(-beta_cap, beta_cap), midpoint=0,
                         name=expression(Delta~"(ASD−Ctrl)")) +
    scale_size_continuous(range=c(2.6,6.8),
                          breaks=c(-log10(0.05),3,4,5),
                          labels=c("0.05","1e−3","1e−4","1e−5"),
                          name=expression(-log[10]*"(p)")) +
    scale_x_discrete(expand=expansion(add=c(1.1, 0.4))) +
    labs(title=title_txt, x=NULL, y=NULL,
         subtitle=sprintf("Only p < %.2f shown | Size = −log10(p) | Color = Δ (cap ±%.2f)", alpha_p, beta_cap)) +
    theme_classic(base_size=10) +
    theme(
      axis.text.x = element_text(angle=60, hjust=1, vjust=1, size=9.5, color="grey15"),
      axis.text.y = element_text(size=9.5, color="grey15"),
      legend.position = "right",
      legend.title = element_text(size=10, face="bold"),
      legend.text  = element_text(size=9),
      plot.title   = element_text(face="bold", size=14),
      plot.subtitle= element_text(size=10, color="grey35", margin=margin(b=6)),
      panel.border = element_rect(fill=NA, color="black", linewidth=0.3),
      plot.margin  = margin(10, 70, 10, 10)
    )

  ggsave(file, p, width=7.6, height=6.2, device=cairo_pdf, dpi=600)
  p
}

# res6A: (Region, cluster, gene, beta, p, n)
genes6A_order <- c("AHSG","C1RL","C3","CFHR3","COL6A1","POSTN","VASN","VNN1")  # display order
p6A <- bubble_with_sidebars(
  df = res6A,
  title_txt = "Risk proteins (cell-level, CRSE p)",
  file = "riskprotein_bubble_annot.pdf",
  alpha_p = 0.05, beta_cap = 0.5,
  gene_levels = genes6A_order
)


plot6A <- res6A %>%
  mutate(mlog10p = -log10(p),
         gene = factor(gene, levels = sort(unique(gene))),
         rowlab = paste(Region, cluster, sep=" — ")) %>%
  filter(p < alpha_p)

p6A <- ggplot(plot6A,
              aes(x = gene, y = rowlab)) +
  geom_point(aes(size = mlog10p, fill = pmax(pmin(beta, 0.5), -0.5)),
             shape = 21, color = "grey20", stroke = 0.25) +
  scale_fill_gradient2(low="#2B6CB0", mid="white", high="#C53030",
                       midpoint=0, limits=c(-0.5,0.5), name=expression(Delta~"(ASD−Ctrl)")) +
  scale_size_continuous(range=c(2.8,6.8), name=expression(-log[10]*" p")) +
  labs(title="Risk proteins (cell-level, CRSE p)",
       subtitle=sprintf("Only p < %.2f shown", alpha_p),
       x=NULL, y=NULL) +
  theme_classic(base_size=10) +
  theme(axis.text.x=element_text(angle=60,hjust=1),
        legend.position="right")

ggsave("riskprotein_bubble_celllevel.pdf", p6A, width=7.5, height=5.5, device=cairo_pdf)

## ---- 7) Receptor statistics: delta and p per Region x cluster x gene ----
res6B_genes <- df6B %>%
  dplyr::group_by(Region, cluster, gene) %>%
  dplyr::reframe(fit_cell_crse_df(dplyr::pick(everything()))) %>%  # pick() replaces the deprecated cur_data()
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(beta), !is.na(p))

## ---- 8) Receptor heatmap (p<alpha) ----
# gene-to-axis labels
gene2axis <- stack(receptors) %>% rename(gene=values, axis=ind)

ggsave("receptor_readiness_heatmap.pdf", p6B, width=8.5, height=6.0, device=cairo_pdf)


# receptor results with gene-to-axis labels
gene2axis <- stack(receptors) %>% dplyr::rename(gene=values, axis=ind)
res6B_ready <- res6B_genes %>%
  dplyr::inner_join(gene2axis, by="gene") %>%
  dplyr::mutate(gene = factor(gene))  # order is set automatically or specified manually

# (optional) prefix axis name to gene labels for visual grouping
res6B_ready$gene <- factor(paste0(as.character(res6B_ready$axis), "·", as.character(res6B_ready$gene)))

p6B <- bubble_with_sidebars(
  df = res6B_ready %>% dplyr::select(Region, cluster, gene, beta, p, n),
  title_txt = "Receptor readiness (cell-level, CRSE p)",
  file = "receptor_readiness_bubble_annot.pdf",
  alpha_p = 0.05, beta_cap = 0.4
)


ggsave("receptor_readiness_heatmap.pdf", p6B, width=8.5, height=6.0, device=cairo_pdf)

library(patchwork)

# shrink left panel (0.92) and expand right panel (1.08)
p_side_tuned <- (p6A | p6B) +
  plot_layout(widths = c(0.92, 1.08), guides = "collect") &
  theme(legend.position = "right")

# panel tags
p_side_tuned <- p_side_tuned + plot_annotation(tag_prefix = "6", tag_levels = "A")

# save (keep landscape aspect ratio)
ggsave("riskprotein_readiness_side.pdf", p_side_tuned,
       width = 14, height = 6.8, device = cairo_pdf, dpi = 600)


## ---- 9) Build ligand-receptor edges (Integrin/ECM and Complement axes) ----
# sender: risk-protein genes with delta>0 and p<alpha
sender_genes <- res6A %>%
  filter(gene %in% names(lig_axis)) %>%
  mutate(axis = lig_axis[gene]) %>%
  filter(!is.na(axis), beta > 0, p < alpha_p) %>%
  transmute(Region, axis, cluster_s = cluster, ligand = gene,
            delta_s = beta, p_s = p)

# receiver: receptor genes with delta>0 and p<alpha
receiver_genes <- res6B_genes %>%
  inner_join(gene2axis, by="gene") %>%
  filter(axis %in% c("Integrin_ECM","Complement"),
         beta > 0, p < alpha_p) %>%
  transmute(Region, axis, cluster_r = cluster, receptor = gene,
            delta_r = beta, p_r = p)

# join to form edges
edges <- inner_join(sender_genes, receiver_genes,
                    by=c("Region","axis"), relationship="many-to-many") %>%
  mutate(weight   = delta_s * delta_r,
         strength = ((-log10(p_s)) + (-log10(p_r)))/2) %>%
  arrange(Region, axis, desc(weight))

# save
write.csv(edges, "ligand_receptor_edges.csv", row.names = FALSE)
suppressPackageStartupMessages({library(dplyr); library(scales); library(CellChat)})

plot_circle_cc_edges2 <- function(
  edges, ax, rg,
  topN = 25L,
  file = NULL,
  sign = c("pos","both","neg"),   # select positive / negative / both
  use_abs = FALSE,                # if TRUE, use |weight| for width and ranking
  per_node_k = NULL,              # keep at least top-k per node (limits hub dominance)
  fix_sectors = TRUE,             # fixed sector padding per node
  width = 8.0, height = 8.0,
  seed = 1234                     # reproducible tie-breaking
){
  stopifnot(is.data.frame(edges))
  sign <- match.arg(sign)

  E <- edges
  if (all(c("cluster_s","cluster_r") %in% names(E))) {
    E <- dplyr::rename(E, sender_cluster = cluster_s, receiver_cluster = cluster_r)
  }
  stopifnot(all(c("Region","axis","sender_cluster","receiver_cluster") %in% names(E)))

  if (!("weight" %in% names(E)) && all(c("delta_s","delta_r") %in% names(E))) {
    E <- dplyr::mutate(E, weight = delta_s * delta_r)
  }

  # basic filter
  E <- E %>%
    dplyr::filter(axis == ax, Region == rg, is.finite(weight)) %>%
    dplyr::mutate(
      w_use = if (use_abs) abs(weight) else weight,
      sgn   = ifelse(weight >= 0, "pos", "neg")
    )

  if (sign == "pos") E <- dplyr::filter(E, sgn == "pos")
  if (sign == "neg") E <- dplyr::filter(E, sgn == "neg")
  if (nrow(E) == 0L) { message("[skip] no edges for ", ax, "/", rg); return(invisible(NULL)) }

  # per-node top-K (optional)
  if (!is.null(per_node_k) && per_node_k > 0) {
    set.seed(seed)
    Es <- E %>% group_by(sender_cluster)  %>% slice_max(order_by = w_use, n = per_node_k, with_ties = FALSE)
    Er <- E %>% group_by(receiver_cluster)%>% slice_max(order_by = w_use, n = per_node_k, with_ties = FALSE)
    E  <- bind_rows(Es, Er) %>% distinct()
  }

  # global top-N (reproducible ties)
  set.seed(seed)
  E <- E %>%
    arrange(desc(w_use), sender_cluster, receiver_cluster) %>%
    slice_head(n = as.integer(topN))

  if (nrow(E) == 0L) { message("[skip] all filtered out"); return(invisible(NULL)) }

  # adjacency matrix
  M <- xtabs(w_use ~ sender_cluster + receiver_cluster, data = E)

  # pad to square (fixed sectors)
  labs <- sort(unique(c(rownames(M), colnames(M))))
  if (fix_sectors) {
    miss_r <- setdiff(labs, rownames(M)); if (length(miss_r)) M <- rbind(M, matrix(0, length(miss_r), ncol(M), dimnames = list(miss_r, colnames(M))))
    miss_c <- setdiff(labs, colnames(M)); if (length(miss_c)) M <- cbind(M, matrix(0, nrow(M), length(miss_c), dimnames = list(rownames(M), miss_c)))
    M <- M[labs, labs, drop = FALSE]
  } else {
    # ensure the matrix is square
    miss_r <- setdiff(labs, rownames(M)); if (length(miss_r)) M <- rbind(M, matrix(0, length(miss_r), ncol(M), dimnames = list(miss_r, colnames(M))))
    miss_c <- setdiff(labs, colnames(M)); if (length(miss_c)) M <- cbind(M, matrix(0, nrow(M), length(miss_c), dimnames = list(rownames(M), miss_c)))
    M <- M[labs, labs, drop = FALSE]
  }
  if (all(M == 0)) { message("[skip] all-zero after selection"); return(invisible(NULL)) }

  # scale
  vsize  <- scales::rescale(rowSums(M) + colSums(M), to = c(5, 14))
  nz     <- M[M > 0]
  ew_max <- sqrt(max(nz))

  if (is.null(file)) file <- sprintf("circle_%s_%s_top%d.pdf", rg, ax, topN)

  pdf(file, width = width, height = height)
  CellChat::netVisual_circle(
    net = M,
    weight.scale = FALSE,
    vertex.weight = vsize,
    edge.weight.max = ew_max,
    edge.width.max  = 8,
    alpha.edge      = 0.60,
    arrow.size      = 0.50,
    vertex.label.cex   = 0.95,
    vertex.label.color = "#1F1F1F",
    label.edge         = FALSE,
    title.name = sprintf("Axis=%s  Region=%s  Top%d (%s, %s)",
                         ax, rg, topN, sign, ifelse(use_abs,"|w|","w"))
  )
  dev.off(); message("Saved: ", file)
  invisible(list(file=file, edges_used=E))
}

plot_circle_cc_edges2(edges, "Integrin_ECM", "PFC",
  topN=200, sign="both", use_abs=TRUE, per_node_k=6, fix_sectors=TRUE, width = 6, height = 6,
  file="circle_PFC_Integrin_dense.pdf")

plot_circle_cc_edges2(edges, "Complement", "PFC",
  topN=200, sign="both", use_abs=TRUE, per_node_k=6, fix_sectors=TRUE, width = 6, height = 6,
  file="fig_circle_PFC_Complement_dense.pdf")


# ==== Ligand->Receptor coupling heatmap (dC = d_s x d_r) ====
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
})

## parameters (adjust values as needed)
alpha_p_send <- 0.10     # sender p cutoff
alpha_p_recv <- 0.10     # receiver p cutoff
val_cap      <- 0.40     # heatmap colour cap (+/-)
min_edges    <- 1        # minimum number of contributing genes per pair (within axis and region)
lab_topN     <- 12       # number of top pairs to label (by |weight|)
out_dir      <- "output/coupling_heatmap"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1) build edges: Region x axis x (sender_cluster, receiver_cluster)
build_edges_axis_region <- function(axis_target, region_target) {
  # sender: risk-axes ligands (Δ_s, p_s)
  S <- res6A %>%
    filter(gene %in% names(lig_axis)) %>%
    mutate(axis = lig_axis[gene]) %>%
    filter(axis == axis_target, Region == region_target,
           is.finite(beta), is.finite(p), p <= alpha_p_send) %>%
    transmute(Region, axis, sender_cluster = cluster,
              ligand = gene, delta_s = beta, p_s = p)

  # receiver: receptor genes (Δ_r, p_r)
  R <- res6B_genes %>%
    inner_join(gene2axis, by = "gene") %>%
    filter(axis == axis_target, Region == region_target,
           is.finite(beta), is.finite(p), p <= alpha_p_recv) %>%
    transmute(Region, axis, receiver_cluster = cluster,
              receptor = gene, delta_r = beta, p_r = p)

  if (nrow(S) == 0 || nrow(R) == 0) return(tibble())

  # join -> gene-level coupling -> aggregate per pair
  E_gene <- inner_join(S, R, by = c("Region","axis"), relationship = "many-to-many") %>%
    mutate(weight = delta_s * delta_r,
           abs_w  = abs(weight),
           # combined p (conservative): p_max = max(p_s, p_r); Fisher's method optionally available
           p_max  = pmax(p_s, p_r))

  if (nrow(E_gene) == 0) return(tibble())

  # aggregate per pair (sum multiple L-R contributions within the same sender/receiver)
  E_pair <- E_gene %>%
    group_by(Region, axis, sender_cluster, receiver_cluster) %>%
    summarise(
      n_lr   = n(),
      w_sum  = sum(weight, na.rm = TRUE),     # default: sum
      w_mean = mean(weight, na.rm = TRUE),    # reference: mean
      p_min  = suppressWarnings(min(p_max, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    filter(n_lr >= min_edges) %>%
    mutate(
      value   = w_sum,               # value shown in the heatmap (sum by default)
      value_c = pmax(pmin(value,  val_cap), -val_cap),   # colour cap
      sign    = ifelse(value >= 0, "pos", "neg"),
      mlabel  = sprintf("%s→%s", sender_cluster, receiver_cluster)
    )

  E_pair
}

# 2) draw heatmap
plot_heatmap_6C <- function(axis_target, region_target,
                            use = c("sum","mean"),
                            file = NULL, w = 7.5, h = 6.0) {
  use <- match.arg(use)
  E <- build_edges_axis_region(axis_target, region_target)
  if (nrow(E) == 0) {
    message("[6C] No edges for ", region_target, " / ", axis_target); return(invisible(NULL))
  }

  # choose the value to plot (sum or mean)
  E$val  <- if (use == "sum") E$w_sum else E$w_mean
  E$valc <- pmax(pmin(E$val, val_cap), -val_cap)

  # tile matrix (include all combinations so empty cells are visible)
  senders   <- sort(unique(E$sender_cluster))
  receivers <- sort(unique(E$receiver_cluster))
  grid <- expand.grid(sender_cluster = senders, receiver_cluster = receivers,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
    as_tibble()

  H <- grid %>%
    left_join(E %>% select(sender_cluster, receiver_cluster, val, valc, n_lr, p_min),
              by = c("sender_cluster","receiver_cluster")) %>%
    replace_na(list(val = 0, valc = 0, n_lr = 0, p_min = 1)) %>%
    mutate(sender_cluster = factor(sender_cluster, levels = senders),
           receiver_cluster = factor(receiver_cluster, levels = receivers))

  # label only the top-N pairs
  lab_df <- E %>%
    arrange(desc(abs(val))) %>%
    slice_head(n = lab_topN) %>%
    transmute(sender_cluster, receiver_cluster, lab = "*")

  p <- ggplot(H, aes(x = receiver_cluster, y = sender_cluster, fill = valc)) +
    geom_tile(color = "grey92", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2B6CB0", mid = "white", high = "#C53030",
                         limits = c(-val_cap, val_cap),
                         oob = squish, name = expression(Delta*C == Delta[s] %.% Delta[r])) +
    geom_text(data = lab_df, aes(x = receiver_cluster, y = sender_cluster, label = lab),
              color = "black", size = 3.2, fontface = 2, inherit.aes = FALSE) +
    coord_fixed() +
    labs(title = sprintf("Coupling heatmap (%s, %s)", axis_target, region_target),
         subtitle = sprintf("Value = %s of Δ_s×Δ_r per (sender, receiver); p_s≤%.2g, p_r≤%.2g; cap ±%.2f",
                            use, alpha_p_send, alpha_p_recv, val_cap),
         x = "Receiver cluster", y = "Sender cluster") +
    theme_minimal(base_size = 10.5) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30")
    )

  if (is.null(file)) {
    file <- file.path(out_dir, sprintf("coupling_heatmap_%s_%s_%s.pdf",
                                       axis_target, region_target, use))
  }
  ggsave(file, p, width = w, height = h, device = cairo_pdf, dpi = 300)
  message("Saved: ", file)

  # (optional) also keep the numeric table
  write.csv(E %>% arrange(desc(abs(val))),
            file.path(out_dir, sprintf("coupling_edges_table_%s_%s_%s.csv",
                                       axis_target, region_target, use)),
            row.names = FALSE)

  invisible(list(plot = p, table = E))
}

# 3) example runs (Integrin/Complement x PFC/ACC, using sum and mean)
axes    <- c("Integrin_ECM","Complement")
regions <- c("PFC","ACC")
for (ax in axes) for (rg in regions) {
  plot_heatmap_6C(ax, rg, use = "sum")
  plot_heatmap_6C(ax, rg, use = "mean")
}


# ==== Downstream module activation heatmap ====
suppressPackageStartupMessages({ library(ggplot2); library(scales); library(ggnewscale) })

# Short x-axis labels
module_label_short <- c(
  Integrin_downstream = "Integrin/FAK downstream",
  Complement_response = "Complement response"
)

# Prepare plotting values
res <- res %>%
  dplyr::mutate(
    rowlab = paste(Region, "—", cluster),
    valc   = pmax(pmin(beta, cap_delta), -cap_delta)
  )

# Order rows (ACC before PFC by default; adjust if needed)
row_order <- res %>%
  dplyr::distinct(Region, cluster, rowlab) %>%
  dplyr::arrange(factor(Region, levels = c("ACC","PFC")), cluster) %>%
  dplyr::pull(rowlab)

res$rowlab <- factor(res$rowlab, levels = row_order)
res$module <- factor(res$module, levels = names(module_label_short))

# Auto height based on number of rows
n_rows <- length(unique(res$rowlab))
plot_h <- max(4.8, 0.45 * n_rows + 1.6)

p <- ggplot(res, aes(x = module, y = rowlab, fill = valc)) +
  geom_tile(width = 0.92, height = 0.86, color = "white", linewidth = 0.4) +
  scale_x_discrete(labels = module_label_short, expand = expansion(add = 0.25)) +
  scale_y_discrete(limits = rev(levels(res$rowlab))) +
  scale_fill_gradient2(
    low = "#2B6CB0", mid = "white", high = "#C53030",
    limits = c(-cap_delta, cap_delta), oob = squish,
    name = expression(Delta~"(ASD−Ctrl)")
  ) +
  ggnewscale::new_scale_colour() +
  geom_point(
    data = subset(res, q < alpha_fdr),
    aes(x = module, y = rowlab),
    shape = 21, stroke = 0.6, size = 2.6, colour = "black", fill = NA
  ) +
  labs(
    title = "Downstream module activation (cell-level, donor FE + CRSE)",
    subtitle = "FDR < 0.05 marked; cap ±0.4 by default",
    x = NULL, y = NULL,
    caption = paste(
      "Downstream signals:",
      "• Integrin/FAK downstream = focal adhesion / ECM–integrin targets (e.g., PTK2/FAK, SRC, PXN, TLN1, VCL, ITGA5/ITGB1, DDR1/2, LRP1, SDC4)",
      "• Complement response = complement receptor–responsive genes (e.g., C3AR1, C5AR1, NFKBIA, TNFAIP3, JUN/FOS, SERPING1, C1QB/C1QC)",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 11.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1, margin = margin(t = 2)),
    axis.text.y = element_text(color = "#222222"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key.height = unit(3.8, "mm"),
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = 6)),
    plot.subtitle = element_text(color = "grey30", size = 10.5, margin = margin(b = 8)),
    plot.caption = element_text(color = "grey35", size = 8.8, hjust = 1, margin = margin(t = 8))
  )

ggsave("downstream_module_selected.pdf", p,
       width = 3.8, height = plot_h, device = cairo_pdf, dpi = 600)

ggplot2::ggsave(
  filename = "mr_heatmap_cis_trans.pdf",
  plot = p,
  device = cairo_pdf,   # requires system Cairo
  width = 8, height = 3, units = "in"  # figure size
)
# ===== Continuous coupling score (ranknorm(risk) x ranknorm(readiness)) =====
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(glue)
  library(lmtest); library(sandwich)
})

# --- 0) input: df_base (cell-level)
# df_base cols: cell, Region, cluster, donor, group, risk, r_int, r_cmp
stopifnot(all(c("cell","Region","cluster","donor","group","risk","r_int","r_cmp") %in% names(df_base)))

# --- 1) helper functions ---
ranknorm <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(is.finite(x)))

make_coupling <- function(df_base, axis = c("Integrin","Complement"),
                          receivers, region = c("PFC","ACC")) {
  ax <- match.arg(axis); reg <- match.arg(region)
  ax_col <- if (ax == "Integrin") "r_int" else "r_cmp"
  df_base |>
    filter(Region == reg, cluster %in% receivers) |>
    select(cell, Region, cluster, donor, group, risk, readiness = .data[[ax_col]]) |>
    filter(is.finite(risk), is.finite(readiness)) |>
    mutate(couple = ranknorm(risk) * ranknorm(readiness))
}

test_coupling <- function(cells){
  dd <- cells |>
    group_by(Region, cluster, donor, group) |>
    summarise(couple = mean(couple), .groups = "drop")

  fit_crse <- function(d){
    if (n_distinct(d$group) < 2 || n_distinct(d$donor) < 2)
      return(tibble(beta=NA_real_, p_two=NA_real_, p_one=NA_real_, n_donor=n_distinct(d$donor)))
    m  <- lm(couple ~ group, data = d)
    ct <- try(coeftest(m, vcov.=vcovCL, cluster=~donor), silent=TRUE)
    if (inherits(ct,"try-error") || !"groupASD" %in% rownames(ct))
      return(tibble(beta=NA_real_, p_two=NA_real_, p_one=NA_real_, n_donor=n_distinct(d$donor)))
    beta <- unname(coef(m)["groupASD"])
    p2   <- unname(ct["groupASD","Pr(>|t|)"])
    p1   <- if (is.finite(beta) && is.finite(p2)) if (beta >= 0) p2/2 else 1 - p2/2 else NA_real_
    tibble(beta=beta, p_two=p2, p_one=p1, n_donor=n_distinct(d$donor))
  }

  fit_wcx <- function(d){
    if (n_distinct(d$group) < 2) return(tibble(p_wcx=NA_real_, p_wcx_one=NA_real_))
    a  <- d$couple[d$group=="ASD"]; c0 <- d$couple[d$group=="Control"]
    p2 <- suppressWarnings(wilcox.test(a, c0, exact=FALSE, alternative="two.sided")$p.value)
    p1 <- suppressWarnings(wilcox.test(a, c0, exact=FALSE, alternative="greater")$p.value) # ASD > Ctrl
    tibble(p_wcx=p2, p_wcx_one=p1)
  }

  est <- dd |>
    group_by(Region, cluster) |>
    do(bind_cols(fit_crse(.), fit_wcx(.))) |>
    ungroup() |>
    mutate(q_one     = p.adjust(p_one,     method="BH"),
           q_wcx_one = p.adjust(p_wcx_one, method="BH"))
  list(donor = dd, est = est)
}

plot_coupling <- function(donor_df, est_df, axis_name, region,
                          file, width=4.8, height=3.6) {
  ee <- est_df |> filter(Region == region)

  subt <- function(cl){
    e <- ee |> filter(cluster == cl)
    if (nrow(e) == 0) return("")
    pdisp <- if (is.finite(e$p_one)) glue("p={formatC(e$p_one, format='e', digits=2)}")
            else if (is.finite(e$p_wcx_one)) glue("p={formatC(e$p_wcx_one, format='e', digits=2)}")
            else "p n/a"
    glue("{cl}: β={formatC(e$beta, format='f', digits=3)}; {pdisp}")
  }

  dd <- donor_df |> filter(Region == region)
  stopifnot(nrow(dd) > 0)

  p <- ggplot(dd, aes(group, couple, fill=group)) +
    geom_boxplot(width=.55, outlier.shape=NA, alpha=.55, colour="grey30") +
    geom_jitter(width=.08, height=0, size=1.9, alpha=.9, colour="grey15") +
    scale_fill_manual(values=c("Control"="#8da0cb","ASD"="#fc8d62"), guide="none") +
    facet_wrap(~ cluster, nrow=1, scales="free_y",
               labeller = as_labeller(\(x) sapply(x, subt))) +
    labs(title = glue("Continuous coupling ({axis_name}, {region})"),
         x = NULL, y = "Mean coupling per donor (ranknorm(risk) × ranknorm(readiness))") +
    theme_classic(base_size=11) +
    theme(strip.background = element_rect(fill="grey95", colour=NA),
          plot.title = element_text(face="bold"))
  ggsave(file, p, width=width, height=height, device=cairo_pdf, dpi=300)
  message("Saved: ", file)
  p
}

# --- 2) receiver list (with label corrections) ---
recv_integrin_all   <- c("Oligodendrocytes", "AST-PP")
recv_complement_all <- c("AST-PP","OPC")
recv_integrin   <- intersect(recv_integrin_all,   unique(df_base$cluster))
recv_complement <- intersect(recv_complement_all, unique(df_base$cluster))

# --- 3) run for PFC (repeat with region="ACC" if needed) ---
# Integrin (receiver = neurons)
cells_I_PFC <- make_coupling(df_base, "Integrin",   receivers=recv_integrin,   region="PFC")
res_I_PFC   <- test_coupling(cells_I_PFC)
write.csv(res_I_PFC$est, "continuous_coupling_PFC_Integrin.est.csv", row.names=FALSE)
plot_coupling(res_I_PFC$donor, res_I_PFC$est, "Integrin", "PFC",
              file="continuous_PFC_Integrin.pdf", width=4.8, height=3.6)

# Complement (receiver = AST-PP, OPC)
cells_C_PFC <- make_coupling(df_base, "Complement", receivers=recv_complement, region="PFC")
res_C_PFC   <- test_coupling(cells_C_PFC)
write.csv(res_C_PFC$est, "continuous_PFC_Complement.est.csv", row.names=FALSE)
plot_coupling(res_C_PFC$donor, res_C_PFC$est, "Complement", "PFC",
              file="continuous_PFC_Complement.pdf", width=4.8, height=3.6)


# --- (optional) ACC ---
# cells_I_ACC <- make_coupling(df_base, "Integrin",   recv_integrin,   "ACC")
# plot_coupling(test_coupling(cells_I_ACC)$donor, test_coupling(cells_I_ACC)$est,
#               "Integrin","ACC","continuous_ACC_Integrin.pdf", 4.8, 3.6)
# cells_C_ACC <- make_coupling(df_base, "Complement", recv_complement, "ACC")
# plot_coupling(test_coupling(cells_C_ACC)$donor, test_coupling(cells_C_ACC)$est,
#               "Complement","ACC","continuous_ACC_Complement.pdf", 4.8, 3.6)


# ====== Bubble plots + left sidebars (Region=PFC, custom broad palette) ======
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(scales); library(ggnewscale)
  library(patchwork)
})

# ---------------- Common options ----------------
alpha_p  <- 0.05
beta_cap <- 0.5
genes6A_order <- c("AHSG","C1RL","C3","CFHR3","COL6A1","POSTN","VASN","VNN1")

# -------------- Loader for semicolon CSV --------------
load_sc_csv <- function(path, numeric_cols = c("beta","p")) {
  df1 <- readr::read_delim(path, delim = ";", show_col_types = FALSE,
                           locale = readr::locale(decimal_mark = ".", grouping_mark = ","))
  need_retry <- any(numeric_cols %in% names(df1)) &&
    any(sapply(df1[numeric_cols[numeric_cols %in% names(df1)]],
               \(x) is.character(x) || all(!is.finite(suppressWarnings(as.numeric(x))))))
  if (need_retry) {
    df1 <- readr::read_delim(path, delim = ";", show_col_types = FALSE,
                             locale = readr::locale(decimal_mark = ",", grouping_mark = "."))
  }
  df1
}

# -------------- Broad cell-type mapper (7 categories) --------------

infer_broad <- function(x){
  s <- tolower(as.character(x))
  dplyr::case_when(
    grepl("^astro|^ast|astro", s)                     ~ "Astrocyte",
    grepl("endothel|bbb_endo|\\bendo\\b", s)          ~ "Endothelial",
    grepl("^ext|^neu|excit|^l[2-6]\\b", s)            ~ "Excitatory neuron",
    grepl("^int|interneuron|gaba", s)                 ~ "Interneuron",
    grepl("oligo|^odc|oligodend", s)               ~ "Oligodendrocyte",
    grepl("^opc|oligoprecursor", s)                ~ "OPC",
    grepl("^bbb|pericyte|\\bperi\\b|smooth|\\bsmc\\b", s) ~ "BBB",       # can be split into Pericyte/SMC if needed
    grepl("^mg\\d*|microgl", s)                       ~ "Microglia",
    TRUE                                              ~ "Other"
  )
}

# -------------- Palettes --------------
pal_region_vals <- c("PFC" = "#F58518")
pal_broad_vals <- c(
  # seven-colour palette
  "Astrocyte"         = "#E8B9C5",
  "Endothelial"       = "#8BC34A",
  "Excitatory neuron" = "#A3D0C2",
  "Interneuron"       = "#F5B14C",
  "Oligodendrocyte"   = "#4C78A8",
  "OPC"               = "#00B5AD",
  "BBB"        = "#C0CA33",
  "Microglia"      = "#d93316",
  "Other"      = "#B0B0B0"
  #"Immune (non-microglia)" = "#AF7AA1",
  #"Neu-mat"             = "#90A4AE",
  #"BBB (misc)"          = "#BDBDBD"
)

# -------------- Column picker --------------
`%||%` <- function(a,b) if(!is.null(a) && length(a)>0 && !is.na(a)) a else b
pick_col <- function(cols, cands){
  hit <- cands[cands %in% cols][1] %||% {
    low <- tolower(cols); names(low) <- cols
    had <- names(low)[match(tolower(cands), low)]
    had[!is.na(had)][1]
  }
  if (is.null(hit) || is.na(hit)) stop("Required column not found: ", paste(cands, collapse="/"))
  hit
}

# -------------- Main plotter: bubble + 2 sidebars --------------
bubble_with_sidebars <- function(df, title_txt, file,
                                 alpha_p = 0.05, beta_cap = 0.3,
                                 gene_levels = NULL) {

  cols <- colnames(df)
  col_region <- pick_col(cols, c("Region","region"))  # Region is already fixed to "PFC" by the caller
  col_cell   <- pick_col(cols, c("cluster","celltype","cell_type","cell_type_broad","label","Cell Type"))
  col_gene   <- pick_col(cols, c("gene","symbol","Gene","Gene Name"))
  col_delta  <- pick_col(cols, c("beta","delta","effect","coef","logFC","log2FC","ASD_minus_Ctrl","logFC (ASD vs CTL)"))
  col_p      <- pick_col(cols, c("p","pval","p_value","pvalue","Pvalue","p_adj","padj","FDR","q","qval","p-value (FDR)"))

  dat <- df %>% dplyr::transmute(
    Region = .data[[col_region]],
    cluster= .data[[col_cell]],
    gene   = as.character(.data[[col_gene]]),
    beta   = suppressWarnings(as.numeric(.data[[col_delta]])),
    p      = suppressWarnings(as.numeric(.data[[col_p]]))
  ) %>%
    mutate(
      broad   = infer_broad(cluster),
      p       = pmax(p, 1e-300),
      mlog10p = -log10(p),
      rowlab  = paste(Region, "\u2014", cluster)  # en dash
    ) %>%
    filter(is.finite(beta), is.finite(p), p < alpha_p)

  if (nrow(dat) == 0) stop("No points pass the p < alpha threshold.")

  # x order (genes)
  if (!is.null(gene_levels)) {
    dat$gene <- factor(dat$gene, levels = gene_levels)
  } else {
    dat$gene <- factor(dat$gene, levels = sort(unique(dat$gene)))
  }

  # y order: Region(PFC) → broad → cluster
  ord_df <- dat %>%
    group_by(Region, broad, cluster) %>%
    summarise(.groups="drop", m = median(beta, na.rm=TRUE)) %>%
    arrange(Region, broad, cluster) %>%
    mutate(rowlab = paste(Region, "\u2014", cluster))
  y_order <- ord_df$rowlab
  dat$rowlab <- factor(dat$rowlab, levels = y_order)

  # left stripes
  strip_region_label <- "_region"
  strip_broad_label  <- "_broad"
  x_levels <- c(strip_region_label, strip_broad_label, levels(dat$gene))

  strip_region <- dat %>% distinct(rowlab, Region) %>% mutate(x = strip_region_label)
  strip_broad  <- dat %>% distinct(rowlab, broad)  %>% mutate(x = strip_broad_label)

  p <- ggplot() +
    # (1) Region stripe
    geom_tile(data = strip_region,
              aes(x = x, y = rowlab, fill = Region),
              width = 0.9, height = 0.9) +
    scale_fill_manual(name = "Region", values = pal_region_vals, drop = FALSE) +

    ggnewscale::new_scale_fill() +
    # (2) Broad stripe
    geom_tile(data = strip_broad,
              aes(x = x, y = rowlab, fill = broad),
              width = 0.9, height = 0.9) +
    scale_fill_manual(name = "Broad cell type", values = pal_broad_vals, drop = FALSE) +

    ggnewscale::new_scale_color() +
    # (3) Bubbles (Δ color, -log10 p size)
    geom_point(data = dat,
               aes(x = gene, y = rowlab,
                   size = pmin(mlog10p, quantile(mlog10p, 0.99, na.rm=TRUE)),
                   color = pmax(pmin(beta, beta_cap), -beta_cap)),
               alpha = 0.9, shape = 16) +
    scale_color_gradient2(
      low="#2B6CB0", mid="white", high="#C53030",
      midpoint=0, limits=c(-beta_cap, beta_cap), oob=scales::squish,
      name=expression(Delta~"(ASD−Ctrl)")
    ) +
    scale_size_continuous(name=expression(-log[10]*" p"), range=c(1.6, 7.2)) +

    scale_x_discrete(limits = x_levels,
                     labels = c("","", as.character(levels(dat$gene)))) +
    scale_y_discrete(limits = rev(levels(dat$rowlab))) +
    labs(title = title_txt, x=NULL, y=NULL) +
    theme_bw(base_size = 11.5) +
    theme(
      panel.grid.major = element_line(size = 0.2, colour = "grey90"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 55, hjust = 1),
      axis.ticks = element_blank(),
      legend.box = "vertical",
      plot.margin = margin(6, 16, 6, 6)
    )

  if (!is.null(file)) {
    ggsave(file, p, width = 8.5, height = 6.0, device = cairo_pdf, dpi = 600)
  }
  return(p)
}

# ===================== Load & map your external CSVs =====================
# risk proteins - Region fixed to PFC
df6A_raw <- load_sc_csv("single_cell/33ASD/deg_results.csv")
res6A <- df6A_raw %>%
  transmute(
    Region  = "PFC",
    cluster = `Cell Type`,
    gene    = `Gene Name`,
    beta    = `logFC (ASD vs CTL)`,
    p       = `p-value (FDR)`
  )

# receptors - Region fixed to PFC
df6B_raw <- load_sc_csv("single_cell/33ASD/deg_results_Receptor.csv")
res6B <- df6B_raw %>%
  transmute(
    Region  = "PFC",
    cluster = `Cell Type`,
    gene    = `Gene Name`,
    beta    = `logFC (ASD vs CTL)`,
    p       = `p-value (FDR)`
  )

# ===================== Draw & Save =====================
p6A <- bubble_with_sidebars(
  df = res6A,
  title_txt = "Risk proteins (cell-level, CRSE p)",
  file = "riskprotein_bubble_annot.pdf",
  alpha_p = alpha_p, beta_cap = 0.2,
  gene_levels = genes6A_order
)

p6B <- bubble_with_sidebars(
  df = res6B,
  title_txt = "Receptor readiness (cell-level, CRSE p)",
  file = "receptor_readiness_bubble_annot.pdf",
  alpha_p = 0.05, beta_cap = 0.2
)

# combine side by side (narrower left panel, wider right panel)
p_side_tuned <- (p6A + theme(legend.position = "right")) |
                (p6B + theme(legend.position = "right"))

p_side_tuned <- p_side_tuned +
  plot_layout(widths = c(0.92, 1.08), guides = "collect")

ggsave("riskprotein_readiness_side.pdf", p_side_tuned,
       width = 14, height = 6.8, device = cairo_pdf, dpi = 600)

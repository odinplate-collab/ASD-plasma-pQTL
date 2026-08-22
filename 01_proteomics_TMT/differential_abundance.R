# =====================================================================
# ASD plasma proteomics: differential abundance, family-shared expression
# exclusion, and neutrality-weighted priority score.
# Self-contained (no external lab helper scripts). Input boundary:
# the normalized/batch-corrected protein expression matrix.
#
# Inputs (set ASD_ROOT or edit ROOT):
#   - normalized expression matrix (genes x samples, z-scored)   [averaged_mat]
#   - sample metadata (Sample, Family, Age, Sex)                  [Sample_info]
#   - list of significant cis/trans pQTL proteins                 [sig_pQTL]
# =====================================================================
suppressMessages({library(tidyverse); library(limma); library(scales)})
ROOT <- Sys.getenv("ASD_ROOT", ".")
OUT  <- Sys.getenv("ASD_OUT",  file.path(ROOT,"output")); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

## ---- inputs ----
load(file.path(ROOT,"TMT/250429/ASD_250701_norm_averaged_mat.r"))   # -> averaged_mat, sample.info
sinfo <- read.table(file.path(ROOT,"TMT/Sample_info_no_replicate.txt"), sep="\t", header=TRUE, quote="")
sig_pQTL <- as.vector(unlist(read.table(file.path(ROOT,"TMT/250429/ASD_sig_Cis_Trans_pQTL_protein.txt"),
                                        sep="\t", header=FALSE, quote="")))
sub_mat <- averaged_mat[rownames(averaged_mat) %in% sig_pQTL, , drop=FALSE]
cat("expression matrix:", paste(dim(averaged_mat),collapse=" x "),
    "| sig-pQTL proteins used:", nrow(sub_mat), "\n")

## ---- parameters (pre-specified) ----
thr <- 0.255            # minimum ASD deviation (mean ASD - mean baseline)
spread_thr <- 0.782     # maximum baseline spread (range across non-ASD roles)
baseline_roles <- c("Father","Mother","Sibling","non-ASD"); role_order <- c(baseline_roles,"ASD")

meta <- sinfo %>%
  transmute(Sample,
            role = if_else(Family %in% c("sibiling","Sibling"), "Sibling", Family),
            Age  = as.numeric(Age), Sex = factor(Sex, levels=c("M","F"))) %>%
  filter(role %in% role_order) %>%
  mutate(group = factor(if_else(role=="ASD","ASD","Others"), levels=c("Others","ASD")))
cat("roles present:", paste(names(table(meta$role)),table(meta$role),sep="=",collapse=", "), "\n")

## ---- (1) family-shared expression exclusion ----
expr_long <- sub_mat %>% as.data.frame() %>% rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to="Sample", values_to="z") %>% inner_join(meta, by="Sample")
role_means <- expr_long %>% group_by(gene,role) %>%
  summarise(mu=mean(z,na.rm=TRUE), .groups="drop") %>% pivot_wider(names_from=role, values_from=mu)
gene_stats <- role_means %>% rowwise() %>% mutate(
  n_base = sum(!is.na(c_across(any_of(baseline_roles)))),
  baseline_mean   = if(n_base>0) mean(c_across(any_of(baseline_roles)), na.rm=TRUE) else NA_real_,
  baseline_spread = if(n_base>1) diff(range(c_across(any_of(baseline_roles)), na.rm=TRUE)) else 0,
  d_asd = ASD - baseline_mean) %>% ungroup()
gene_cluster <- gene_stats %>% mutate(cluster = case_when(
  !is.na(ASD) & n_base>0 & baseline_spread<=spread_thr & d_asd >=  thr ~ "ASD_up_only",
  !is.na(ASD) & n_base>0 & baseline_spread<=spread_thr & d_asd <= -thr ~ "ASD_down_only",
  TRUE ~ "Mixed")) %>% select(gene, cluster)
n_up <- sum(gene_cluster$cluster=="ASD_up_only"); n_dn <- sum(gene_cluster$cluster=="ASD_down_only")
cat(sprintf(">> ASD-only UP = %d ,  ASD-only DOWN = %d\n", n_up, n_dn))

## ---- (2) differential abundance (limma: ASD vs Others + Age + Sex) ----
samples_use <- intersect(colnames(sub_mat), meta$Sample)
X <- sub_mat[, samples_use, drop=FALSE]; X <- X[rowSums(is.na(X))==0, , drop=FALSE]
meta2 <- meta %>% filter(Sample %in% samples_use) %>% arrange(match(Sample, colnames(X)))
stopifnot(identical(colnames(X), meta2$Sample))
design <- model.matrix(~ group + scale(Age) + Sex, data=meta2)
fit <- eBayes(lmFit(X, design)); cn <- colnames(design)
res_all <- tibble(gene=rownames(X),
  beta_asd=fit$coefficients[,"groupASD"], p_asd=fit$p.value[,"groupASD"],
  beta_age=fit$coefficients[,grep("^scale\\(Age\\)",cn,value=TRUE)], p_age=fit$p.value[,grep("^scale\\(Age\\)",cn,value=TRUE)],
  beta_sex=fit$coefficients[,grep("^Sex",cn,value=TRUE)], p_sex=fit$p.value[,grep("^Sex",cn,value=TRUE)]) %>%
  mutate(fdr_asd = p.adjust(p_asd, method="BH"))

## ---- (3) neutrality-weighted priority score (0-1) ----
eps <- 1e-12; w_neutral <- 0.5
rank_df <- res_all %>% mutate(
  direction    = if_else(beta_asd>0,"ASD_up","ASD_down"),
  neutrality   = (pmin(p_age,1)*pmin(p_sex,1)) / (1+abs(beta_age)+abs(beta_sex)),
  strength_raw = -log10(pmax(fdr_asd,eps)) * abs(beta_asd),
  neutrality_n = rescale(neutrality,   to=c(0,1), from=range(neutrality,   na.rm=TRUE)),
  strength_n   = rescale(strength_raw, to=c(0,1), from=range(strength_raw, na.rm=TRUE)),
  priority     = w_neutral*neutrality_n + (1-w_neutral)*strength_n)
asd_up_genes <- gene_cluster %>% filter(cluster=="ASD_up_only") %>% pull(gene)
rank_up <- rank_df %>% filter(gene %in% asd_up_genes, direction=="ASD_up") %>%
  arrange(desc(priority)) %>% mutate(rank=row_number())
cat(">> top-5 priority (ASD-up-only):\n"); print(as.data.frame(head(rank_up[,c("rank","gene","priority","fdr_asd")],5)))
cat(">> MBL2 priority rank =", which(rank_up$gene=="MBL2"), "\n")

## ---- outputs ----
write.table(gene_cluster, file.path(OUT,"family_exclusion_labels.txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(rank_df,      file.path(OUT,"priority_scores.txt"),         sep="\t", row.names=FALSE, quote=FALSE)
cat("[done] wrote family_exclusion_labels.txt, priority_scores.txt\n")

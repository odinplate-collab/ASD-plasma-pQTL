# Protein Fisher test: genotype (0/1/2) x diagnosis (ASD vs TD) Fisher's exact test for ASD pQTL variants, combined with ASD-vs-TD abundance change and pQTL direction
# Criterion (i) of the prioritization (05_prioritization/prioritization_sets.R).
# Inputs (relative to ASD_ROOT):
#   pQTL/pQTL_data/cis_pQTL_caseonly_FDR_exprcount_rm29903_20220413.tsv    (ASD cis-pQTLs, FDR < 0.05)
#   pQTL/pQTL_data/trans_pQTL_caseonly_FDR_exprcount_rm29903_20220413.tsv  (ASD trans-pQTLs, FDR < 0.05; ASD-only rows used)
#   pQTL/pQTL_data/matrixeQTL_input/ASD_pQTL_imputated_sample90_rm29903_n_alt.tsv  (alt-allele counts, 48 ASD + 42 TD)
#   TMT/250429/Normalized_expression_w_stat_90_pqtl_ttest.xlsx  (per-protein ASD/TD fold change and t-test P, 90 genotyped individuals)
# Outputs (relative to ASD_ROOT):
#   pQTL/protein_fisher_test/variant_fisher_test.tsv        (per-variant Fisher P and BH-FDR)
#   pQTL/protein_fisher_test/protein_fisher_variants.tsv    (variant-level rows passing the test, up and down)
#   pQTL/protein_fisher_test/protein_fisher_up.txt          (proteins passing, elevated in ASD; 93)
#   pQTL/protein_fisher_test/protein_fisher_down.txt        (proteins passing, decreased in ASD)
#   pQTL/protein_fisher_test/final_sig.rnk                  (up and down proteins with log2 fold change, no header)
# Usage: Rscript 05_prioritization/protein_fisher_test.R

ROOT <- Sys.getenv("ASD_ROOT", ".")

library(dplyr)
library(tidyr)
library(readxl)

out_dir <- file.path(ROOT, "pQTL", "protein_fisher_test")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## 1) ASD pQTL associations (cis: all FDR < 0.05 rows; trans: rows significant in ASD only)
cis_pQTL_case <- read.table(file.path(ROOT, "pQTL/pQTL_data/cis_pQTL_caseonly_FDR_exprcount_rm29903_20220413.tsv"),
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cis_pQTL_case$type <- "cis_case"

trans_pQTL_case <- read.table(file.path(ROOT, "pQTL/pQTL_data/trans_pQTL_caseonly_FDR_exprcount_rm29903_20220413.tsv"),
                              header = TRUE, sep = "\t", stringsAsFactors = FALSE)
trans_pQTL_case <- trans_pQTL_case[trans_pQTL_case$fdr_state == "case_only", ]
trans_pQTL_case$type <- "trans_case"

# One row per SNP-peptide association with chromosome/position parsed from the SNP id
make_pqtl_table <- function(pqtl, pval.col = "p.value", snp.col = "SNP", type.col = "type") {
  snp_split <- strsplit(pqtl[[snp.col]], split = ":")
  out <- data.frame(
    chromosome = sapply(snp_split, function(x) x[1]),
    position   = as.numeric(sapply(snp_split, function(x) x[2])),
    ref        = pqtl[["ref"]],
    alt        = pqtl[["alt"]],
    P.value    = pqtl[[pval.col]],
    FDR        = pqtl[["FDR"]],
    beta       = pqtl[["beta"]],
    type       = pqtl[[type.col]],
    gene_name  = pqtl[["gene_name"]],
    stringsAsFactors = FALSE
  )
  out$chromosome <- factor(out$chromosome, c(1:22, "X"))
  out
}

pqtl_total <- rbind(make_pqtl_table(cis_pQTL_case), make_pqtl_table(trans_pQTL_case))
pqtl_total$variantid <- paste0(pqtl_total$chromosome, ":", pqtl_total$position, ":",
                               pqtl_total$ref, ":", pqtl_total$alt)

## 2) Genotypes of the pQTL variants in the 90 unrelated individuals
geno_90 <- read.delim(file.path(ROOT, "pQTL/pQTL_data/matrixeQTL_input/ASD_pQTL_imputated_sample90_rm29903_n_alt.tsv"),
                      sep = "\t", header = TRUE)
geno_pqtl <- unique(geno_90[match(pqtl_total$variantid, geno_90$variantid, nomatch = 0), ])
rm(geno_90)

# Diagnosis from the WGS sample id prefix (controls are "SNPD_CTRL_*")
sample_info <- data.frame(sample_id = colnames(geno_pqtl)[-1],
                          group = ifelse(grepl("SNPD_CTRL", colnames(geno_pqtl)[-1]), "CTRL", "ASD"))

## 3) Per-variant Fisher's exact test on the 2 x 3 table (group x genotype 0/1/2)
genotype_long <- geno_pqtl %>%
  tidyr::pivot_longer(cols = -variantid, names_to = "sample_id", values_to = "genotype") %>%
  left_join(sample_info %>% select(sample_id, group), by = "sample_id")

variant_fisher_test <- genotype_long %>%
  group_by(variantid) %>%
  group_modify(~ {
    tab <- table(.x$group, .x$genotype)
    if (all(dim(tab) == c(2, 3))) {
      res <- suppressWarnings(fisher.test(tab))
      tibble(p.value = res$p.value)
    } else {
      tibble(p.value = NA)   # variants without all three genotype classes are not tested
    }
  }) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p.value, method = "BH"))

write.table(variant_fisher_test, file.path(out_dir, "variant_fisher_test.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## 4) Combine with ASD vs TD abundance change and pQTL direction
TMT_df <- read_excel(file.path(ROOT, "TMT/250429/Normalized_expression_w_stat_90_pqtl_ttest.xlsx"))

top_variant <- variant_fisher_test %>% dplyr::filter(p.value < 0.05)
# Each variant is linked to its first association in pqtl_total (cis before trans, ordered as in the input files)
sig_variant_df <- pqtl_total[match(top_variant$variantid, pqtl_total$variantid), ]
sig_variant_df$fisher_p   <- variant_fisher_test$p.value[match(sig_variant_df$variantid, variant_fisher_test$variantid)]
sig_variant_df$fisher_FDR <- variant_fisher_test$FDR[match(sig_variant_df$variantid, variant_fisher_test$variantid)]
sig_variant_df$tmt_fc     <- TMT_df$ASD_nonASD_fc[match(sig_variant_df$gene_name, TMT_df$Gene)]
sig_variant_df$tmt_p      <- TMT_df$ASD_nonASD_p_val[match(sig_variant_df$gene_name, TMT_df$Gene)]

# Direction concordance: protein up in ASD with a protein-raising alt allele, or down with a lowering allele
filtered_df_up <- sig_variant_df %>%
  dplyr::filter(tmt_p < 0.05 & fisher_p < 0.05) %>%
  dplyr::filter(tmt_fc > 1 & beta > 0)

filtered_df_down <- sig_variant_df %>%
  dplyr::filter(tmt_p < 0.05 & fisher_p < 0.05) %>%
  dplyr::filter(tmt_fc < 1 & beta < 0)

filtered_df <- rbind(filtered_df_up %>% mutate(direction = "up"),
                     filtered_df_down %>% mutate(direction = "down"))
write.table(filtered_df, file.path(out_dir, "protein_fisher_variants.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

candidate_up   <- unique(filtered_df_up$gene_name)
candidate_down <- unique(filtered_df_down$gene_name)

writeLines(candidate_up,   file.path(out_dir, "protein_fisher_up.txt"))
writeLines(candidate_down, file.path(out_dir, "protein_fisher_down.txt"))

sig_candidate <- c(candidate_up, candidate_down)
final_sig_df <- data.frame(gene = sig_candidate,
                           log2FC = log(sig_variant_df$tmt_fc[match(sig_candidate, sig_variant_df$gene_name)], 2))
write.table(final_sig_df, file.path(out_dir, "final_sig.rnk"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Proteins passing the protein Fisher test: up =", length(candidate_up),
    ", down =", length(candidate_down), "\n")

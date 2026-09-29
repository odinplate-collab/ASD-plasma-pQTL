# Prioritization of risk proteins: intersection of three criteria.
#   (i)   protein Fisher test, up direction              (protein_fisher_test.R; 93 proteins)
#   (ii)  elevated in ASD with an ASD-only cis-pQTL      (ASD_only_Cis_FDR.txt, DEP_up; 39 proteins)
#   (iii) ASD-only after family-shared signal exclusion  (differential_abundance.R; 52 up + 12 down)
# Writes the Venn region counts, the 12 proteins meeting all three criteria and a table of the
# DEPs carrying ASD-only cis/trans pQTLs with the three criteria (prioritization_sets.tsv).
# Inputs (relative to ASD_ROOT): pQTL/pQTL_re/ASD_only_Cis_FDR.txt, pQTL/pQTL_re/ASD_trans_FDR5.tsv
#   (02_pQTL_mapping/pqtl_mapping.R); pQTL/protein_fisher_test/protein_fisher_up.txt
#   (protein_fisher_test.R); output/family_exclusion_labels.txt (01_proteomics_TMT/differential_abundance.R)
# Outputs (relative to ASD_ROOT): output/prioritization/venn_regions.tsv, prioritized_proteins.txt,
#   prioritization_sets.tsv
# Usage: Rscript prioritization_sets.R
suppressMessages(library(data.table))
ROOT <- Sys.getenv("ASD_ROOT", ".")
OUT <- file.path(ROOT, "output", "prioritization")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cis <- fread(file.path(ROOT, "pQTL/pQTL_re/ASD_only_Cis_FDR.txt"))
tra <- fread(file.path(ROOT, "pQTL/pQTL_re/ASD_trans_FDR5.tsv"))
genes_of <- function(d, grp) sort(unique(d[DEP_group == grp & !is.na(GeneSymbol) & GeneSymbol != "", GeneSymbol]))
cis_up <- genes_of(cis, "DEP_up"); cis_dn <- genes_of(cis, "DEP_down")
tra_up <- genes_of(tra, "DEP_up"); tra_dn <- genes_of(tra, "DEP_down")

fisher <- sort(unique(readLines(file.path(ROOT, "pQTL/protein_fisher_test/protein_fisher_up.txt"))))
fam_lab <- fread(file.path(ROOT, "output/family_exclusion_labels.txt"))
family <- sort(unique(fam_lab[cluster %in% c("ASD_up_only", "ASD_down_only"), gene]))

S1 <- fisher; S2 <- family; S3 <- cis_up
regions <- data.table(
  region = c("fisher_only", "family_only", "cis_up_only", "fisher_family", "fisher_cis_up",
             "family_cis_up", "all_three"),
  n = c(length(setdiff(S1, union(S2, S3))), length(setdiff(S2, union(S1, S3))),
        length(setdiff(S3, union(S1, S2))), length(setdiff(intersect(S1, S2), S3)),
        length(setdiff(intersect(S1, S3), S2)), length(setdiff(intersect(S2, S3), S1)),
        length(Reduce(intersect, list(S1, S2, S3)))))
common <- sort(Reduce(intersect, list(S1, S2, S3)))
cat(sprintf("fisher %d | family-excluded %d | cis-up %d | all three %d\n",
            length(S1), length(S2), length(S3), length(common)))
cat("prioritized:", paste(common, collapse = ", "), "\n")
fwrite(regions, file.path(OUT, "venn_regions.tsv"), sep = "\t")
writeLines(common, file.path(OUT, "prioritized_proteins.txt"))

cols <- list(
  Cis_pQTL_up_ASD                       = cis_up,
  Cis_pQTL_down_ASD_DEP                 = cis_dn,
  Total_Cis                             = sort(union(cis_up, cis_dn)),
  Trans_pQTL_up_ASD                     = tra_up,
  Trans_pQTL_down_ASD_DEP               = tra_dn,
  Total_Trans                           = sort(union(tra_up, tra_dn)),
  `Protein fisher test`                 = fisher,
  `After removing family-shared signal` = family,
  Common                                = common)
n <- max(lengths(cols))
TAB <- as.data.table(lapply(cols, function(v) c(v, rep(NA_character_, n - length(v)))))
fwrite(TAB, file.path(OUT, "prioritization_sets.tsv"), sep = "\t", na = "")

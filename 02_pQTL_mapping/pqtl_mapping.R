# Peptide-level cis/trans pQTL mapping (Matrix eQTL) within ASD (n = 48) and within TD (n = 42), and ASD-only / shared / TD-only classification of cis-pQTLs
# Inputs (relative to ASD_ROOT):
#   TMT/Peptide/abundance_peptide_MD.tsv          (peptide log2 abundance from TMT-Integrator; annotation columns + one column per TMT channel)
#   TMT/Sample.information.txt                    (TMT sample sheet: Sample, Batch, ...)
#   TMT/Sample_info_no_replicate.txt              (per-individual sheet: Sample, Group, pQTL_90, WGS_ID)
#   reference/2025-08-20-decoys-contam-UP000005640_9606_one_protein_per_gene.fasta.fas  (FragPipe search database, UniProt human one protein per gene)
#   pQTL/pQTL_re/ASD_pQTL_imputated_case_rm29903_n_alt.tsv  (alt-allele counts, SNP x ASD sample)
#   pQTL/pQTL_re/ASD_pQTL_imputated_ctrl_rm29903_n_alt.tsv  (alt-allele counts, SNP x TD sample)
#   pQTL/pQTL_re/ASD_WGS_imputated_variantloc_rm29903.tsv   (variant positions, GRCh37)
#   pQTL/pQTL_re/pQTL_case48_cov_imputated_rm29903_20220411.txt, pQTL/pQTL_re/pQTL_ctrl42_cov_imputated_rm29903_20220411.txt
#                                                  (source of the Age and IsFemale covariate rows)
#   TMT/250429/90_DEP.txt                         (differentially abundant proteins, ASD vs TD, 90 genotyped individuals)
# Outputs (relative to ASD_ROOT):
#   pQTL/pQTL_re/pep_unique_averaged.tsv          (ComBat-corrected, unique-peptide, replicate-averaged matrix)
#   pQTL/pQTL_re/pqtl_peptide_locus.txt           (peptide genomic positions, GRCh37)
#   pQTL/pQTL_re/pep_unique_ASD.txt, pQTL/pQTL_re/pep_unique_UHC.txt   (Matrix eQTL expression input; UHC = TD)
#   pQTL/pQTL_re/pc_score_ASD.txt, pQTL/pQTL_re/pc_score_UHC.txt       (Matrix eQTL covariates: PC1-5, Age, IsFemale)
#   pQTL/pQTL_re/Results/cis_pQTL_ASD_result.tsv, trans_pQTL_ASD_result.tsv, cis_pQTL_UHC_result.tsv, trans_pQTL_UHC_result.tsv
#   pQTL/pQTL_re/ASD_only_Cis_FDR.txt             (cis-pQTLs at FDR < 0.05 in ASD and not in TD)
#   pQTL/pQTL_re/ASD_trans_FDR5.tsv               (trans-pQTLs at FDR < 0.05 in ASD, with DEP annotation)
#   pQTL/pQTL_re/pQTL_overlap_counts.tsv          (cis ASD-only / shared / TD-only counts and trans totals)
# Usage: Rscript 02_pQTL_mapping/pqtl_mapping.R

ROOT <- Sys.getenv("ASD_ROOT", ".")

library(dplyr)
library(tidyverse)
library(readr)
library(vsn)
library(DEP)
library(SummarizedExperiment)
library(S4Vectors)
library(sva)
library(Biostrings)
library(ensembldb)
library(EnsDb.Hsapiens.v75)
library(AnnotationDbi)
library(GenomeInfoDb)
library(MatrixEQTL)

base.dir   <- file.path(ROOT, "pQTL", "pQTL_re")
fasta_path <- file.path(ROOT, "reference", "2025-08-20-decoys-contam-UP000005640_9606_one_protein_per_gene.fasta.fas")
dir.create(file.path(base.dir, "Results"), recursive = TRUE, showWarnings = FALSE)

############################################################
## 1) Peptide matrix: VSN, MinProb imputation, ComBat on TMT batch
############################################################
sample.info <- read.table(file.path(ROOT, "TMT/Sample.information.txt"), sep = "\t", header = TRUE, quote = "")

peptide_df <- read.table(file.path(ROOT, "TMT/Peptide/abundance_peptide_MD.tsv"),
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Channel columns = all columns except the peptide annotation columns
anno_cols <- c("Index", "Gene", "ProteinID", "Peptide", "SequenceWindow", "Start", "End",
               "pep_length", "MaxPepProb", "ReferenceIntensity")
quant_mat <- as.matrix(peptide_df[, !(colnames(peptide_df) %in% anno_cols)])
rownames(quant_mat) <- peptide_df$Index

mat <- as.matrix(2^quant_mat)

cn <- colnames(mat)
colData <- DataFrame(sample_id = cn, row.names = cn)

# Row ids are "<UniProt>_<PEPTIDE>"
rid <- rownames(mat)
has_sep    <- grepl("_", rid)
protein_id <- ifelse(has_sep, sub("_.*$", "", rid), rid)
peptide    <- ifelse(has_sep, sub("^[^_]*_", "", rid), rid)
name_uniq  <- make.unique(peptide, sep = "_")

rowData <- DataFrame(name = name_uniq, ID = protein_id, row.names = rid)

se <- SummarizedExperiment(assays = list(intensity = mat), colData = colData, rowData = rowData)

# Keep peptides observed in at least 50% of channels
keep <- rowSums(!is.na(assay(se, "intensity"))) >= ceiling(ncol(se) * 0.50)
se   <- se[keep, ]

# VSN normalization
fit <- vsn2(assay(se, "intensity"), verbose = FALSE)
assay(se, "intensity") <- predict(fit, assay(se, "intensity"))

# Left-censored (MNAR) imputation. MinProb draws random values; the original run did not set a seed.
se_minprob <- impute(se, fun = "MinProb", q = 0.01)
norm_imp_mat <- assay(se_minprob, "intensity")

# ComBat on TMT batch (no model covariates)
tmt_batch <- as.factor(sample.info$Batch[match(colnames(norm_imp_mat), sample.info$Sample)])
combat_edata1 <- ComBat(dat = na.omit(norm_imp_mat), batch = tmt_batch, mod = NULL,
                        par.prior = TRUE, prior.plots = FALSE)

############################################################
## 2) Keep peptides unique to one protein in the search database
############################################################
peptide <- gsub("^.*_", "", rownames(combat_edata1))

# Count the database proteins containing each peptide (I/L treated as identical, case-insensitive)
check_peptide_uniqueness <- function(peptides,
                                     fasta_path,
                                     treat_IL_same = FALSE,
                                     ignore_case   = TRUE,
                                     decoy_regex   = "^(REV_|DECOY_)") {
  stopifnot(length(peptides) > 0)
  pep_raw <- peptides
  peptides <- trimws(peptides)

  aa <- readAAStringSet(fasta_path)
  ids <- names(aa)

  if (!is.null(decoy_regex) && nzchar(decoy_regex)) {
    keep <- !grepl(decoy_regex, ids)
    aa   <- aa[keep]
    ids  <- ids[keep]
  }

  if (treat_IL_same) {
    aa  <- AAStringSet(chartr("I", "L", as.character(aa)))
    peptides <- chartr("I", "L", peptides)
  }

  if (ignore_case) {
    aa  <- AAStringSet(toupper(as.character(aa)))
    peptides <- toupper(peptides)
  }

  res_list <- lapply(peptides, function(pep) {
    if (is.na(pep) || pep == "") {
      return(tibble(Peptide = pep, Match_Proteins = NA_character_, Match_Count = 0L, Is_Unique = FALSE))
    }
    cnts <- vcountPattern(pep, aa, fixed = TRUE)
    hit_idx <- which(cnts > 0L)
    hit_ids <- if (length(hit_idx) > 0) ids[hit_idx] else character(0)
    tibble(
      Peptide        = pep,
      Match_Proteins = if (length(hit_ids) == 0) NA_character_ else paste(hit_ids, collapse = ";"),
      Match_Count    = length(hit_ids),
      Is_Unique      = length(hit_ids) == 1L
    )
  })

  out <- bind_rows(res_list)
  out <- out %>% mutate(Input = pep_raw, .before = Peptide)
  out
}

ans <- check_peptide_uniqueness(peptide, fasta_path = fasta_path, treat_IL_same = TRUE, ignore_case = TRUE)

combat_edata1_filtered <- combat_edata1[ans$Is_Unique, ]

############################################################
## 3) Average technical replicates (columns "<sample>_1" / "<sample>_2")
############################################################
base_names   <- gsub("(_1|_2)$", "", colnames(combat_edata1_filtered))
grouped_cols <- split(seq_along(base_names), base_names)

averaged_mat <- sapply(grouped_cols, function(idxs) {
  if (length(idxs) == 1) {
    combat_edata1_filtered[, idxs]
  } else {
    rowMeans(combat_edata1_filtered[, idxs, drop = FALSE], na.rm = TRUE)
  }
})

write.table(averaged_mat, file = file.path(base.dir, "pep_unique_averaged.tsv"), sep = "\t", col.names = TRUE, row.names = TRUE)

averaged_mat <- as.data.frame(averaged_mat)

############################################################
## 4) ASD and TD individuals of the pQTL set (48 + 42), renamed to WGS ids
############################################################
sample.info_modified <- read.table(file.path(ROOT, "TMT/Sample_info_no_replicate.txt"), sep = "\t", header = TRUE, quote = "")

ASD_pqtl_samples <- sample.info_modified$Sample[sample.info_modified$pQTL_90 == 1 & sample.info_modified$Group == "ASD"]
UHC_pqtl_samples <- sample.info_modified$Sample[sample.info_modified$pQTL_90 == 1 & sample.info_modified$Group == "non-ASD"]

averaged_mat_ASD <- averaged_mat[, ASD_pqtl_samples[ASD_pqtl_samples %in% colnames(averaged_mat)], drop = FALSE]
colnames(averaged_mat_ASD) <- sample.info_modified$WGS_ID[match(colnames(averaged_mat_ASD), sample.info_modified$Sample)]
ASD_pep_df <- cbind(data.frame(id = rownames(averaged_mat_ASD), data.frame(averaged_mat_ASD)))

averaged_mat_UHC <- averaged_mat[, UHC_pqtl_samples[UHC_pqtl_samples %in% colnames(averaged_mat)], drop = FALSE]
colnames(averaged_mat_UHC) <- sample.info_modified$WGS_ID[match(colnames(averaged_mat_UHC), sample.info_modified$Sample)]
UHC_pep_df <- cbind(data.frame(id = rownames(averaged_mat_UHC), data.frame(averaged_mat_UHC)))

############################################################
## 5) Peptide genomic positions (Ensembl GRCh37, EnsDb v75)
############################################################
edb <- EnsDb.Hsapiens.v75

# Map "<UniProt>_<PEPTIDE>" ids to the genome through the first Ensembl protein containing the peptide
simple_pep2genome <- function(ids, edb, treat_IL_same = TRUE) {
  stopifnot(is.character(ids), inherits(edb, "EnsDb"))
  if (!any(grepl("_", ids)))
    stop("ids must be of the form 'UniProt_PEPTIDE', e.g. 'Q9Y6R7_LDSLVAQQLQSK'")

  pr <- proteins(edb, columns = c("protein_id", "protein_sequence"))
  pr <- as.data.frame(pr)
  pr <- pr[!is.na(pr$protein_sequence), ]
  aa <- AAStringSet(pr$protein_sequence); names(aa) <- pr$protein_id

  aa_str <- toupper(as.character(aa))
  if (treat_IL_same) aa_str <- chartr("I", "L", aa_str)
  aa <- AAStringSet(aa_str); names(aa) <- pr$protein_id

  peps <- sub("^[^_]*_", "", ids)
  peps <- toupper(peps)
  if (treat_IL_same) peps <- chartr("I", "L", peps)

  out_list <- lapply(seq_along(peps), function(i) {
    pep <- peps[i]
    if (!nzchar(pep)) return(data.frame(id = ids[i], chr = NA, s1 = NA_integer_, s2 = NA_integer_))

    cnts <- vcountPattern(pep, aa, fixed = TRUE)
    if (!any(cnts > 0)) return(data.frame(id = ids[i], chr = NA, s1 = NA_integer_, s2 = NA_integer_))

    j <- which(cnts > 0)[1]                      # first matching protein
    m <- matchPattern(pep, aa[[j]], fixed = TRUE)
    p_start <- start(m)[1]; p_end <- end(m)[1]   # first match position
    prot_id <- names(aa)[j]

    # proteinToGenome takes an IRanges named by Ensembl protein id and returns a list with one GRanges per range
    grl <- proteinToGenome(IRanges(start = p_start, end = p_end, names = prot_id), edb)
    if (length(grl) == 0 || length(grl[[1]]) == 0) return(data.frame(id = ids[i], chr = NA, s1 = NA_integer_, s2 = NA_integer_))

    gr <- grl[[1]]
    suppressWarnings(seqlevelsStyle(gr) <- "UCSC")
    data.frame(id = ids[i], chr = as.character(seqnames(gr))[1],
               s1 = start(gr)[1], s2 = end(gr)[1])
  })

  do.call(rbind, out_list)
}

result_df <- simple_pep2genome(ASD_pep_df$id, edb = edb, treat_IL_same = TRUE)
# Alternative haplotype / patch contigs -> primary chromosome
result_df$chr[result_df$chr == "chr14_kb021645_fix"] <- "chr14"
result_df$chr[result_df$chr == "HSCHR6_MHC_MANN"] <- "chr6"
result_df$chr[result_df$chr == "HSCHR6_MHC_QBL"] <- "chr6"
result_df$chr[result_df$chr == "HSCHR6_MHC_DBB"] <- "chr6"
result_df$chr[result_df$chr == "HSCHR6_MHC_SSTO"] <- "chr6"

result_df <- result_df[!is.na(result_df$chr), ]

ASD_pep_df <- ASD_pep_df[match(result_df$id, ASD_pep_df$id), ]
UHC_pep_df <- UHC_pep_df[match(result_df$id, UHC_pep_df$id), ]

write.table(result_df, file = file.path(base.dir, "pqtl_peptide_locus.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)

############################################################
## 6) Expression files in genotype-file sample order
############################################################
ASD_WGS   <- read_tsv(file.path(base.dir, "ASD_pQTL_imputated_case_rm29903_n_alt.tsv"), n_max = 1)
ASD_order <- setdiff(colnames(ASD_WGS), "variantid")
ASD_pep_df_ordered <- ASD_pep_df[, c("id", ASD_order)]
write.table(ASD_pep_df_ordered, file = file.path(base.dir, "pep_unique_ASD.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)

UHC_WGS   <- read_tsv(file.path(base.dir, "ASD_pQTL_imputated_ctrl_rm29903_n_alt.tsv"), n_max = 1)
UHC_order <- setdiff(colnames(UHC_WGS), "variantid")
UHC_pep_df_ordered <- UHC_pep_df[, c("id", UHC_order)]
write.table(UHC_pep_df_ordered, file = file.path(base.dir, "pep_unique_UHC.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)

############################################################
## 7) Covariates: top 5 PCs of the peptide matrix (within group) + age + sex
############################################################
# PC1-PC5 are computed from the group's peptide matrix; the Age and IsFemale rows are taken unchanged
# from the group's covariate file of the earlier mapping (pQTL_case48/ctrl42_cov_*.txt, same samples and order).
# Written as a Matrix eQTL covariate file: first column = covariate name, one column per sample.
make_covariates <- function(pep_df_ordered, demo_cov_file) {
  pep_mat <- as.matrix(pep_df_ordered[, -1])
  pca_res <- prcomp(t(pep_mat), center = TRUE, scale. = TRUE)
  pc_scores <- t(as.data.frame(pca_res$x))[1:5, ]

  demo <- read.table(demo_cov_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  demo <- as.matrix(demo[c("Age", "IsFemale"), colnames(pc_scores)])

  cov <- rbind(pc_scores, demo)
  data.frame(id = rownames(cov), cov, check.names = FALSE)
}

cov_ASD <- make_covariates(ASD_pep_df_ordered, file.path(base.dir, "pQTL_case48_cov_imputated_rm29903_20220411.txt"))
cov_UHC <- make_covariates(UHC_pep_df_ordered, file.path(base.dir, "pQTL_ctrl42_cov_imputated_rm29903_20220411.txt"))
write.table(cov_ASD, file = file.path(base.dir, "pc_score_ASD.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(cov_UHC, file = file.path(base.dir, "pc_score_UHC.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

############################################################
## 8) Matrix eQTL (linear model, cis window 1 Mb)
############################################################
run_matrix_eqtl <- function(SNP_file_name, expression_file_name, covariates_file_name,
                            output_file_name_cis, output_file_name_tra) {
  useModel = modelLINEAR
  snps_location_file_name = file.path(base.dir, "ASD_WGS_imputated_variantloc_rm29903.tsv")
  gene_location_file_name = file.path(base.dir, "pqtl_peptide_locus.txt")

  # Only associations significant at this level are saved; FDR is computed over all tests
  pvOutputThreshold_cis = 0.5
  pvOutputThreshold_tra = 1e-5

  errorCovariance = numeric()
  cisDist = 1e6

  snps = SlicedData$new()
  snps$fileDelimiter = "\t"
  snps$fileOmitCharacters = "NA"
  snps$fileSkipRows = 1
  snps$fileSkipColumns = 1
  snps$fileSliceSize = 2000
  snps$LoadFile(SNP_file_name)

  gene = SlicedData$new()
  gene$fileDelimiter = "\t"
  gene$fileOmitCharacters = "NA"
  gene$fileSkipRows = 1
  gene$fileSkipColumns = 1
  gene$fileSliceSize = 2000
  gene$LoadFile(expression_file_name)

  cvrt = SlicedData$new()
  cvrt$fileDelimiter = "\t"
  cvrt$fileOmitCharacters = "NA"
  cvrt$fileSkipRows = 1
  cvrt$fileSkipColumns = 1
  if (length(covariates_file_name) > 0) {
    cvrt$LoadFile(covariates_file_name)
  }

  snpspos = read.table(snps_location_file_name, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  snpspos$chr <- paste("chr", snpspos$chr, sep = "")
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)

  Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name_tra,
    pvOutputThreshold = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = FALSE)
}

fdr_result <- run_matrix_eqtl(
  SNP_file_name        = file.path(base.dir, "ASD_pQTL_imputated_case_rm29903_n_alt.tsv"),
  expression_file_name = file.path(base.dir, "pep_unique_ASD.txt"),
  covariates_file_name = file.path(base.dir, "pc_score_ASD.txt"),
  output_file_name_cis = file.path(base.dir, "Results", "cis_pQTL_ASD_result.tsv"),
  output_file_name_tra = file.path(base.dir, "Results", "trans_pQTL_ASD_result.tsv"))

ASD_cis   <- fdr_result$cis$eqtls
ASD_trans <- fdr_result$trans$eqtls

fdr_result_UHC <- run_matrix_eqtl(
  SNP_file_name        = file.path(base.dir, "ASD_pQTL_imputated_ctrl_rm29903_n_alt.tsv"),
  expression_file_name = file.path(base.dir, "pep_unique_UHC.txt"),
  covariates_file_name = file.path(base.dir, "pc_score_UHC.txt"),
  output_file_name_cis = file.path(base.dir, "Results", "cis_pQTL_UHC_result.tsv"),
  output_file_name_tra = file.path(base.dir, "Results", "trans_pQTL_UHC_result.tsv"))

UHC_cis   <- fdr_result_UHC$cis$eqtls
UHC_trans <- fdr_result_UHC$trans$eqtls

############################################################
## 9) FDR < 0.05 associations; ASD-only / shared / TD-only; DEP annotation
############################################################
ASD_cis_FDR5   <- ASD_cis[ASD_cis$FDR <= 0.05, ]
ASD_trans_FDR5 <- ASD_trans[ASD_trans$FDR <= 0.05, ]
UHC_cis_FDR5   <- UHC_cis[UHC_cis$FDR <= 0.05, ]
UHC_trans_FDR5 <- UHC_trans[UHC_trans$FDR <= 0.05, ]

ASD_cis_FDR5$GeneSymbol   <- peptide_df$Gene[match(ASD_cis_FDR5$gene, peptide_df$Index)]
UHC_cis_FDR5$GeneSymbol   <- peptide_df$Gene[match(UHC_cis_FDR5$gene, peptide_df$Index)]
ASD_trans_FDR5$GeneSymbol <- peptide_df$Gene[match(ASD_trans_FDR5$gene, peptide_df$Index)]
UHC_trans_FDR5$GeneSymbol <- peptide_df$Gene[match(UHC_trans_FDR5$gene, peptide_df$Index)]

# SNP-peptide pair keys
asd_key <- paste(ASD_cis_FDR5$snps, ASD_cis_FDR5$gene, sep = "|")
uhc_key <- paste(UHC_cis_FDR5$snps, UHC_cis_FDR5$gene, sep = "|")

asd_only_key <- setdiff(asd_key, uhc_key)
ASD_only_df  <- ASD_cis_FDR5[asd_key %in% asd_only_key, ]

# Differentially abundant proteins (ASD vs TD) and their direction
DEP_gene <- as.vector(unlist(read.table(file.path(ROOT, "TMT/250429/90_DEP.txt"), sep = "\t", header = FALSE)))

DEP_up_genes <- c(
  "APOL1","ABI3BP","COL6A1","C1QB","MRC2","BCHE","OMD","TNXB","CFHR5","AHSG",
  "PRG4","DPP4","LUM","GGH","IGF1","VASN","GOLM1","CD99","QSOX1","CST3",
  "IGFBP5","LBP","PI16","RARRES2","POSTN","MBL2","APOH","FETUB","MAN1A1","GSN",
  "SSC5D","CD14","CFH","F5","OLFM1","SERPINF1","VTN","CFP","ANPEP","GPLD1",
  "ALDOB","CPN2","C6","C1RL","CAT","FUCA2","C1QA","S100A8","CPN1","VNN1",
  "NCAM1","BTD","NAGLU","PTGDS","KNG1","CHGA","S100A9","IGFALS","C1S","C1R",
  "SERPIND1","APCS","ARMH4","ICAM3","PRDX2","PTPRG","COL6A3","LGALS3BP","CTBS",
  "MASP1","LAMP2","C3","COMP","ICOSLG","SELENOP","MMP2","F11","C4A","AGA",
  "APOC2","SERPINA5","ORM2","RECK","NCAM2","CD109","FAH","CFHR3","PLG","PLXDC2",
  "APOA4","ATRN","CD163","F2","AFM","CD55","LCAT","SAA4","RNASE1","ROBO4",
  "LCP1","PLA2G7","MEGF8","CLEC3B","HABP2","TGOLN2","PCYOX1","PROS1","F7","AMBP",
  "C8B","IGF2R","GNPTG","F9","PTPRJ","C2","ERP44","ADGRF5","HSPG2","OIT3",
  "C1QC","COL1A1","TNC","CA1","COL1A2","LTF","CD44","SLPI"
)

DEP_down_genes <- c(
  "PARK7","SCGB1A1","IGFBP2","ITGB1","APOD","GC","CPB2","IGLV8-61","IGHV3-72",
  "IGHM","PNP","THBS1","ENPP2","C7","AGT","APOF","APOA2","ART3","SERPINC1",
  "LRG1","SERPINA6","SERPINA3","IGHA1","PRDX6","TTR","DBH","LTBP1","MCAM",
  "JCHAIN","CP","APOA1","TALDO1","ALB","F13A1","MSN","SRGN","LDHA","SHBG",
  "SERPINB1","PKM","ADIPOQ","PPBP","TXN","FERMT3","GDI1","TPI1","ARPC4","WDR1",
  "HSPA8","CORO1A","CALR","GAPDH","ENO1","LDHB","ILK","CAVIN2","TAGLN2","PPIA",
  "ARHGDIB","PDLIM1","PLEK","MTPN","CFL1","CAP1","CALD1","VASP","VCL","ACTB",
  "TLN1","YWHAZ","PFN1","FLNA"
)

ASD_only_df$DEP_flag <- ifelse(ASD_only_df$GeneSymbol %in% DEP_gene, "o", "")
ASD_only_df$DEP_group <- case_when(
  ASD_only_df$GeneSymbol %in% DEP_up_genes   ~ "DEP_up",
  ASD_only_df$GeneSymbol %in% DEP_down_genes ~ "DEP_down",
  TRUE ~ "none"
)

ASD_trans_FDR5$DEP_flag <- ifelse(ASD_trans_FDR5$GeneSymbol %in% DEP_gene, "o", "")
ASD_trans_FDR5$DEP_group <- case_when(
  ASD_trans_FDR5$GeneSymbol %in% DEP_up_genes   ~ "DEP_up",
  ASD_trans_FDR5$GeneSymbol %in% DEP_down_genes ~ "DEP_down",
  TRUE ~ "none"
)

# ASD-only cis associations (columns: snps, gene, statistic, pvalue, FDR, beta, GeneSymbol, DEP_flag, DEP_group)
write.table(ASD_only_df, file = file.path(base.dir, "ASD_only_Cis_FDR.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(ASD_trans_FDR5, file = file.path(base.dir, "ASD_trans_FDR5.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Overlap of cis SNP-peptide pairs (FDR < 0.05) between ASD and TD; trans totals
overlap_counts <- data.frame(
  class = c("cis_ASD", "cis_TD", "cis_ASD_only", "cis_shared", "cis_TD_only", "cis_union", "trans_ASD", "trans_TD"),
  n = c(length(unique(asd_key)),
        length(unique(uhc_key)),
        length(setdiff(asd_key, uhc_key)),
        length(intersect(asd_key, uhc_key)),
        length(setdiff(uhc_key, asd_key)),
        length(union(asd_key, uhc_key)),
        nrow(ASD_trans_FDR5),
        nrow(UHC_trans_FDR5))
)
write.table(overlap_counts, file = file.path(base.dir, "pQTL_overlap_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(overlap_counts)
table(ASD_only_df$DEP_group)

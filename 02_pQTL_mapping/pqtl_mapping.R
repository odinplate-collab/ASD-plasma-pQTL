# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
# Paths are relative to the project data root.
# Set the ASD_ROOT environment variable, or run R with that directory as the working directory.
setwd(Sys.getenv("ASD_ROOT", "."))

library(MatrixEQTL)
library(tidyverse)

## Location of the package with the data files.
base.dir = '.'

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste(base.dir, "/tables/ASD_pQTL_imputated_n_alt_20220120.tsv", sep="")
snps_location_file_name = paste(base.dir, '/tables/ASD_WGS_imputated_variantloc_20220120.tsv', sep="")

# Gene expression file name
expression_file_name = paste(base.dir,'/tables/ASD_WGS_imputated_pep_expr_case_20220120.tsv', sep = '')
gene_location_file_name = paste(base.dir, "/tables/pqtl_imputated_gene_loc_20220120.tsv", sep="");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste(base.dir, "/tables/pQTL_sample91_cov_imputated_20220120.txt", sep="")

# Output file name
output_file_name_cis = paste(base.dir, "/Results/cis_pQTL_sample91_imputated_fdr_20220406.tsv", sep="")
output_file_name_tra = paste(base.dir, "/Results/trans_pQTL_sample91_imputated_fdr_20220406.tsv", sep="")

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 0.5
pvOutputThreshold_tra = 1e-5

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name)
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
snpspos$chr <- paste("chr", snpspos$chr, sep="")
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

fdr_result = Matrix_eQTL_main(
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

## Results:
View(imputated_result)
cat('Analysis done in: ', commonvar_result$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(imputated_result$cis$neqtls)
cat('Detected distant eQTLs:', '\n');
show(imputated_result$trans$neqtls)

## Plot the Q-Q plot of local and distant p-values
plot(fdr_result)

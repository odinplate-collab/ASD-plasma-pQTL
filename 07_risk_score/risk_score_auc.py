#### ASD risk score computation and AUC ####

#!/usr/bin/env python3
# 01_sample_matching.py
import pandas as pd

# file names: edit paths if needed
samplesheet_file = "samplesheet_master.csv"   # (if present) or will be created from phenotype
geno_tsv = "ASD_pQTL_imputated_sample90_rm29903_n_alt.tsv"
peptide_tsv = "peptides.tsv"   # OR peptide_longformat.tsv

# 1) load ids
geno = pd.read_csv(geno_tsv, sep='\t', nrows=5)   # just header read
geno_ids = list(geno.columns)[1:]  # first col is variantid

print(f"Genotype sample count (columns): {len(geno_ids)}")

# 2) peptide header (assumes first col=peptide id, rest are sample columns)
peptide_header = pd.read_csv(peptide_tsv, sep='\t', nrows=0)
peptide_ids = list(peptide_header.columns)[1:]
print(f"Peptide sample columns: {len(peptide_ids)}")

# 3) phenotype / master sample sheet
try:
    ss = pd.read_csv(samplesheet_file, dtype=str)
    ss_ids = ss['sample_id'].astype(str).tolist()
    print(f"Master samplesheet loaded: {len(ss_ids)} samples")
except FileNotFoundError:
    ss = None
    ss_ids = []
    print("No master samplesheet found. You can create one and rerun.")

# 4) intersection checks
set_geno = set(geno_ids)
set_pep = set(peptide_ids)
set_master = set(ss_ids)

print("Genotype ∩ Peptide:", len(set_geno & set_pep))
print("Genotype not in Peptide:", sorted(list(set_geno - set_pep))[:10])
print("Peptide not in Genotype:", sorted(list(set_pep - set_geno))[:10])

# produce id_map.csv for manual curation
idmap = pd.DataFrame({
    'geno_id': sorted(list(set_geno)),
    'present_in_peptide': [('yes' if x in set_pep else 'no') for x in sorted(list(set_geno))],
    'present_in_master': [ ('yes' if x in set_master else 'no') for x in sorted(list(set_geno))]
})
idmap.to_csv("id_map_geno_vs_pep_master.csv", index=False)
print("Wrote id_map_geno_vs_pep_master.csv - check and harmonize sample names if needed.")

# ---------------------------------------------------------------------
# Downstream steps were run from the shell using the following commands
# (peptide / pQTL / WGS out-of-fold predictions -> combined AUC report):
#
#
#
#
# python combine_and_test_auc.py \
#   --peptide results/peptide_90/peptide_9genes_90_combined_folds.csv \
#   --pqtl   results/pqtl/combined_auc_folds.csv \
#   --wgs    results/genotype_9genes/gene_scores_zscore.csv \
#   --out    results/combined_auc_summary.txt
#
#
#   python auc_report_template.py \
#   --labels metadata/labels_for_auc.csv \
#   --peptide tmp/peptide_oof.csv \
#   --pqtl    tmp/pqtl_oof.csv \
#   --wgs     tmp/wgs_oof.csv \
#   --out     results/plots/auc_report2.pdf
#

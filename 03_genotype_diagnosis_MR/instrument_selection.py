# Instrument selection for the proteins elevated in ASD with a cis-pQTL:
# lead cis variant (smallest P in the protein-level combined model, n = 90,
# diagnosis-adjusted), F = (beta/SE)^2, MAF, and retention (F > 10 and MAF >= 0.05).
# Inputs (relative to ASD_ROOT): pQTL/pQTL_re/ASD_only_Cis_FDR.txt (ASD-only cis-pQTLs,
#   written by 02_pQTL_mapping/pqtl_mapping.R),
#   derived/12_cis_summary_all_proteins/protein_level_cis_combined_n90_dxadjusted_ALL.csv.gz
#   (written by 02_pQTL_mapping/protein_level_cis.R)
# Outputs (relative to ASD_ROOT): derived/genotype_diagnosis/instrument_selection.tsv
# Usage: python instrument_selection.py
import os
import pandas as pd

ROOT = os.environ.get("ASD_ROOT", ".")
OUT = ROOT + "/derived/genotype_diagnosis"
os.makedirs(OUT, exist_ok=True)

# proteins elevated in ASD (DEP up) that carry an ASD-only cis-pQTL at FDR < 0.05
t3 = pd.read_csv(ROOT + "/pQTL/pQTL_re/ASD_only_Cis_FDR.txt", sep="\t")
cis_up = sorted(t3.loc[t3.DEP_group == "DEP_up", "GeneSymbol"].dropna().unique())
print("cis-pQTL proteins elevated in ASD:", len(cis_up))

C = pd.read_csv(ROOT + "/derived/12_cis_summary_all_proteins/protein_level_cis_combined_n90_dxadjusted_ALL.csv.gz")
C = C[C.protein.isin(cis_up)]
n_cis = C.groupby("protein").size().rename("n_cis_variants")
# lead = smallest P. For five proteins several variants in perfect LD in this sample share
# identical statistics; they are interchangeable for the genotype-diagnosis models, but
# the choice decides which variant is looked up in the outcome GWAS for Mendelian
# randomization, so the variants used are fixed here.
TIED_LEAD = {"AHSG": "3:186338382:G:C", "C1RL": "12:6760367:C:T", "C6": "5:41145722:A:G",
             "GSN": "9:124842798:C:CAA", "PTGDS": "9:138979322:G:A"}
lead = C.sort_values(["p", "chr", "pos"], kind="mergesort").groupby("protein").first()
for prot, snp in TIED_LEAD.items():
    if prot in lead.index:
        row = C[(C.protein == prot) & (C.SNP == snp)].iloc[0]
        assert row.p == lead.loc[prot, "p"], (prot, "fixed variant is not tied with the lead")
        lead.loc[prot, ["SNP", "chr", "pos", "beta", "se", "p", "maf"]] = row[["SNP", "chr", "pos", "beta", "se", "p", "maf"]].values
sel = pd.DataFrame({
    "protein": lead.index, "lead_snp": lead.SNP.values, "chr": lead.chr.values,
    "beta": lead.beta.values, "se": lead.se.values, "p": lead.p.values,
    "F": ((lead.beta / lead.se) ** 2).values, "maf": lead.maf.values,
})
sel["n_cis_variants"] = sel.protein.map(n_cis).values
sel["passes_F10"] = sel.F > 10
sel["passB"] = sel.passes_F10 & (sel.maf >= 0.05)
missing = sorted(set(cis_up) - set(sel.protein))
print("with protein-level cis statistics:", len(sel), "| without:", ", ".join(missing))
print("F > 10:", int(sel.passes_F10.sum()), "| F > 10 and MAF >= 0.05:", int(sel.passB.sum()))
sel.sort_values("protein").to_csv(OUT + "/instrument_selection.tsv", sep="\t", index=False)

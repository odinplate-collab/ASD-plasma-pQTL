# Genotype and diagnosis contributions to plasma abundance, Mendelian randomization with
# this study's instruments, and cis-effect concordance with Niu et al.
#   1. For each retained protein (F > 10, MAF >= 0.05; instrument_selection.py) and for the
#      non-selected comparison set (cis P < 1e-5, MAF >= 0.05, not elevated in ASD):
#      abundance ~ genotype + diagnosis, abundance ~ genotype x diagnosis (OLS, n = 90),
#      abundance ~ genotype within ASD and within TD, variance explained.
#   2. Wald ratio of each retained protein's lead cis variant against ten GWAS; FDR over all tests.
#   3. Per-allele cis effects in Niu et al. (P < 1e-4, pruned at r2 < 0.2 in 1000 Genomes EUR)
#      against this study's combined estimates.
# Genotype is coded as protein-raising allele count; abundance is the mean of peptide-wise
# standardized log intensities across all peptides mapping to the protein.
# Inputs (relative to ASD_ROOT):
#   derived/genotype_diagnosis/instrument_selection.tsv                     (instrument_selection.py)
#   derived/12_cis_summary_all_proteins/protein_level_cis_*_ALL.csv.gz (02_pQTL_mapping/protein_level_cis.R)
#   derived/01_samples/01_samples_covariates_n90.tsv          (sample list and diagnosis)
#   derived/02_proteomics/02b_peptide_abundance_named.tsv.gz  (peptide abundance, 90 individuals)
#   derived/11_protein_maps/measured_proteins.tsv
#   pQTL/pQTL_re/ASD_pQTL_imputated_sample90_rm29903_n_alt.tsv (genotype dosages)
#   derived/13_outcome_gwas_479/<trait>_479regions.tsv.gz     (outcome GWAS)
#   derived/14_ld_reference_479/1000G_EUR_479regions.{bed,bim,fam}; derived/15_niu_windows/Niu_*_cis1Mb.tsv.gz
# Outputs (relative to ASD_ROOT): derived/genotype_diagnosis/retained_decomposition.csv, nonselected_decomposition.csv,
#   mbl2_cells.csv, own_instrument_MR.csv, concordance_niu_all.csv
# Usage: python genotype_diagnosis_mr.py
import os, glob, warnings, time, gc
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats
warnings.filterwarnings("ignore")
ROOT = os.environ.get("ASD_ROOT", ".")
SRC = ROOT + "/derived"; OUT = ROOT + "/derived/genotype_diagnosis"; os.makedirs(OUT, exist_ok=True)
GENO = ROOT + "/pQTL/pQTL_re/ASD_pQTL_imputated_sample90_rm29903_n_alt.tsv"
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}; rc = lambda s: "".join(COMP.get(b, "N") for b in s)

cov = pd.read_csv(f"{SRC}/01_samples/01_samples_covariates_n90.tsv", sep="\t").set_index("sample_id"); samples = cov.index; asd = (cov.diagnosis == "ASD").astype(int)
pep = pd.read_csv(f"{SRC}/02_proteomics/02b_peptide_abundance_named.tsv.gz", sep="\t", index_col=0)
meas = pd.read_csv(f"{SRC}/11_protein_maps/measured_proteins.tsv", sep="\t"); g2u = dict(zip(meas.gene, meas.uniprot))
SEL = pd.read_csv(f"{OUT}/instrument_selection.tsv", sep="\t")
SEL["passB"] = (SEL.F > 10) & (SEL.maf >= 0.05); RETAINED = SEL[SEL.passB].protein.tolist()
C = pd.read_csv(f"{SRC}/12_cis_summary_all_proteins/protein_level_cis_combined_n90_dxadjusted_ALL.csv.gz")
A = pd.read_csv(f"{SRC}/12_cis_summary_all_proteins/protein_level_cis_ASD_only_ALL.csv.gz").drop_duplicates(["protein", "SNP"]).set_index(["protein", "SNP"])
T = pd.read_csv(f"{SRC}/12_cis_summary_all_proteins/protein_level_cis_TD_only_ALL.csv.gz").drop_duplicates(["protein", "SNP"]).set_index(["protein", "SNP"])
# lead cis variant = smallest P; ties among variants in perfect LD (identical statistics)
# are broken by position, which does not change any estimate
LEAD = C.sort_values(["p", "chr", "pos"], kind="mergesort").groupby("protein").first()
# for the proteins elevated in ASD, use the lead recorded in instrument_selection.tsv
Cfix = C.merge(SEL[["protein", "lead_snp"]], left_on=["protein", "SNP"], right_on=["protein", "lead_snp"]).drop(columns="lead_snp").drop_duplicates("protein").set_index("protein")
LEAD = pd.concat([Cfix, LEAD[~LEAD.index.isin(Cfix.index)]])

# non-selected comparison set: cis P < 1e-5, MAF >= 0.05, not among the proteins elevated in ASD
PEPU = pep.index.str.split("_").str[0]
up_set = set(SEL.protein)
NONSEL = LEAD[(LEAD.p < 1e-5) & (LEAD.maf >= 0.05) & (~LEAD.index.isin(up_set))].index.tolist()
NONSEL = [p for p in NONSEL if p in g2u and (PEPU == g2u[p]).sum() > 0]

# dosages (alternative-allele counts) of the lead variants, read from the genotype matrix
def read_dosages(snps):
    want, rows = set(snps), {}
    with open(GENO, encoding="utf-8") as f:
        hdr = [h.replace("-", "_") for h in f.readline().rstrip("\n").split("\t")]
        for line in f:
            k = line[:line.index("\t")]
            if k in want:
                rows[k] = [float(v) if v not in ("", "NA") else np.nan for v in line.rstrip("\n").split("\t")[1:]]
                if len(rows) == len(want): break
    return pd.DataFrame.from_dict(rows, orient="index", columns=hdr[1:])
GEN = read_dosages(set(LEAD.loc[RETAINED + NONSEL, "SNP"]) | {"10:54536839:T:G"})

def protein_level(gene):
    sub = pep.loc[PEPU == g2u[gene], samples]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0); return z.mean(axis=0), sub.shape[0]

# ---------------- per-protein decomposition (exact when dosage available, else summary-statistic approximation)
def decompose(gene):
    r = LEAD.loc[gene]; y, npep = protein_level(gene); flip = r.beta < 0
    out = dict(protein=gene, snp=r.SNP, n_pep=npep, instrument_p=r.p, F=(r.beta / r.se) ** 2, maf=r.maf)
    if r.SNP in GEN.index:
        g = GEN.loc[r.SNP, samples].astype(float); g = 2 - g if flip else g
        d = pd.DataFrame({"value": y, "genotype": g, "asd": asd})
        m0 = smf.ols("value ~ genotype + asd", d).fit(); m1 = smf.ols("value ~ genotype * asd", d).fit(); mg = smf.ols("value ~ genotype", d).fit()
        sl = {}
        for k, grp in [(1, "ASD"), (0, "TD")]:
            s = d[d.asd == k]; f = smf.ols("value ~ genotype", s).fit(); sl[grp] = (f.params.genotype, f.bse.genotype, int(s.genotype.value_counts().min()))
        out.update(exact=True, g=m0.params.genotype, g_se=m0.bse.genotype, dx=m0.params.asd, dx_se=m0.bse.asd, dx_p=m0.pvalues.asd,
                   gx=m1.params["genotype:asd"], gx_se=m1.bse["genotype:asd"], gx_p=m1.pvalues["genotype:asd"], r2_g=mg.rsquared, r2_dx=m0.rsquared - mg.rsquared,
                   b_asd=sl["ASD"][0], se_asd=sl["ASD"][1], min_asd=sl["ASD"][2], b_td=sl["TD"][0], se_td=sl["TD"][1], min_td=sl["TD"][2])
    else:
        a, t = A.loc[(gene, r.SNP)], T.loc[(gene, r.SNP)]; ya, yt = y[asd == 1], y[asd == 0]
        raw = ya.mean() - yt.mean(); raw_se = np.sqrt(ya.var(ddof=1) / len(ya) + yt.var(ddof=1) / len(yt)); dg = 2 * (a.eaf - t.eaf)
        adj = raw - r.beta * dg; adj_se = np.sqrt(raw_se ** 2 + (dg * r.se) ** 2); s = -1 if flip else 1
        gx = s * (a.beta - t.beta); gx_se = np.sqrt(a.se ** 2 + t.se ** 2); tt = (r.beta / r.se) ** 2
        out.update(exact=False, g=s * r.beta, g_se=r.se, dx=adj, dx_se=adj_se, dx_p=2 * stats.norm.sf(abs(adj / adj_se)), gx=gx, gx_se=gx_se, gx_p=2 * stats.norm.sf(abs(gx / gx_se)),
                   r2_g=tt / (tt + 88), r2_dx=adj ** 2 * asd.var() / y.var(), b_asd=s * a.beta, se_asd=a.se, min_asd=np.nan, b_td=s * t.beta, se_td=t.se, min_td=np.nan)
    return out

ST = pd.DataFrame([decompose(p) for p in RETAINED]).set_index("protein")
ST["slope_diff"] = ST.b_asd - ST.b_td; ST["slope_diff_se"] = np.sqrt(ST.se_asd ** 2 + ST.se_td ** 2); ST["slope_diff_p"] = 2 * stats.norm.sf(abs(ST.slope_diff / ST.slope_diff_se))
ST.to_csv(f"{OUT}/retained_decomposition.csv"); print("retained set:", len(ST), "| exact", int(ST.exact.sum()))

NS = pd.DataFrame([decompose(p) for p in NONSEL]).set_index("protein")
NS["slope_diff"] = NS.b_asd - NS.b_td; NS["slope_diff_se"] = np.sqrt(NS.se_asd ** 2 + NS.se_td ** 2); NS["slope_diff_p"] = 2 * stats.norm.sf(abs(NS.slope_diff / NS.slope_diff_se))
NS.to_csv(f"{OUT}/nonselected_decomposition.csv"); print("non-selected set:", len(NS))

# MBL2 per-individual values (rs7899547), used for the within-genotype comparisons
y, _ = protein_level("MBL2"); g = GEN.loc["10:54536839:T:G", samples].astype(float)
pd.DataFrame({"value": y, "genotype": g, "group": cov.diagnosis}).to_csv(f"{OUT}/mbl2_cells.csv")

# ---------------- Mendelian randomization with this study's instruments (Wald ratio)
def load_outcome(t):
    d = pd.read_csv(f"{SRC}/13_outcome_gwas_479/{t}_479regions.tsv.gz", sep="\t", low_memory=False); c = {k.lower(): k for k in d.columns}
    def col(*names):
        for n in names:
            if n.lower() in c: return d[c[n.lower()]]
    o = pd.DataFrame({"chr": pd.to_numeric(col("CHR", "#CHROM", "CHROM", "chromosome").astype(str).str.replace("chr", ""), errors="coerce"), "pos": pd.to_numeric(col("BP", "POS", "base_pair_location"), errors="coerce"),
                      "ea": col("A1", "EA", "effect_allele").astype(str).str.upper(), "oa": col("A2", "NEA", "other_allele").astype(str).str.upper()})
    beta = col("BETA", "Beta", "stdBeta", "beta"); o["beta"] = beta.astype(float) if beta is not None else np.log(col("OR").astype(float))
    o["se"] = col("SE", "standard_error").astype(float); o["p"] = col("P", "PVAL", "Pval", "p_value").astype(float)
    eaf = col("FCON", "FRQ_U_186843", "EAF", "EAF_HRC", "effect_allele_frequency"); o["eaf"] = eaf.astype(float) if eaf is not None else np.nan
    o = o.dropna(subset=["chr", "pos", "beta", "se"]); o["key"] = o.chr.astype(int).astype(str) + ":" + o.pos.astype(int).astype(str)
    return o.drop_duplicates("key").set_index("key")[["ea", "oa", "beta", "se", "p", "eaf"]]
TRAITS = ["ASD", "ADHD", "SCZ_EAS", "SCZ_EUR", "BIP", "MDD", "XDIS", "AD", "IQ", "EDU"]
def harm1(ea, oa, fx, o):
    """Allele harmonization; palindromic variants with MAF > 0.42 are dropped."""
    if o is None: return np.nan
    flip = 1 if (ea, oa) == (o.ea, o.oa) else -1 if (ea, oa) == (o.oa, o.ea) else 1 if (rc(ea), rc(oa)) == (o.ea, o.oa) else -1 if (rc(ea), rc(oa)) == (o.oa, o.ea) else np.nan
    if np.isnan(flip): return np.nan
    if (ea, oa) in {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}:
        if min(fx, 1 - fx) > 0.42: return np.nan
        if pd.notna(o.eaf):
            fy = o.eaf if flip == 1 else 1 - o.eaf
            if (fx < 0.5) != (fy < 0.5): flip = -flip
    return flip
KEYS = {p: f"{int(LEAD.loc[p].chr)}:{int(LEAD.loc[p].pos)}" for p in RETAINED}
rows = []; t0 = time.time()
for t in TRAITS:
    O = load_outcome(t); O = O[O.index.isin(set(KEYS.values()))].copy(); print(t, len(O), f"{time.time()-t0:.0f}s", flush=True)
    for p in RETAINED:
        r = LEAD.loc[p]; key = KEYS[p]; ea, oa = r.effect_allele, r.other_allele; fx = r.eaf if "eaf" in r else np.nan
        o = O.loc[key] if key in O.index else None
        fl = harm1(ea, oa, fx, o) if o is not None else np.nan
        if o is None or np.isnan(fl): rows.append(dict(protein=p, trait=t, b=np.nan, se=np.nan, p=np.nan)); continue
        by = fl * o.beta; b = by / r.beta; se = o.se / abs(r.beta); rows.append(dict(protein=p, trait=t, b=b, se=se, p=2 * stats.norm.sf(abs(b / se))))
    del O; gc.collect()
OWN = pd.DataFrame(rows); ok = OWN.p.notna(); OWN.loc[ok, "fdr"] = stats.false_discovery_control(OWN.p[ok].values); OWN.to_csv(f"{OUT}/own_instrument_MR.csv", index=False)
print("own-instrument MR tests:", int(ok.sum()), " FDR<0.05:", int((OWN.fdr < 0.05).sum()))

# ---------------- cis-effect concordance with Niu et al. across all instrumented loci
from bed_reader import open_bed
bed = open_bed(f"{SRC}/14_ld_reference_479/1000G_EUR_479regions.bed")
bim = pd.read_csv(f"{SRC}/14_ld_reference_479/1000G_EUR_479regions.bim", sep="\t", header=None, names=["chr", "rs", "cm", "pos", "a1", "a2"])
bim["key"] = bim.chr.astype(str) + ":" + bim.pos.astype(str); bim["idx"] = np.arange(len(bim)); BIM = bim.drop_duplicates("key").set_index("key")
def r2_with(lead_key, keys):
    keys = [k for k in keys if k in BIM.index]
    if lead_key not in BIM.index or not keys: return pd.Series(dtype=float)
    idx = [int(BIM.loc[lead_key, "idx"])] + [int(BIM.loc[k, "idx"]) for k in keys]
    G = bed.read(index=np.s_[:, idx], dtype="float32"); G = np.where(np.isnan(G), np.nanmean(G, axis=0), G); G = G - G.mean(axis=0)
    num = (G[:, 1:] * G[:, [0]]).sum(axis=0); den = np.sqrt((G[:, 1:] ** 2).sum(axis=0) * (G[:, 0] ** 2).sum()); return pd.Series((num / den) ** 2, index=keys)
own_all = C.assign(key=C.chr.astype(str) + ":" + C.pos.astype(str))
rows = []; t0 = time.time()
for f in sorted(glob.glob(f"{SRC}/15_niu_windows/Niu_*_cis1Mb.tsv.gz")):
    prot = os.path.basename(f).split("_")[1]
    if prot not in set(own_all.protein): continue
    d = pd.read_csv(f, sep="\t", usecols=["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "p_value"])
    d = d[d.p_value < 1e-4]
    if len(d) == 0: continue
    d["key"] = d.chromosome.astype(str) + ":" + d.base_pair_location.astype(str); d = d.sort_values("p_value").drop_duplicates("key").set_index("key")
    keep = []
    for k in d.index[:400]:
        if k not in BIM.index: continue
        if not keep or (r2_with(k, keep) < 0.2).all(): keep.append(k)
    o = own_all[own_all.protein == prot].drop_duplicates("key").set_index("key")
    for k in keep:
        if k not in o.index: continue
        e, w = d.loc[k], o.loc[k]; ea, oa = str(e.effect_allele).upper(), str(e.other_allele).upper()
        s = 1 if (ea, oa) == (w.effect_allele, w.other_allele) else -1 if (ea, oa) == (w.other_allele, w.effect_allele) else 1 if (rc(ea), rc(oa)) == (w.effect_allele, w.other_allele) else -1 if (rc(ea), rc(oa)) == (w.other_allele, w.effect_allele) else np.nan
        if np.isnan(s): continue
        rows.append(dict(protein=prot, key=k, beta_niu=e.beta, p_niu=e.p_value, beta_own=s * w.beta, se_own=w.se, p_own=w.p))
CONC = pd.DataFrame(rows); CONC.to_csv(f"{OUT}/concordance_niu_all.csv", index=False)
print("concordance points:", len(CONC), "loci:", CONC.protein.nunique(), f"{time.time()-t0:.0f}s")

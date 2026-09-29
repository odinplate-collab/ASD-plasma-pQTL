# Proteome-wide two-sample Mendelian randomization and colocalization with external instruments.
# Exposures: Niu et al. 2025 plasma pQTL (n = 1,909) for every measured protein with a cis
#   instrument, and UKB-PPP for six proteins (reference). Instruments: cis window (gene +/- 1 Mb),
#   P < 5e-8, greedy LD clumping at r2 < 0.001 in 1000 Genomes EUR.
# Estimators: Wald ratio (one instrument), inverse-variance weighted with multiplicative random
#   effects (two or more), weighted median and MR-Egger (three or more). FDR across all Niu tests
#   and within each outcome. coloc.abf (p1 = p2 = 1e-4, p12 = 1e-5) on lead variant +/- 250 kb.
# Outcomes: ten GWAS, cut to the measured-protein windows (GRCh37).
# Inputs (relative to ASD_ROOT): derived/13_outcome_gwas_479/, derived/14_ld_reference_479/,
#   derived/11_protein_maps/measured_proteins_regions_1Mb.tsv, derived/15_niu_windows/,
#   derived/08_external_pqtl/ukbppp/
# Outputs (relative to ASD_ROOT): derived/proteome_wide_MR/pw_MR_coloc.csv, pw_instruments.csv, pw_coverage.csv
# Usage: python niu_mr_colocalization.py
import os, sys, glob, warnings, time
import numpy as np, pandas as pd
from scipy import stats
from bed_reader import open_bed
warnings.filterwarnings("ignore")
ROOT = os.environ.get("ASD_ROOT", ".")
SRC = ROOT + "/derived"; OUT = ROOT + "/derived/proteome_wide_MR"; os.makedirs(OUT, exist_ok=True)
P_INST, R2_CLUMP, COLOC_WIN = 5e-8, 0.001, 250_000
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
def rc(s): return "".join(COMP.get(b, "N") for b in s)

TRAITS = {"ASD": ("ASD", "cc"), "ADHD": ("ADHD", "cc"), "SCZ_EAS": ("SCZ (EAS)", "cc"), "SCZ_EUR": ("SCZ (EUR)", "cc"), "BIP": ("Bipolar", "cc"),
          "MDD": ("Depression", "cc"), "XDIS": ("Cross-disorder", "quant"), "AD": ("Alzheimer's", "cc"), "IQ": ("IQ", "quant"), "EDU": ("Education", "quant")}

def load_outcome(t):
    d = pd.read_csv(f"{SRC}/13_outcome_gwas_479/{t}_479regions.tsv.gz", sep="\t", low_memory=False)
    c = {k.lower(): k for k in d.columns}
    def col(*names):
        for n in names:
            if n.lower() in c: return d[c[n.lower()]]
        return None
    o = pd.DataFrame({"chr": pd.to_numeric(col("CHR", "#CHROM", "CHROM", "chromosome").astype(str).str.replace("chr", ""), errors="coerce"),
                      "pos": pd.to_numeric(col("BP", "POS", "base_pair_location"), errors="coerce"),
                      "ea": col("A1", "EA", "effect_allele").astype(str).str.upper(), "oa": col("A2", "NEA", "other_allele").astype(str).str.upper()})
    beta = col("BETA", "Beta", "stdBeta", "beta"); o["beta"] = beta.astype(float) if beta is not None else np.log(col("OR").astype(float))
    o["se"] = col("SE", "standard_error").astype(float); o["p"] = col("P", "PVAL", "Pval", "p_value").astype(float)
    eaf = col("FCON", "FRQ_U_186843", "EAF", "EAF_HRC", "effect_allele_frequency"); o["eaf"] = eaf.astype(float) if eaf is not None else np.nan
    o = o.dropna(subset=["chr", "pos", "beta", "se"]); o["chr"] = o.chr.astype(int); o["pos"] = o.pos.astype(int)
    o["key"] = o.chr.astype(str) + ":" + o.pos.astype(str)
    return o.drop_duplicates("key").set_index("key")[["ea", "oa", "beta", "se", "p", "eaf"]]
t0 = time.time(); OUTC = {t: load_outcome(t) for t in TRAITS}; print("outcomes loaded", {t: len(v) for t, v in OUTC.items()}, f"{time.time()-t0:.0f}s", flush=True)

bed = open_bed(f"{SRC}/14_ld_reference_479/1000G_EUR_479regions.bed")
bim = pd.read_csv(f"{SRC}/14_ld_reference_479/1000G_EUR_479regions.bim", sep="\t", header=None, names=["chr", "rs", "cm", "pos", "a1", "a2"])
bim["key"] = bim.chr.astype(str) + ":" + bim.pos.astype(str); bim["idx"] = np.arange(len(bim)); BIM = bim.drop_duplicates("key").set_index("key")
def ld_r2(keys):
    have = [k for k in keys if k in BIM.index]; idx = [int(BIM.loc[k, "idx"]) for k in have]
    if len(idx) == 0: return pd.DataFrame()
    G = bed.read(index=np.s_[:, idx], dtype="float32"); G = np.where(np.isnan(G), np.nanmean(G, axis=0), G)
    R = np.corrcoef(G.T) ** 2 if len(idx) > 1 else np.ones((1, 1)); return pd.DataFrame(R, index=have, columns=have)
def clump(e):
    cand = e[e.p < P_INST].sort_values("p"); cand = cand[cand.index.isin(BIM.index)]
    if len(cand) == 0: return cand
    if len(cand) > 3000: cand = cand.iloc[:3000]
    R = ld_r2(cand.index.tolist()); keep = []
    for k in cand.index:
        if all(R.loc[k, j] < R2_CLUMP for j in keep): keep.append(k)
    return cand.loc[keep]

REG = pd.read_csv(f"{SRC}/11_protein_maps/measured_proteins_regions_1Mb.tsv", sep="\t").set_index("protein")
def load_niu(f):
    d = pd.read_csv(f, sep="\t", usecols=["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value"])
    e = pd.DataFrame({"chr": d.chromosome.astype(int), "pos": d.base_pair_location.astype(int), "ea": d.effect_allele.str.upper(), "oa": d.other_allele.str.upper(),
                      "beta": d.beta, "se": d.standard_error, "p": d.p_value, "eaf": d.effect_allele_frequency})
    e = e[e.ea.str.match("^[ACGT]+$") & e.oa.str.match("^[ACGT]+$")]; e["key"] = e.chr.astype(str) + ":" + e.pos.astype(str)
    return e.sort_values("p").drop_duplicates("key").set_index("key")
def load_ukb(p):
    d = pd.read_csv(f"{SRC}/08_external_pqtl/ukbppp/UKBPPP_{p}_cis1Mb.tsv.gz", sep="\t")
    e = pd.DataFrame({"chr": d.chr.astype(int), "pos": d.pos_grch37.astype(int), "ea": d.effect_allele.str.upper(), "oa": d.other_allele.str.upper(), "beta": d.beta, "se": d.se, "p": d.p, "eaf": d.eaf})
    e = e[e.ea.str.match("^[ACGT]+$") & e.oa.str.match("^[ACGT]+$")]; e["key"] = e.chr.astype(str) + ":" + e.pos.astype(str)
    return e.sort_values("p").drop_duplicates("key").set_index("key")

def harmonise(ex, out):
    m = ex.join(out, how="inner", lsuffix="_x", rsuffix="_y")
    if len(m) == 0: return m
    ea, oa, yea, yoa = m.ea_x.values, m.oa_x.values, m.ea_y.values, m.oa_y.values
    cea = np.array([rc(s) for s in ea]); coa = np.array([rc(s) for s in oa])
    flip = np.full(len(m), np.nan)
    flip[(ea == yea) & (oa == yoa)] = 1; flip[(ea == yoa) & (oa == yea)] = -1
    un = np.isnan(flip); flip[un & (cea == yea) & (coa == yoa)] = 1; flip[un & (cea == yoa) & (coa == yea)] = -1
    pal = np.array([(a, b) in {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")} for a, b in zip(ea, oa)])
    fx, fy = m.eaf_x.values, m.eaf_y.values
    amb = pal & (np.minimum(fx, 1 - fx) > 0.42); flip[amb] = np.nan
    both = pal & ~amb & ~np.isnan(fy) & ~np.isnan(fx)
    fy_al = np.where(flip == 1, fy, 1 - fy); sw = both & ((fx < 0.5) != (fy_al < 0.5)); flip[sw] = -flip[sw]
    m = m.assign(flip=flip).dropna(subset=["flip"])
    return pd.DataFrame({"bx": m.beta_x, "sx": m.se_x, "px": m.p_x, "fx": m.eaf_x, "by": m.flip * m.beta_y, "sy": m.se_y, "py": m.p_y}, index=m.index)

def mr(h):
    k = len(h); res = dict(n_snp=k)
    if k == 0: return res
    bx, by, sy = h.bx.values, h.by.values, h.sy.values
    if k == 1: res.update(b=by[0] / bx[0], se=sy[0] / abs(bx[0]), method="Wald")
    else:
        w = 1 / (sy / bx) ** 2; ratio = by / bx; b = (w * ratio).sum() / w.sum(); se = np.sqrt(1 / w.sum())
        Q = (w * (ratio - b) ** 2).sum(); phi = max(1, Q / (k - 1)); res.update(b=b, se=se * np.sqrt(phi), method="IVW-MRE", Q=Q, Q_p=stats.chi2.sf(Q, k - 1))
        if k >= 3:
            o = np.argsort(ratio); r_s, w_s = ratio[o], w[o]; cw = np.cumsum(w_s) / w_s.sum(); res["b_wmed"] = np.interp(0.5, cw - w_s / w_s.sum() / 2, r_s)
            X = np.column_stack([np.ones(k), bx]); W = np.diag(1 / sy ** 2); bhat = np.linalg.solve(X.T @ W @ X, X.T @ W @ by); resid = by - X @ bhat
            s2 = max(1, (resid ** 2 / sy ** 2).sum() / (k - 2)); cov = s2 * np.linalg.inv(X.T @ W @ X)
            res.update(b_egger=bhat[1], se_egger=np.sqrt(cov[1, 1]), egger_int=bhat[0], egger_int_p=2 * stats.norm.sf(abs(bhat[0] / np.sqrt(cov[0, 0]))))
    res["p"] = 2 * stats.norm.sf(abs(res["b"] / res["se"])); return res

def logsum(x): m = np.max(x); return m + np.log(np.sum(np.exp(x - m)))
def abf(beta, varb, sd_prior): r = sd_prior ** 2 / (sd_prior ** 2 + varb); return 0.5 * (np.log(1 - r) + r * beta ** 2 / varb)
def coloc(h, ttype, p1=1e-4, p2=1e-4, p12=1e-5):
    l1 = abf(h.bx.values, h.sx.values ** 2, 0.15); l2 = abf(h.by.values, h.sy.values ** 2, 0.2 if ttype == "cc" else 0.15)
    lH1 = np.log(p1) + logsum(l1); lH2 = np.log(p2) + logsum(l2); lH4 = np.log(p12) + logsum(l1 + l2)
    lH3 = np.log(p1) + np.log(p2) + (logsum(np.add.outer(l1, l2)[~np.eye(len(l1), dtype=bool)]) if len(l1) > 1 else -np.inf)
    ls = np.array([0.0, lH1, lH2, lH3, lH4]); pp = np.exp(ls - logsum(ls))
    return dict(PP_H0=pp[0], PP_H1=pp[1], PP_H2=pp[2], PP_H3=pp[3], PP_H4=pp[4], n_coloc=len(h))

rows, inst_rows, cov_rows = [], [], []
sources = [("Niu", os.path.basename(f).split("_")[1], f) for f in sorted(glob.glob(f"{SRC}/15_niu_windows/Niu_*_cis1Mb.tsv.gz"))]
sources += [("UKB-PPP", p, None) for p in ["AHSG", "C1RL", "C3", "CFH", "F11", "MBL2"]]
t0 = time.time()
for i, (src, prot, f) in enumerate(sources):
    try:
        e = load_niu(f) if src == "Niu" else load_ukb(prot)
    except Exception as ex:
        print("skip", src, prot, ex); continue
    ins = clump(e); n_gw = int((e.p < P_INST).sum())
    cov_rows.append(dict(source=src, protein=prot, n_cis_variants=len(e), n_gw_sig=n_gw, n_instruments=len(ins), min_p=e.p.min()))
    if len(ins) == 0: continue
    for k, r in ins.iterrows():
        inst_rows.append(dict(source=src, protein=prot, key=k, rs=BIM.loc[k, "rs"], ea=r.ea, oa=r.oa, beta=r.beta, se=r.se, p=r.p, eaf=r.eaf, F=(r.beta / r.se) ** 2))
    lead = e.p.idxmin(); lpos = e.loc[lead, "pos"]; reg = e[e.pos.between(lpos - COLOC_WIN, lpos + COLOC_WIN)]
    for t, (lab, ttype) in TRAITS.items():
        h = harmonise(ins, OUTC[t]); m = mr(h); hc = harmonise(reg, OUTC[t]); c = coloc(hc, ttype) if len(hc) > 50 else {}
        rows.append(dict(source=src, protein=prot, trait=t, label=lab, lead=lead, min_p_gwas_window=hc.py.min() if len(hc) else np.nan, **m, **c))
    if i % 25 == 0: print(f"{i}/{len(sources)} {prot} instruments={len(ins)} {time.time()-t0:.0f}s", flush=True)
R = pd.DataFrame(rows); I = pd.DataFrame(inst_rows); CV = pd.DataFrame(cov_rows)
for src in ["Niu", "UKB-PPP"]:
    sub = R[(R.source == src) & R.p.notna()]
    if len(sub): R.loc[sub.index, "fdr"] = stats.false_discovery_control(sub.p.values)
    for t in TRAITS:
        s2 = sub[sub.trait == t]
        if len(s2): R.loc[s2.index, "fdr_within_trait"] = stats.false_discovery_control(s2.p.values)
R.to_csv(f"{OUT}/pw_MR_coloc.csv", index=False); I.to_csv(f"{OUT}/pw_instruments.csv", index=False); CV.to_csv(f"{OUT}/pw_coverage.csv", index=False)
print("done", len(R), "tests;", CV.groupby("source").n_instruments.apply(lambda s: (s > 0).sum()).to_dict(), "proteins with instruments")

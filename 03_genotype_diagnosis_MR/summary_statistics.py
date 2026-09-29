# Summary statistics of the genotype and diagnosis analyses:
#   - DerSimonian-Laird random-effects pooling of the genotype, diagnosis and interaction
#     estimates across the retained and the non-selected proteins
#   - variance explained by genotype and diagnosis (medians; MBL2)
#   - agreement of genotype effects estimated within ASD and within TD (fixed-effect pooled
#     difference, correlation)
#   - MBL2: ASD minus TD difference within each rs7899547 genotype (Welch t-test)
#   - Mendelian randomization summaries (this study's and Niu et al. instruments)
#   - cis-effect concordance with Niu et al. by effect-size bin
# Inputs (relative to ASD_ROOT): derived/genotype_diagnosis/instrument_selection.tsv, retained_decomposition.csv,
#   nonselected_decomposition.csv, mbl2_cells.csv, own_instrument_MR.csv, concordance_niu_all.csv
#   (genotype_diagnosis_mr.py); derived/proteome_wide_MR/pw_MR_coloc.csv (niu_mr_colocalization.py)
# Outputs (relative to ASD_ROOT): derived/genotype_diagnosis/summary_statistics.txt, derived/genotype_diagnosis/mbl2_within_genotype.csv
# Usage: python summary_statistics.py
import os
import numpy as np, pandas as pd
from scipy import stats

ROOT = os.environ.get("ASD_ROOT", ".")
DAT = ROOT + "/derived/genotype_diagnosis"


def pool_re(b, se):
    """DerSimonian-Laird random-effects pooled estimate with I2."""
    b, se = np.asarray(b, float), np.asarray(se, float); w = 1 / se ** 2; bf = (w * b).sum() / w.sum(); Q = (w * (b - bf) ** 2).sum(); df = len(b) - 1
    tau2 = max(0, (Q - df) / (w.sum() - (w ** 2).sum() / w.sum())); wr = 1 / (se ** 2 + tau2); br = (wr * b).sum() / wr.sum(); ser = np.sqrt(1 / wr.sum())
    return dict(b=br, se=ser, lo=br - 1.96 * ser, hi=br + 1.96 * ser, p=2 * stats.norm.sf(abs(br / ser)), I2=max(0, (Q - df) / Q) * 100 if Q > 0 else 0)


SEL = pd.read_csv(f"{DAT}/instrument_selection.tsv", sep="\t"); SEL["passB"] = (SEL.F > 10) & (SEL.maf >= 0.05)
ST = pd.read_csv(f"{DAT}/retained_decomposition.csv", index_col=0); NS = pd.read_csv(f"{DAT}/nonselected_decomposition.csv", index_col=0)
MB = pd.read_csv(f"{DAT}/mbl2_cells.csv", index_col=0); OWN = pd.read_csv(f"{DAT}/own_instrument_MR.csv")
CONC = pd.read_csv(f"{DAT}/concordance_niu_all.csv")
PW = pd.read_csv(ROOT + "/derived/proteome_wide_MR/pw_MR_coloc.csv"); NIU = PW[PW.source == "Niu"]
ORDER = ST.sort_values("g", ascending=False).index.tolist()
POOL = {k: pool_re(ST[k], ST[k + "_se"]) for k in ["g", "dx", "gx"]}; POOLNS = {k: pool_re(NS[k], NS[k + "_se"]) for k in ["g", "dx", "gx"]}

# MBL2: ASD minus TD within each genotype (0, 1, 2 protein-raising G alleles at rs7899547)
rows = []
for g in [0, 1, 2]:
    a, t = MB[(MB.genotype == g) & (MB.group == "ASD")].value, MB[(MB.genotype == g) & (MB.group == "TD")].value
    tt = stats.ttest_ind(a, t, equal_var=False) if len(a) > 1 and len(t) > 1 else None
    rows.append(dict(genotype=g, n_asd=len(a), n_td=len(t), diff=a.mean() - t.mean(), welch_p=tt.pvalue if tt else np.nan))
pd.DataFrame(rows).to_csv(f"{DAT}/mbl2_within_genotype.csv", index=False)

with open(f"{DAT}/summary_statistics.txt", "w") as fh:
    fh.write(f"selection: {len(SEL)} proteins, F>10: {int((SEL.F>10).sum())}, F>10 & MAF>=0.05: {int(SEL.passB.sum())}; dropped: {', '.join(SEL[~SEL.passB].protein)}\n")
    fh.write(f"retained set exact fits: {int(ST.exact.sum())}, summary-statistic approximation: {', '.join(ST[~ST.exact].index)}\n")
    for k, v in POOL.items(): fh.write(f"retained ({len(ST)}) {k}: {v['b']:+.3f} ({v['lo']:+.3f}, {v['hi']:+.3f}) P={v['p']:.2g} I2={v['I2']:.0f}%\n")
    for k, v in POOLNS.items(): fh.write(f"non-selected ({len(NS)}) {k}: {v['b']:+.3f} ({v['lo']:+.3f}, {v['hi']:+.3f}) P={v['p']:.2g} I2={v['I2']:.0f}%\n")
    fh.write(f"variance explained median: genotype {ST.r2_g.median()*100:.1f}% diagnosis {ST.r2_dx.median()*100:.1f}%; MBL2 total {ST.loc['MBL2', ['r2_g','r2_dx']].sum()*100:.0f}%\n")
    w = 1 / ST.slope_diff_se ** 2; fh.write(f"slope diff retained (FE): {(w*ST.slope_diff).sum()/w.sum():+.3f} SE {np.sqrt(1/w.sum()):.3f}; r(ASD,TD slopes) = {np.corrcoef(ST.b_asd, ST.b_td)[0,1]:.2f}; P<0.05: {int((ST.slope_diff_p<0.05).sum())}\n")
    w = 1 / NS.slope_diff_se ** 2; fh.write(f"slope diff non-selected (FE): {(w*NS.slope_diff).sum()/w.sum():+.3f} SE {np.sqrt(1/w.sum()):.3f}; r = {np.corrcoef(NS.b_asd, NS.b_td)[0,1]:.2f}\n")
    for r in rows: fh.write(f"MBL2 genotype {r['genotype']}: ASD-TD {r['diff']:+.2f} (n {r['n_asd']}/{r['n_td']}), Welch P = {r['welch_p']:.2g}\n")
    ok = OWN.p.notna(); fh.write(f"own-instrument MR: {int(ok.sum())} tests, nominal {int((OWN.p<0.05).sum())} (expected {0.05*ok.sum():.1f}), FDR<0.05: {int((OWN.fdr<0.05).sum())}: " + "; ".join(f"{r.protein}-{r.trait} b={r.b:+.3f} P={r.p:.1e}" for _, r in OWN[OWN.fdr < 0.05].iterrows()) + "\n")
    nn = NIU[NIU.protein.isin(ORDER) & NIU.p.notna()]; fh.write(f"Niu-instrument MR on the retained set: {nn.protein.nunique()} proteins, {len(nn)} tests, nominal {int((nn.p<0.05).sum())}, FDR(all Niu tests)<0.05: {int((nn.fdr<0.05).sum())}: " + "; ".join(f"{r.protein}-{r.trait} b={r.b:+.3f} P={r.p:.1e} H4={r.PP_H4:.2f}" for _, r in nn[nn.fdr < 0.05].iterrows()) + "\n")
    a = NIU[(NIU.trait == "ASD") & NIU.p.notna()]; fh.write(f"ASD, Niu instruments: {len(a)} proteins, nominal {int((a.p<0.05).sum())}, min P {a.p.min():.1e} ({a.sort_values('p').protein.iloc[0]}), within-trait FDR<0.05: {int((a.fdr_within_trait<0.05).sum())}\n")
    for lo, hi in [(0, 0.2), (0.2, 0.4), (0.4, 9)]:
        s = CONC[(CONC.beta_niu.abs() >= lo) & (CONC.beta_niu.abs() < hi)]; fh.write(f"concordance |beta_niu| {lo}-{hi}: n={len(s)} sign {(np.sign(s.beta_niu)==np.sign(s.beta_own)).mean():.2f}\n")
    big = CONC[CONC.beta_niu.abs() >= 0.2]
    fh.write(f"concordance |beta_niu| >= 0.2: n={len(big)}, r={np.corrcoef(big.beta_niu, big.beta_own)[0,1]:.2f}, sign {(np.sign(big.beta_niu)==np.sign(big.beta_own)).mean()*100:.0f}%\n")
    fh.write(f"concordance loci {CONC.protein.nunique()}, points {len(CONC)}\n")

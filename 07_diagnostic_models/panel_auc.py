# Diagnostic performance of the 12-protein panel in the 90 unrelated individuals.
#   e  single-marker balanced AUC (best peptide per protein, 400 class-balanced subsamples)
#   f  12-peptide panel: stratified 5-fold CV x 20 repeats, peptide selection, scaling and
#      L2 logistic regression (C = 1) inside training folds only
#   g  Peptide, pQTL and WGS feature blocks and their combinations, same CV
# Inputs (relative to ASD_ROOT): pQTL/pQTL_re/pep_unique_ASD.txt, pep_unique_UHC.txt,
#   panel12/data/panel12_lead_instruments.tsv, genotype_12leads.tsv.gz, genotype_12regions.tsv.gz,
#   derived/11_protein_maps/measured_proteins_regions_1Mb.tsv (see extract_12_regions.py)
# Outputs (relative to ASD_ROOT): panel12/data/panel_e_single_marker_auc.tsv, panel_e_roc.tsv,
#   panel_f_panel_cv.json, panel_f_roc.tsv, panel_g_triomics_auc.tsv
# Usage: python panel_auc.py
import os
import csv
import gzip
import json
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve

ROOT = os.environ.get('ASD_ROOT', '.')
RE = ROOT + '/pQTL/pQTL_re'
D = ROOT + '/panel12/data'
SEED = 20260922
P12 = ['AHSG', 'ANPEP', 'BCHE', 'BTD', 'C1RL', 'C3', 'CLEC3B', 'IGFBP5',
       'MBL2', 'POSTN', 'PTGDS', 'QSOX1']
UNIPROT = {'AHSG': 'P02765', 'ANPEP': 'P15144', 'BCHE': 'P06276',
           'BTD': 'P43251', 'C1RL': 'Q9NZP8', 'C3': 'P01024',
           'CLEC3B': 'P05452', 'IGFBP5': 'P24593', 'MBL2': 'P11226',
           'POSTN': 'Q15063', 'PTGDS': 'P41222', 'QSOX1': 'O00391'}
TAB = chr(9)
NL = chr(10)
NA_STRINGS = ('', 'NA', chr(34) + 'NA' + chr(34))


def read_matrix(path):
    rows, ids = [], []
    with open(path, encoding='utf-8') as f:
        hdr = [c.strip(chr(34)) for c in f.readline().rstrip(NL).split(TAB)]
        for line in f:
            p = line.rstrip(NL).split(TAB)
            ids.append(p[0].strip(chr(34)))
            rows.append([np.nan if x in NA_STRINGS else float(x) for x in p[1:]])
    return hdr[1:], ids, np.array(rows, float)


sa, pid_a, A = read_matrix(RE + '/pep_unique_ASD.txt')
sb, pid_b, B = read_matrix(RE + '/pep_unique_UHC.txt')
assert pid_a == pid_b
samples = sa + sb
y = np.array([1] * len(sa) + [0] * len(sb))
print('peptides %d | ASD %d | TD %d' % (len(pid_a), y.sum(), (1 - y).sum()), flush=True)

# peptide-wise common z-score across both groups (as in the pQTL pipeline)
X_all = np.hstack([A, B])
mu = np.nanmean(X_all, 1, keepdims=True)
sd = np.nanstd(X_all, 1, keepdims=True)
Z = (X_all - mu) / np.where(sd > 0, sd, np.nan)

pep_of = {g: [i for i, p in enumerate(pid_a) if p.split('_')[0] == UNIPROT[g]]
          for g in P12}
for g in P12:
    print('  %-9s %d peptides' % (g, len(pep_of[g])), flush=True)
assert all(pep_of[g] for g in P12)


def balanced_auc(score, lab, n_rep=200, rng=None):
    """Mean AUC over repeated class-balanced subsamples (42 vs 42)."""
    rng = rng or np.random.default_rng(SEED)
    ok = np.isfinite(score)
    s, l = score[ok], lab[ok]
    n = min((l == 1).sum(), (l == 0).sum())
    i1, i0 = np.where(l == 1)[0], np.where(l == 0)[0]
    out = []
    for _ in range(n_rep):
        k = np.concatenate([rng.choice(i1, n, False), rng.choice(i0, n, False)])
        out.append(roc_auc_score(l[k], s[k]))
    return float(np.mean(out)), float(np.std(out))


# ---------------- panel e: best peptide per protein ----------------------
best, e_rows = {}, []
for g in P12:
    cand = []
    for i in pep_of[g]:
        s = Z[i]
        if np.isfinite(s).sum() < 60:
            continue
        a, _ = balanced_auc(s, y, 80, np.random.default_rng(SEED + i))
        cand.append((max(a, 1 - a), a, i))
    cand.sort(reverse=True)
    _, a, i = cand[0]
    sign = 1.0 if a >= 0.5 else -1.0          # orient so higher = more ASD-like
    ab, sdv = balanced_auc(sign * Z[i], y, 400, np.random.default_rng(SEED + 7))
    best[g] = (i, sign)
    ok = np.isfinite(Z[i])
    e_rows.append(dict(protein=g, peptide=pid_a[i], balanced_auc=round(ab, 4),
                       sd=round(sdv, 4),
                       plain_auc=round(roc_auc_score(y[ok], (sign * Z[i])[ok]), 4),
                       direction='up in ASD' if sign > 0 else 'down in ASD',
                       n_peptides_considered=len(cand)))
    print('  e  %-9s %-28s balanced AUC %.3f' % (g, pid_a[i], ab), flush=True)

with open(D + '/panel_e_single_marker_auc.tsv', 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, list(e_rows[0].keys()), delimiter=TAB, lineterminator=NL)
    w.writeheader()
    w.writerows(e_rows)


# ---------------- feature blocks -----------------------------------------
def impute(M):
    M = M.copy()
    for j in range(M.shape[1]):
        c = M[:, j]
        if np.isnan(c).any():
            c[np.isnan(c)] = np.nanmean(c) if np.isfinite(c).any() else 0.0
    return M


PEP = impute(np.array([best[g][1] * Z[best[g][0]] for g in P12]).T)   # 90 x 12


def read_geno(path):
    with gzip.open(path, 'rt') as f:
        cols = f.readline().rstrip(NL).split(TAB)[1:]
        ids, rows = [], []
        for line in f:
            p = line.rstrip(NL).split(TAB)
            ids.append(p[0])
            rows.append([float(x) for x in p[1:]])
    return cols, ids, np.array(rows, float)


gcols, lead_ids, LEAD = read_geno(D + '/genotype_12leads.tsv.gz')
idx = [gcols.index(s) for s in samples]
PQTL = LEAD[:, idx].T                                                  # 90 x 12

lead_of = {}
for r in csv.DictReader(open(D + '/panel12_lead_instruments.tsv', encoding='utf-8'),
                        delimiter=TAB):
    lead_of[r['variant_id']] = r['protein']
order = [lead_of[v] for v in lead_ids]
PQTL = PQTL[:, [order.index(g) for g in P12]]

regions = {}
REGF = (ROOT + '/derived/11_protein_maps/'
        'measured_proteins_regions_1Mb.tsv')
for r in csv.DictReader(open(REGF, encoding='utf-8'), delimiter=TAB):
    if r['protein'] in P12:
        regions[r['protein']] = (r['chr'], int(r['lo']), int(r['hi']))

burden = {g: np.zeros(len(samples)) for g in P12}
nvar = {g: 0 for g in P12}
with gzip.open(D + '/genotype_12regions.tsv.gz', 'rt') as f:
    cols = f.readline().rstrip(NL).split(TAB)[1:]
    ix = [cols.index(s) for s in samples]
    for line in f:
        p = line.rstrip(NL).split(TAB)
        parts = p[0].split(':')
        c, pos = parts[0], int(parts[1])
        vals = p[1:]
        d = np.array([float(vals[i]) for i in ix])
        maf = d.mean() / 2
        maf = min(maf, 1 - maf)
        if maf < 0.01:
            continue
        for g, (ch, lo, hi) in regions.items():
            if c == ch and lo <= pos <= hi:
                burden[g] += d
                nvar[g] += 1
                break

WGS = np.array([(burden[g] / max(nvar[g], 1)) for g in P12]).T
WGS = (WGS - WGS.mean(0)) / np.where(WGS.std(0) > 0, WGS.std(0), 1)
print('WGS burden variants per gene:', {g: nvar[g] for g in P12}, flush=True)

BLOCKS = {'Peptide': PEP, 'pQTL': PQTL, 'WGS': WGS}
SETS = [('Peptide', ['Peptide']),
        ('pQTL', ['pQTL']),
        ('WGS', ['WGS']),
        ('pQTL+Peptide', ['pQTL', 'Peptide']),
        ('pQTL+WGS', ['pQTL', 'WGS']),
        ('WGS+Peptide', ['WGS', 'Peptide']),
        ('Tri-omics', ['pQTL', 'WGS', 'Peptide'])]

N_SPLIT, N_REP = 5, 20


def peptide_reselect(tr):
    """Pick the best peptide per protein using the TRAINING fold only."""
    cols = []
    for g in P12:
        bi, ba, bs = None, -1.0, 1.0
        for i in pep_of[g]:
            s = Z[i][tr]
            ok = np.isfinite(s)
            if ok.sum() < 30 or len(set(y[tr][ok])) < 2:
                continue
            a = roc_auc_score(y[tr][ok], s[ok])
            if max(a, 1 - a) > ba:
                ba, bi, bs = max(a, 1 - a), i, (1.0 if a >= 0.5 else -1.0)
        if bi is None:
            bi, bs = best[g]
        cols.append(bs * Z[bi])
    return impute(np.array(cols).T)


def cv_auc(X_fixed, reselect=None):
    """Mean out-of-fold AUC over N_REP repeats of stratified N_SPLIT-fold CV."""
    cv = RepeatedStratifiedKFold(n_splits=N_SPLIT, n_repeats=N_REP,
                                 random_state=SEED)
    per_rep, oof, cur = [], np.zeros(len(y)), 0
    for tr, te in cv.split(np.zeros(len(y)), y):
        Xf = reselect(tr) if reselect is not None else X_fixed
        m = Xf[tr].mean(0)
        s = Xf[tr].std(0)
        s[s == 0] = 1.0
        clf = LogisticRegression(penalty='l2', C=1.0, max_iter=5000)
        clf.fit((Xf[tr] - m) / s, y[tr])
        oof[te] = clf.predict_proba((Xf[te] - m) / s)[:, 1]
        cur += 1
        if cur == N_SPLIT:
            per_rep.append(roc_auc_score(y, oof))
            cur = 0
            oof = np.zeros(len(y))
    a = np.array(per_rep)
    return float(a.mean()), float(a.std()), np.percentile(a, [2.5, 97.5]).tolist()


results = []
for name, blocks in SETS:
    X = np.hstack([BLOCKS[b] for b in blocks])
    if 'Peptide' in blocks:
        others = [BLOCKS[b] for b in blocks if b != 'Peptide']

        def mk(tr, o=others):
            P = peptide_reselect(tr)
            return np.hstack(o + [P]) if o else P

        m, s, ci = cv_auc(X, reselect=mk)
    else:
        m, s, ci = cv_auc(X)
    results.append(dict(feature_set=name, n_features=X.shape[1],
                        cv_auc=round(m, 4), sd=round(s, 4),
                        ci_lo=round(ci[0], 4), ci_hi=round(ci[1], 4)))
    print('  g  %-14s AUC %.3f (SD %.3f, 95%% %.3f-%.3f)'
          % (name, m, s, ci[0], ci[1]), flush=True)

with open(D + '/panel_g_triomics_auc.tsv', 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, list(results[0].keys()), delimiter=TAB, lineterminator=NL)
    w.writeheader()
    w.writerows(results)

pep_only = [r for r in results if r['feature_set'] == 'Peptide'][0]
with open(D + '/panel_f_panel_cv.json', 'w', encoding='utf-8') as f:
    json.dump(dict(n_proteins=12, n_ASD=int(y.sum()), n_TD=int((1 - y).sum()),
                   n_splits=N_SPLIT, n_repeats=N_REP, seed=SEED,
                   cv_auc=pep_only['cv_auc'], sd=pep_only['sd'],
                   ci=[pep_only['ci_lo'], pep_only['ci_hi']]), f, indent=2)

# ROC curve for the panel (pooled out-of-fold of a single 5-fold split)
cv1 = RepeatedStratifiedKFold(n_splits=N_SPLIT, n_repeats=1, random_state=SEED)
oof = np.zeros(len(y))
for tr, te in cv1.split(np.zeros(len(y)), y):
    Xf = peptide_reselect(tr)
    m = Xf[tr].mean(0)
    s = Xf[tr].std(0)
    s[s == 0] = 1.0
    clf = LogisticRegression(max_iter=5000).fit((Xf[tr] - m) / s, y[tr])
    oof[te] = clf.predict_proba((Xf[te] - m) / s)[:, 1]
fpr, tpr, _ = roc_curve(y, oof)
with open(D + '/panel_f_roc.tsv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, delimiter=TAB, lineterminator=NL)
    w.writerow(['fpr', 'tpr'])
    w.writerows(zip(fpr, tpr))
print('panel f pooled out-of-fold AUC (single split): %.4f'
      % roc_auc_score(y, oof), flush=True)

with open(D + '/panel_e_roc.tsv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, delimiter=TAB, lineterminator=NL)
    w.writerow(['protein', 'peptide', 'fpr', 'tpr'])
    for g in P12:
        i, sg = best[g]
        s = sg * Z[i]
        ok = np.isfinite(s)
        fp, tp, _ = roc_curve(y[ok], s[ok])
        for a, b in zip(fp, tp):
            w.writerow([g, pid_a[i], a, b])
print('done', flush=True)

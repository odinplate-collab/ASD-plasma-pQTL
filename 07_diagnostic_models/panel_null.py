# Reference distributions for the 12-protein panel AUC: 200 panels of 12 proteins drawn
# at random from the measured proteome, and 100 permutations of the ASD/TD labels, each run
# through the panel pipeline of panel_auc.py (stratified 5-fold CV x 4 repeats).
# Inputs (relative to ASD_ROOT): pQTL/pQTL_re/pep_unique_ASD.txt, pep_unique_UHC.txt
# Outputs (relative to ASD_ROOT): panel12/data/panel_null_summary.json, panel_null_draws.tsv
# Usage: python panel_null.py
import os
import csv
import json
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score

ROOT = os.environ.get('ASD_ROOT', '.')
RE = ROOT + '/pQTL/pQTL_re'
D = ROOT + '/panel12/data'
SEED = 20260922
TAB = chr(9)
NL = chr(10)
NA_STRINGS = ('', 'NA', chr(34) + 'NA' + chr(34))
P12 = ['AHSG', 'ANPEP', 'BCHE', 'BTD', 'C1RL', 'C3', 'CLEC3B', 'IGFBP5',
       'MBL2', 'POSTN', 'PTGDS', 'QSOX1']
UNIPROT = {'AHSG': 'P02765', 'ANPEP': 'P15144', 'BCHE': 'P06276',
           'BTD': 'P43251', 'C1RL': 'Q9NZP8', 'C3': 'P01024',
           'CLEC3B': 'P05452', 'IGFBP5': 'P24593', 'MBL2': 'P11226',
           'POSTN': 'Q15063', 'PTGDS': 'P41222', 'QSOX1': 'O00391'}
N_DRAWS = 200
N_SPLIT, N_REP = 5, 4


def read_matrix(path):
    rows, ids = [], []
    with open(path, encoding='utf-8') as f:
        hdr = [c.strip(chr(34)) for c in f.readline().rstrip(NL).split(TAB)]
        for line in f:
            p = line.rstrip(NL).split(TAB)
            ids.append(p[0].strip(chr(34)))
            rows.append([np.nan if x in NA_STRINGS else float(x) for x in p[1:]])
    return hdr[1:], ids, np.array(rows, float)


sa, pid, A = read_matrix(RE + '/pep_unique_ASD.txt')
sb, _, B = read_matrix(RE + '/pep_unique_UHC.txt')
y0 = np.array([1] * len(sa) + [0] * len(sb))
X = np.hstack([A, B])
mu = np.nanmean(X, 1, keepdims=True)
sd = np.nanstd(X, 1, keepdims=True)
Z = (X - mu) / np.where(sd > 0, sd, np.nan)

by_prot = {}
for i, p in enumerate(pid):
    if np.isfinite(Z[i]).sum() >= 60:
        by_prot.setdefault(p.split('_')[0], []).append(i)
real = [UNIPROT[g] for g in P12]
pool = sorted(by_prot)
print('proteins with usable peptides: %d' % len(pool), flush=True)


def run_panel(prot_ids, y):
    def reselect(tr):
        cols = []
        for u in prot_ids:
            bi, ba, bs = None, -1.0, 1.0
            for i in by_prot[u]:
                s = Z[i][tr]
                ok = np.isfinite(s)
                if ok.sum() < 30 or len(set(y[tr][ok])) < 2:
                    continue
                a = roc_auc_score(y[tr][ok], s[ok])
                if max(a, 1 - a) > ba:
                    ba, bi, bs = max(a, 1 - a), i, (1.0 if a >= 0.5 else -1.0)
            if bi is None:
                bi, bs = by_prot[u][0], 1.0
            cols.append(bs * Z[bi])
        M = np.array(cols).T.copy()
        for j in range(M.shape[1]):
            c = M[:, j]
            if np.isnan(c).any():
                c[np.isnan(c)] = np.nanmean(c) if np.isfinite(c).any() else 0.0
        return M

    cv = RepeatedStratifiedKFold(n_splits=N_SPLIT, n_repeats=N_REP,
                                 random_state=SEED)
    per, oof, cur = [], np.zeros(len(y)), 0
    for tr, te in cv.split(np.zeros(len(y)), y):
        M = reselect(tr)
        m = M[tr].mean(0)
        s = M[tr].std(0)
        s[s == 0] = 1.0
        clf = LogisticRegression(max_iter=5000).fit((M[tr] - m) / s, y[tr])
        oof[te] = clf.predict_proba((M[te] - m) / s)[:, 1]
        cur += 1
        if cur == N_SPLIT:
            per.append(roc_auc_score(y, oof))
            cur = 0
            oof = np.zeros(len(y))
    return float(np.mean(per))


obs = run_panel(real, y0)
print('observed 12-protein panel AUC (same settings): %.4f' % obs, flush=True)

rng = np.random.default_rng(SEED)
null_panels = []
for k in range(N_DRAWS):
    draw = list(rng.choice(pool, 12, replace=False))
    null_panels.append(run_panel(draw, y0))
    if (k + 1) % 25 == 0:
        a = np.array(null_panels)
        print('  random panels %3d/%d  mean %.3f  95%% %.3f-%.3f'
              % (k + 1, N_DRAWS, a.mean(), *np.percentile(a, [2.5, 97.5])),
              flush=True)

null_perm = []
for k in range(N_DRAWS // 2):
    yp = rng.permutation(y0)
    null_perm.append(run_panel(real, yp))
    if (k + 1) % 25 == 0:
        a = np.array(null_perm)
        print('  permuted labels %3d/%d  mean %.3f' % (k + 1, N_DRAWS // 2, a.mean()),
              flush=True)

np_ = np.array(null_panels)
pp = np.array(null_perm)
summary = dict(
    observed_auc=round(obs, 4),
    random_panel_null=dict(n=len(np_), mean=round(float(np_.mean()), 4),
                           sd=round(float(np_.std()), 4),
                           ci=[round(float(x), 4) for x in
                               np.percentile(np_, [2.5, 97.5])],
                           p_empirical=round(float((np_ >= obs).mean()), 4)),
    permuted_label_null=dict(n=len(pp), mean=round(float(pp.mean()), 4),
                             sd=round(float(pp.std()), 4),
                             ci=[round(float(x), 4) for x in
                                 np.percentile(pp, [2.5, 97.5])],
                             p_empirical=round(float((pp >= obs).mean()), 4)))
print(json.dumps(summary, indent=2), flush=True)
with open(D + '/panel_null_summary.json', 'w', encoding='utf-8') as f:
    json.dump(summary, f, indent=2)
with open(D + '/panel_null_draws.tsv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, delimiter=TAB, lineterminator=NL)
    w.writerow(['null_type', 'auc'])
    for v in np_:
        w.writerow(['random_protein_panel', round(float(v), 5)])
    for v in pp:
        w.writerow(['permuted_labels', round(float(v), 5)])
print('done', flush=True)

# Fixed-effect inverse-variance meta-analysis of the 12 panel genes across the two snRNA-seq datasets, per cell type.
# Inputs (relative to ASD_ROOT): brain/data/lr_singlecell/panel12_by_celltype_both_datasets.tsv  (sc_de_both.R)
# Outputs (relative to ASD_ROOT): brain/data/lr_singlecell/panel12_meta_two_datasets.tsv
# Usage: python panel12_meta.py
# Meta-analysed only where the gene is expressed in that cell type in BOTH datasets (column `expressed`
# of sc_de_both.R). SE = |logFC / t|; heterogeneity by Cochran's Q (1 df).
import csv
import math
import os

from scipy.stats import chi2, norm

ROOT = os.environ.get("ASD_ROOT", ".")
D = os.path.join(ROOT, 'brain', 'data', 'lr_singlecell')
TAB, NL = chr(9), chr(10)
GENES = ['C3', 'AHSG', 'MBL2', 'C1RL', 'POSTN', 'BCHE', 'BTD', 'ANPEP', 'CLEC3B',
         'IGFBP5', 'PTGDS', 'QSOX1']
CTS = ['EXN', 'INN', 'AST', 'ODC', 'OPC', 'MG', 'END']

R = list(csv.DictReader(open(D + '/panel12_by_celltype_both_datasets.tsv', encoding='utf-8'),
                        delimiter=TAB))
E = {(r['dataset'], r['celltype'], r['gene']): r for r in R}


def ok(r):
    return r is not None and r['expressed'] == 'TRUE' and r['logFC'] not in ('', 'NA')


out = []
for g in GENES:
    for ct in CTS:
        row = dict(gene=g, celltype=ct)
        est = []
        for lab, ds in [('W', 'Wamsley'), ('V', 'Velmeshev')]:
            r = E.get((ds, ct, g))
            if ok(r):
                b = float(r['logFC'])
                est.append((b, abs(b / float(r['t']))))
                row[lab + '_logFC'], row[lab + '_P'] = round(b, 3), float(r['P.Value'])
            else:
                row[lab + '_logFC'], row[lab + '_P'] = '', ''
        if len(est) == 2:
            w = [1 / s ** 2 for _, s in est]
            b = sum(x * wi for (x, _), wi in zip(est, w)) / sum(w)
            se = math.sqrt(1 / sum(w))
            q = sum(wi * (x - b) ** 2 for (x, _), wi in zip(est, w))
            row.update(meta_logFC=round(b, 3), meta_P=2 * norm.sf(abs(b / se)),
                       het_P=round(1 - chi2.cdf(q, 1), 3), same_sign=int(est[0][0] * est[1][0] > 0))
        else:
            row.update(meta_logFC='', meta_P='', het_P='', same_sign='')
        out.append(row)

tested = [r for r in out if r['meta_P'] != '']
for r in out:
    r['meta_P_bonf'] = min(1.0, r['meta_P'] * len(tested)) if r['meta_P'] != '' else ''
with open(D + '/panel12_meta_two_datasets.tsv', 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, list(out[0].keys()), delimiter=TAB, lineterminator=NL)
    w.writeheader()
    w.writerows(out)
print('testable in both: %d | same sign: %d' % (len(tested), sum(r['same_sign'] for r in tested)))
for r in sorted(tested, key=lambda r: r['meta_P'])[:6]:
    print('%-7s %-4s W %+.2f  V %+.2f  meta %+.2f  P %.4f  Bonferroni %.3f  het P %.2f'
          % (r['gene'], r['celltype'], r['W_logFC'], r['V_logFC'], r['meta_logFC'],
             r['meta_P'], r['meta_P_bonf'], r['het_P']))

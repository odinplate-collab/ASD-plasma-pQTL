# Velmeshev pseudobulk exported by gene symbol, plus counts and detecting nuclei for the ligand-receptor genes, in the Wamsley layout.
# Inputs (relative to ASD_ROOT):
#   brain/data/velmeshev/pb_counts_full.npz, pb_detect_full.npz, pb_genes.tsv   (velmeshev_pseudobulk.py)
#   brain/data/lr_singlecell/wamsley_lr_counts.tsv                              (wamsley_lr_fetch.py; gene list)
# Outputs (relative to ASD_ROOT):
#   brain/data/velmeshev/pb_counts_symbol.tsv.gz   genes (symbol; duplicate symbols summed; genes with
#                                                  zero total dropped) x groups
#   brain/data/lr_singlecell/velmeshev_lr_counts.tsv, velmeshev_lr_detect.tsv
# Usage: python velmeshev_export.py
import csv
import gzip
import os
from collections import OrderedDict

import numpy as np

ROOT = os.environ.get("ASD_ROOT", ".")
HERE = os.path.join(ROOT, 'brain')
V = HERE + '/data/velmeshev'
OUT = HERE + '/data/lr_singlecell'
TAB, NL = chr(9), chr(10)

C = np.load(V + '/pb_counts_full.npz')['counts']
Dt = np.load(V + '/pb_detect_full.npz')['detect']
genes = [r['gene'] for r in csv.DictReader(open(V + '/pb_genes.tsv', encoding='utf-8'), delimiter=TAB)]
NG = C.shape[1]

rows = OrderedDict()
for i, g in enumerate(genes):
    rows.setdefault(g, []).append(i)
kept = 0
with gzip.open(V + '/pb_counts_symbol.tsv.gz', 'wt', newline='') as f:
    w = csv.writer(f, delimiter=TAB, lineterminator=NL)
    w.writerow(['gene'] + ['g%d' % i for i in range(NG)])
    for g, ix in rows.items():
        v = C[ix].sum(axis=0)
        if v.sum() == 0:
            continue
        w.writerow([g] + ['%d' % x for x in v])
        kept += 1
print('symbols written %d (of %d unique, %d rows)' % (kept, len(rows), len(genes)))

lr = [l.split(TAB)[0] for l in open(OUT + '/wamsley_lr_counts.tsv', encoding='utf-8')][1:]
lr = ['MBL2', 'IFITM1', 'IGF2'] + lr   # not in the Wamsley matrix, present here
lr = list(dict.fromkeys(lr))
for arr, name in [(C, 'velmeshev_lr_counts.tsv'), (Dt, 'velmeshev_lr_detect.tsv')]:
    with open(OUT + '/' + name, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter=TAB, lineterminator=NL)
        w.writerow(['gene'] + ['g%d' % i for i in range(NG)])
        for g in lr:
            if g not in rows:
                print('not in Velmeshev:', g)
                continue
            w.writerow([g] + ['%d' % x for x in arr[rows[g]].sum(axis=0)])
print('LR genes written')

# Wamsley et al. 2024 snRNA-seq: summed UMI counts and number of detecting nuclei per pseudobulk group for the 12 panel genes and their curated receptors.
# Inputs (relative to ASD_ROOT):
#   brain/data/lr_pairs_panel12_cellchat.tsv   (lr_pairs_cellchat.R)
#   brain/data/lr_pairs_panel12_omnipath.tsv   (lr_pairs_omnipath.py)
#   brain/data/wamsley/meta.tsv, exprMatrix.json, pb_groups.tsv   (wamsley_pseudobulk.py)
#   https://cells.ucsc.edu/asd-psychencode/exprMatrix.bin   (HTTP range requests)
# Outputs (relative to ASD_ROOT):
#   brain/data/lr_singlecell/wamsley_lr_counts.tsv   summed UMI counts, genes x groups
#   brain/data/lr_singlecell/wamsley_lr_detect.tsv   number of nuclei with count > 0, genes x groups
# Usage: python wamsley_lr_fetch.py
# Genes: the 12 panel genes, every receptor subunit from CellChatDB.human and OmniPath, plus IGF1R and
# IGF2R (receptors of the IGFs bound by IGFBP5) and VSIG4 and C5AR1 (microglial complement receptors).
# Group order is that of pb_groups.tsv (rebuilt from meta.tsv the same way and checked).
import csv
import json
import os
import struct
import time
import urllib.request
import zlib

import numpy as np

ROOT = os.environ.get("ASD_ROOT", ".")
HERE = os.path.join(ROOT, 'brain')
W = HERE + '/data/wamsley'
OUT = HERE + '/data/lr_singlecell'
BIN = 'https://cells.ucsc.edu/asd-psychencode/exprMatrix.bin'
TAB, NL = chr(9), chr(10)

P12 = ['AHSG', 'ANPEP', 'BCHE', 'BTD', 'C1RL', 'C3', 'CLEC3B', 'IGFBP5',
       'MBL2', 'POSTN', 'PTGDS', 'QSOX1']
rec = set()
for fn, col in [('lr_pairs_panel12_cellchat.tsv', 'receptor_genes'),
                ('lr_pairs_panel12_omnipath.tsv', 'receptor')]:
    for r in csv.DictReader(open(HERE + '/data/' + fn, encoding='utf-8'), delimiter=TAB):
        rec.update(r[col].split('+'))
EXTRA = ['IGF1R', 'IGF2R', 'VSIG4', 'C5AR1']
want = P12 + sorted(rec) + EXTRA

IDX = json.load(open(W + '/exprMatrix.json', encoding='utf-8'))
missing = [g for g in want if g not in IDX]
want = [g for g in dict.fromkeys(want) if g in IDX]
print('genes to fetch %d | not in Wamsley matrix: %s' % (len(want), ', '.join(missing)), flush=True)

keys, key_of, grp = [], {}, []
with open(W + '/meta.tsv', encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter=TAB):
        k = (r['individual_ID'], r['Brain_Region'], r['annotation'])
        if k not in key_of:
            key_of[k] = len(keys)
            keys.append(k)
        grp.append(key_of[k])
grp = np.array(grp, dtype=np.int64)
NCELL, NG = len(grp), len(keys)
ref = [(r['individual_ID'], r['region'], r['celltype'])
       for r in csv.DictReader(open(W + '/pb_groups.tsv', encoding='utf-8'), delimiter=TAB)]
assert ref == keys, 'group order differs from pb_groups.tsv'


def fetch(g, tries=6):
    off, ln = IDX[g]
    for t in range(tries):
        try:
            req = urllib.request.Request(BIN, headers={
                'Range': 'bytes=%d-%d' % (off, off + ln - 1), 'User-Agent': 'curl/8'})
            raw = urllib.request.urlopen(req, timeout=180).read()
            if len(raw) != ln:
                raise IOError('short read')
            return zlib.decompress(raw)
        except Exception:
            if t == tries - 1:
                raise
            time.sleep(5 * (t + 1))


CNT = np.zeros((len(want), NG))
DET = np.zeros((len(want), NG), dtype=np.int64)
for j, g in enumerate(want):
    dec = fetch(g)
    nlen = struct.unpack('<H', dec[:2])[0]
    assert dec[2:2 + nlen].decode('utf-8') == g
    v = np.frombuffer(dec[2 + nlen:], dtype='<u4')
    assert len(v) == NCELL
    CNT[j] = np.bincount(grp, weights=v.astype(np.float64), minlength=NG)
    DET[j] = np.bincount(grp, weights=(v > 0).astype(np.float64), minlength=NG)
    print('  %2d/%d %s' % (j + 1, len(want), g), flush=True)

os.makedirs(OUT, exist_ok=True)
for arr, name in [(CNT, 'wamsley_lr_counts.tsv'), (DET, 'wamsley_lr_detect.tsv')]:
    with open(OUT + '/' + name, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter=TAB, lineterminator=NL)
        w.writerow(['gene'] + ['g%d' % i for i in range(NG)])
        for g, row in zip(want, arr):
            w.writerow([g] + ['%d' % x for x in row])
print('saved %d genes x %d groups' % CNT.shape, flush=True)

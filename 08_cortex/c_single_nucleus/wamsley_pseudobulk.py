# Wamsley et al. 2024 snRNA-seq (UCSC Cell Browser 'asd-psychencode', 591,043 nuclei): raw UMI counts summed per donor x region x cell type, gene by gene.
# Inputs (relative to ASD_ROOT):
#   brain/data/wamsley/meta.tsv          cell metadata, https://cells.ucsc.edu/asd-psychencode/meta.tsv
#   brain/data/wamsley/exprMatrix.json   gene -> byte range index, https://cells.ucsc.edu/asd-psychencode/exprMatrix.json
#   https://cells.ucsc.edu/asd-psychencode/exprMatrix.bin   (read by HTTP range requests, one gene at a time)
#   panel12/data/gene_modules.tsv        (../b_cortex_transcriptome/module_gene_sets.py)
# Outputs (relative to ASD_ROOT):
#   brain/data/wamsley/pb_groups.tsv     one row per pseudobulk group (donor x region x cell type)
#   brain/data/wamsley/pb_counts.npy, pb_counts.tsv.gz   summed UMI counts, genes x groups
#   brain/data/wamsley/pb_genes.tsv      fetched genes (module or background)
# Usage: python wamsley_pseudobulk.py
# Genes: all pathway-module genes present in the matrix plus a systematic background sample
# (every 26th gene of the sorted gene index, about 1,200 genes). Library size = summed nCount_RNA.
# Each byte range is a zlib stream that decodes to uint16 name_length | name | uint32 x n_cells.
import csv
import gzip
import json
import os
import struct
import time
import urllib.request
import zlib
from collections import defaultdict

import numpy as np

ROOT = os.environ.get("ASD_ROOT", ".")
D = os.path.join(ROOT, 'brain', 'data', 'wamsley')
MODULES = os.path.join(ROOT, 'panel12', 'data')
BIN = 'https://cells.ucsc.edu/asd-psychencode/exprMatrix.bin'
TAB, NL = chr(9), chr(10)
BG_STEP = 26            # 31,113 / 26 -> about 1,200 background genes

# ---- metadata -> group index per cell ---------------------------------------
print('reading meta.tsv ...', flush=True)
keys, key_of = [], {}
cell_group = []
lib = defaultdict(float)
info = {}
with open(D + '/meta.tsv', encoding='utf-8') as f:
    rd = csv.DictReader(f, delimiter=TAB)
    for r in rd:
        k = (r['individual_ID'], r['Brain_Region'], r['annotation'])
        if k not in key_of:
            key_of[k] = len(keys)
            keys.append(k)
        gi = key_of[k]
        cell_group.append(gi)
        try:
            lib[gi] += float(r['nCount_RNA'])
        except ValueError:
            pass
        info[r['individual_ID']] = (r['Diagnosis'], r['Age'], r['Sex_Chromosome'],
                                    r['Brain_Brank'])
cell_group = np.array(cell_group, dtype=np.int64)
NCELL, NG = len(cell_group), len(keys)
print('cells %d | pseudobulk groups %d | donors %d' % (NCELL, NG, len(info)), flush=True)

with open(D + '/pb_groups.tsv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, delimiter=TAB, lineterminator=NL)
    w.writerow(['group', 'individual_ID', 'region', 'celltype', 'n_cells', 'lib_size',
                'diagnosis', 'age', 'sex', 'brain_bank'])
    ncell = np.bincount(cell_group, minlength=NG)
    for i, (ind, reg, ct) in enumerate(keys):
        dx, age, sex, bank = info[ind]
        w.writerow([i, ind, reg, ct, int(ncell[i]), int(lib[i]), dx, age, sex, bank])

# ---- gene list --------------------------------------------------------------
IDX = json.load(open(D + '/exprMatrix.json', encoding='utf-8'))
mod_genes = set()
for r in csv.DictReader(open(MODULES + '/gene_modules.tsv', encoding='utf-8'), delimiter=TAB):
    mod_genes.add(r['gene'])
allg = sorted(IDX)
bg = [g for i, g in enumerate(allg) if i % BG_STEP == 0 and g not in mod_genes]
want = [g for g in sorted(mod_genes) if g in IDX] + bg
print('module genes present %d | background %d | total %d'
      % (sum(g in IDX for g in mod_genes), len(bg), len(want)), flush=True)


def fetch(g, tries=6):
    off, ln = IDX[g]
    for t in range(tries):
        try:
            req = urllib.request.Request(BIN, headers={
                'Range': 'bytes=%d-%d' % (off, off + ln - 1), 'User-Agent': 'curl/8'})
            raw = urllib.request.urlopen(req, timeout=180).read()
            if len(raw) != ln:
                raise IOError('short read %d/%d' % (len(raw), ln))
            return zlib.decompress(raw)
        except Exception as e:
            if t == tries - 1:
                raise
            time.sleep(5 * (t + 1))


# ---- fetch and aggregate ----------------------------------------------------
PB = np.zeros((len(want), NG), dtype=np.float64)
done_names = []
t0 = time.time()
for j, g in enumerate(want):
    dec = fetch(g)
    nlen = struct.unpack('<H', dec[:2])[0]
    name = dec[2:2 + nlen].decode('utf-8')
    vals = np.frombuffer(dec[2 + nlen:], dtype='<u4')
    if len(vals) != NCELL:
        raise SystemExit('cell count mismatch for %s: %d vs %d' % (g, len(vals), NCELL))
    if name != g:
        raise SystemExit('name mismatch: asked %s got %s' % (g, name))
    PB[j] = np.bincount(cell_group, weights=vals.astype(np.float64), minlength=NG)
    done_names.append(g)
    if (j + 1) % 50 == 0 or j < 3:
        el = time.time() - t0
        print('  %4d/%d  %-10s  %.0fs elapsed, ~%.0f min left'
              % (j + 1, len(want), g, el, el / (j + 1) * (len(want) - j - 1) / 60), flush=True)

np.save(D + '/pb_counts.npy', PB)
with gzip.open(D + '/pb_counts.tsv.gz', 'wt', newline='') as f:
    w = csv.writer(f, delimiter=TAB, lineterminator=NL)
    w.writerow(['gene'] + ['g%d' % i for i in range(NG)])
    for g, row in zip(done_names, PB):
        w.writerow([g] + ['%d' % v for v in row])
with open(D + '/pb_genes.tsv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, delimiter=TAB, lineterminator=NL)
    w.writerow(['gene', 'set'])
    for g in done_names:
        w.writerow([g, 'module' if g in mod_genes else 'background'])
print('saved pseudobulk: %d genes x %d groups' % PB.shape, flush=True)

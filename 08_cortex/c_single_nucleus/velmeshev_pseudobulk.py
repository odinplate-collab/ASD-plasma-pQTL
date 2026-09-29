# Velmeshev et al. 2019 snRNA-seq (104,559 nuclei, raw UMI counts): summed counts and detecting nuclei per sample x cluster pseudobulk, all genes.
# Inputs (relative to ASD_ROOT):
#   single_cell/rawMatrix/matrix.mtx, genes.tsv, barcodes.tsv, meta.txt
#     (raw count matrix and cell metadata from the UCSC Cell Browser, dataset 'autism')
# Outputs (relative to ASD_ROOT):
#   brain/data/velmeshev/pb_groups.tsv       one row per sample x cluster pseudobulk
#   brain/data/velmeshev/pb_genes.tsv        ensembl id, symbol (row order of the matrices)
#   brain/data/velmeshev/pb_counts_full.npz  summed UMI counts, all genes x all groups
#   brain/data/velmeshev/pb_detect_full.npz  number of nuclei with count > 0, same shape
# Usage: python velmeshev_pseudobulk.py
# The 3.3 GB matrix.mtx is streamed in chunks; each chunk is summed into the pseudobulk and dropped.
import csv
import os
import time
from collections import defaultdict

import numpy as np
import pandas as pd

ROOT = os.environ.get("ASD_ROOT", ".")
RAW = os.path.join(ROOT, 'single_cell', 'rawMatrix')
OUT = os.path.join(ROOT, 'brain', 'data', 'velmeshev')
TAB, NL = chr(9), chr(10)
CHUNK = 20_000_000

# broad classes matched to the seven Wamsley annotations; NRGN and maturing neurons have no
# Wamsley counterpart and are kept as their own classes
BROAD = {'L2/3': 'EXN', 'L4': 'EXN', 'L5/6': 'EXN', 'L5/6-CC': 'EXN',
         'IN-PV': 'INN', 'IN-SST': 'INN', 'IN-SV2C': 'INN', 'IN-VIP': 'INN',
         'AST-PP': 'AST', 'AST-FB': 'AST', 'Oligodendrocytes': 'ODC', 'OPC': 'OPC',
         'Microglia': 'MG', 'Endothelial': 'END',
         'Neu-NRGN-I': 'NRGN', 'Neu-NRGN-II': 'NRGN', 'Neu-mat': 'NEUMAT'}

genes = [l.rstrip(NL).split(TAB) for l in open(RAW + '/genes.tsv', encoding='utf-8')]
bcs = [l.strip() for l in open(RAW + '/barcodes.tsv', encoding='utf-8')]
meta = list(csv.DictReader(open(RAW + '/meta.txt', encoding='utf-8'), delimiter=TAB))
assert [m['cell'] for m in meta] == bcs, 'barcode order differs from meta'
NGENE, NCELL = len(genes), len(bcs)

keys, key_of, info = [], {}, {}
cell_group = np.empty(NCELL, dtype=np.int64)
umi = defaultdict(float)
for i, m in enumerate(meta):
    k = (m['sample'], m['cluster'])
    if k not in key_of:
        key_of[k] = len(keys)
        keys.append(k)
        info[k] = m
    cell_group[i] = key_of[k]
    umi[key_of[k]] += float(m['UMIs'])
NG = len(keys)
print('genes %d | nuclei %d | pseudobulk groups %d' % (NGENE, NCELL, NG), flush=True)

acc_c = np.zeros(NGENE * NG, dtype=np.float64)
acc_d = np.zeros(NGENE * NG, dtype=np.int64)
t0, nnz = time.time(), 0
rd = pd.read_csv(RAW + '/matrix.mtx', sep=' ', header=None, skiprows=3,
                 usecols=[0, 1, 2], dtype=np.int64, chunksize=CHUNK, engine='c')
for ch in rd:
    a = ch.to_numpy()
    idx = (a[:, 0] - 1) * NG + cell_group[a[:, 1] - 1]
    acc_c += np.bincount(idx, weights=a[:, 2].astype(np.float64), minlength=NGENE * NG)
    acc_d += np.bincount(idx, minlength=NGENE * NG)
    nnz += len(a)
    print('  %11d entries  %.0fs' % (nnz, time.time() - t0), flush=True)
print('total entries %d' % nnz, flush=True)

C = acc_c.reshape(NGENE, NG)
Dt = acc_d.reshape(NGENE, NG)
colsum = C.sum(axis=0)
lib = np.array([umi[i] for i in range(NG)])
print('matrix column sums vs meta UMIs: max abs diff %.1f, max rel diff %.2e'
      % (np.abs(colsum - lib).max(), (np.abs(colsum - lib) / lib).max()), flush=True)

os.makedirs(OUT, exist_ok=True)
np.savez_compressed(OUT + '/pb_counts_full.npz', counts=C.astype(np.int64))
np.savez_compressed(OUT + '/pb_detect_full.npz', detect=Dt.astype(np.int32))
with open(OUT + '/pb_genes.tsv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, delimiter=TAB, lineterminator=NL)
    w.writerow(['ensembl', 'gene'])
    for g in genes:
        w.writerow(g[:2])
ncell = np.bincount(cell_group, minlength=NG)
with open(OUT + '/pb_groups.tsv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, delimiter=TAB, lineterminator=NL)
    w.writerow(['group', 'sample', 'individual_ID', 'region', 'cluster', 'celltype',
                'n_cells', 'lib_size', 'lib_size_matrix', 'diagnosis', 'age', 'sex',
                'capbatch', 'seqbatch', 'pmi', 'rin'])
    for i, (s, cl) in enumerate(keys):
        m = info[(s, cl)]
        w.writerow([i, s, m['individual'], m['region'], cl, BROAD[cl], int(ncell[i]),
                    int(lib[i]), int(colsum[i]), m['diagnosis'], m['age'], m['sex'],
                    m['Capbatch'], m['Seqbatch'],
                    m['post-mortem interval (hours)'], m['RNA Integrity Number']])
print('saved. %.0fs' % (time.time() - t0), flush=True)

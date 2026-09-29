# Donor-level QC of the two snRNA-seq datasets: sex from XIST / Y-gene expression, Dup15q status and donor overlap from Gandal et al. 2022.
# Inputs (relative to ASD_ROOT):
#   brain/data/asd_cortex/Gandal2022_MOESM3.xlsx   Gandal et al. 2022 subject metadata (sheet datMeta_datSeq)
#   brain/data/wamsley/meta.tsv, exprMatrix.json, pb_groups.tsv   (wamsley_pseudobulk.py)
#   https://cells.ucsc.edu/asd-psychencode/exprMatrix.bin          (XIST, RPS4Y1, UTY by HTTP range requests)
#   brain/data/velmeshev/pb_counts_full.npz, pb_genes.tsv, pb_groups.tsv   (velmeshev_pseudobulk.py)
# Outputs (relative to ASD_ROOT):
#   brain/data/lr_singlecell/donor_qc_wamsley.tsv, donor_qc_velmeshev.tsv
#     per donor: recorded sex, XIST / RPS4Y1 / UTY CPM, sex from expression, match in Gandal 2022
#     (subject, diagnosis, sex), dup15q flag, presence in the other snRNA-seq dataset
# Usage: python donor_qc.py
# The Wamsley browser metadata records diagnosis only as ASD / CTL; donors listed as Dup15q in
# Gandal 2022 (same donor number, age within 1.5 years) are flagged and excluded in sc_de_both.R.
import csv
import json
import os
import re
import struct
import urllib.request
import zlib
from collections import Counter, defaultdict

import numpy as np
import openpyxl

ROOT = os.environ.get("ASD_ROOT", ".")
D = os.path.join(ROOT, 'brain', 'data')
OUT = D + '/lr_singlecell'
TAB, NL = chr(9), chr(10)
SEXG = ['XIST', 'RPS4Y1', 'UTY']

# ---- Gandal 2022 subjects ---------------------------------------------------
wb = openpyxl.load_workbook(D + '/asd_cortex/Gandal2022_MOESM3.xlsx', read_only=True)
rows = list(wb['datMeta_datSeq'].iter_rows(values_only=True))
h = list(rows[0])
GAN = {}
for r in rows[1:]:
    s = str(r[h.index('subject')])
    num = re.findall(r'\d+', s)
    GAN.setdefault(num[-1].lstrip('0') if num else s,
                   dict(subject=s, dx=r[h.index('Diagnosis')], sex=r[h.index('Sex')],
                        age=r[h.index('Age')]))


def gandal(donor, age):
    g = GAN.get(donor.lstrip('0'))
    if not g:
        return None
    try:
        if abs(float(g['age']) - float(age)) > 1.5:
            return None      # same digits, different person
    except (TypeError, ValueError):
        return None
    return g


def sex_call(x, y, u):
    # x, y, u: donor-level CPM of XIST, RPS4Y1, UTY
    if x > 100 and u < 20:
        return 'XX'
    if x < 50 and u > 50:
        return 'XY'
    return 'ambiguous'


# ---- Wamsley: fetch the sex genes -------------------------------------------
IDX = json.load(open(D + '/wamsley/exprMatrix.json', encoding='utf-8'))
don, info = [], {}
with open(D + '/wamsley/meta.tsv', encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter=TAB):
        don.append(r['individual_ID'])
        info[r['individual_ID']] = (r['Diagnosis'], r['Sex_Chromosome'], r['Age'])
don = np.array(don)
lib = Counter()
for r in csv.DictReader(open(D + '/wamsley/pb_groups.tsv', encoding='utf-8'), delimiter=TAB):
    lib[r['individual_ID']] += float(r['lib_size'])
W = {}
for g in SEXG:
    off, ln = IDX[g]
    req = urllib.request.Request('https://cells.ucsc.edu/asd-psychencode/exprMatrix.bin',
                                 headers={'Range': 'bytes=%d-%d' % (off, off + ln - 1),
                                          'User-Agent': 'curl/8'})
    dec = zlib.decompress(urllib.request.urlopen(req, timeout=180).read())
    n = struct.unpack('<H', dec[:2])[0]
    v = np.frombuffer(dec[2 + n:], dtype='<u4').astype(float)
    W[g] = {d: v[don == d].sum() / lib[d] * 1e6 for d in lib}

# ---- Velmeshev: from the full pseudobulk ------------------------------------
C = np.load(D + '/velmeshev/pb_counts_full.npz')['counts']
vg = [r['gene'] for r in csv.DictReader(open(D + '/velmeshev/pb_genes.tsv', encoding='utf-8'), delimiter=TAB)]
vgr = [r for r in csv.DictReader(open(D + '/velmeshev/pb_groups.tsv', encoding='utf-8'), delimiter=TAB)]
vinfo, vcol = {}, defaultdict(list)
for i, r in enumerate(vgr):
    vinfo[r['individual_ID']] = (r['diagnosis'], r['sex'], r['age'])
    vcol[r['individual_ID']].append(i)
V = {g: {} for g in SEXG}
for d, ix in vcol.items():
    tot = C[:, ix].sum()
    for g in SEXG:
        V[g][d] = C[[k for k, s in enumerate(vg) if s == g], :][:, ix].sum() / tot * 1e6

# ---- tables -------------------------------------------------------------------
os.makedirs(OUT, exist_ok=True)
for name, info_, E, other in [('wamsley', info, W, set(vinfo)), ('velmeshev', vinfo, V, set(info))]:
    out = []
    for d in sorted(info_):
        dx, sx, age = info_[d]
        g = gandal(d, age)
        se = sex_call(E['XIST'][d], E['RPS4Y1'][d], E['UTY'][d])
        sm = {'M': 'XY', 'F': 'XX'}.get(sx, sx)
        sm = 'XY' if sm == 'XYY' else sm
        out.append(dict(donor=d, dx_meta=dx, age=age, sex_meta=info_[d][1],
                        XIST_cpm=round(E['XIST'][d], 1), RPS4Y1_cpm=round(E['RPS4Y1'][d], 1),
                        UTY_cpm=round(E['UTY'][d], 1), sex_expr=se,
                        sex_mismatch=int(se != 'ambiguous' and se != sm),
                        in_gandal2022=g['subject'] if g else '',
                        gandal_dx=g['dx'] if g else '', gandal_sex=g['sex'] if g else '',
                        dup15q=int(bool(g) and g['dx'] == 'Dup15q'),
                        in_other_sn_dataset=int(d in other)))
    with open(OUT + '/donor_qc_%s.tsv' % name, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, list(out[0].keys()), delimiter=TAB, lineterminator=NL)
        w.writeheader()
        w.writerows(out)
    print('%s: donors %d | sex mismatch %d | ambiguous %d | Dup15q %d | in Gandal2022 %d | in other sn %d'
          % (name, len(out), sum(o['sex_mismatch'] for o in out),
             sum(o['sex_expr'] == 'ambiguous' for o in out), sum(o['dup15q'] for o in out),
             sum(bool(o['in_gandal2022']) for o in out), sum(o['in_other_sn_dataset'] for o in out)))
    for o in out:
        if o['sex_mismatch'] or o['sex_expr'] == 'ambiguous' or o['dup15q'] or o['in_other_sn_dataset']:
            print('   ', {k: o[k] for k in ['donor', 'dx_meta', 'age', 'sex_meta', 'XIST_cpm', 'UTY_cpm',
                                          'sex_expr', 'gandal_dx', 'gandal_sex', 'in_other_sn_dataset']})

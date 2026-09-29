# Ligand-receptor pairs from OmniPath (ligrecextra) with one of the 12 panel proteins as the source (ligand).
# Inputs (relative to ASD_ROOT): none (queries https://omnipathdb.org; the table used was retrieved on 2026-09-23)
# Outputs (relative to ASD_ROOT): brain/data/lr_pairs_panel12_omnipath.tsv  (ligand, receptor, n_resources)
# Usage: python lr_pairs_omnipath.py
# n_resources = number of ';'-separated entries in the OmniPath 'sources' field of the interaction.
import collections
import csv
import io
import os
import urllib.request

ROOT = os.environ.get("ASD_ROOT", ".")
OUT = os.path.join(ROOT, 'brain', 'data', 'lr_pairs_panel12_omnipath.tsv')
P12 = ['AHSG', 'ANPEP', 'BCHE', 'BTD', 'C1RL', 'C3', 'CLEC3B', 'IGFBP5', 'MBL2', 'POSTN', 'PTGDS', 'QSOX1']
u = ('https://omnipathdb.org/interactions?genesymbols=yes&datasets=ligrecextra'
     '&sources=' + ','.join(P12) + '&fields=sources&format=tsv')
txt = urllib.request.urlopen(u, timeout=180).read().decode()
rows = list(csv.DictReader(io.StringIO(txt), delimiter='\t'))
by = collections.defaultdict(list)
for r in rows:
    by[r['source_genesymbol']].append((r['target_genesymbol'], len(r.get('sources', '').split(';'))))

# ligands in order of first appearance in the response, receptors sorted within each ligand
with open(OUT, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['ligand', 'receptor', 'n_resources'])
    for lig, pairs in by.items():
        for rec, n in sorted(pairs):
            w.writerow([lig, rec, n])
print('interactions %d | ligands %d' % (len(rows), len(by)))

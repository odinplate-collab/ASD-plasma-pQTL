# Per panel protein, the plasma-CSF and plasma-brain colocalization test with the highest PP.H4.
# Inputs (relative to ASD_ROOT):
#   brain/data/coloc_plasma_vs_csf.tsv          (coloc_plasma_csf.R)
#   brain/data/coloc_plasma_vs_brain_eqtl.tsv   (coloc_plasma_brain.R)
# Outputs (relative to ASD_ROOT): brain/results/coloc_best_per_gene.tsv
# Usage: python coloc_best_per_gene.py
import csv
import math
import os

ROOT = os.environ.get("ASD_ROOT", ".")
BL = os.path.join(ROOT, 'brain')
OUT = os.path.join(BL, 'results')
TAB, NL = chr(9), chr(10)
P12 = ['C3', 'AHSG', 'MBL2', 'C1RL', 'BCHE', 'BTD', 'ANPEP', 'IGFBP5', 'POSTN', 'CLEC3B',
       'QSOX1', 'PTGDS']


def write(name, rows):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, list(rows[0].keys()), delimiter=TAB, lineterminator=NL)
        w.writeheader()
        w.writerows(rows)
    print('wrote %-30s %4d rows' % (name, len(rows)))


def rd(path):
    return list(csv.DictReader(open(path, encoding='utf-8'), delimiter=TAB))


csf = rd(os.path.join(BL, 'data', 'coloc_plasma_vs_csf.tsv'))
eq = rd(os.path.join(BL, 'data', 'coloc_plasma_vs_brain_eqtl.tsv'))
rows = []
for g in P12:
    # CSF: several assays per protein possible; keep the test with the highest PP.H4
    c = [r for r in csf if r['gene'] == g]
    if c:
        best = max(c, key=lambda r: float(r['PP.H4']))
        rows.append(dict(gene=g, comparison='CSF pQTL', n_tests=len(c), best_test=best['accession'],
                         plasma_max_nlp=best['max_nlp_plasma'], other_max_nlp=best['max_nlp_csf'],
                         **{'PP.H%d' % k: float(best['PP.H%d' % k]) for k in range(5)}))
    else:
        rows.append(dict(gene=g, comparison='CSF pQTL', n_tests=0, best_test='', plasma_max_nlp='',
                         other_max_nlp='', **{'PP.H%d' % k: '' for k in range(5)}))
    # brain eQTL: up to 20 datasets per gene; keep the test with the highest PP.H4
    e = [r for r in eq if r['gene'] == g]
    if e:
        best = max(e, key=lambda r: float(r['PP.H4']))
        rows.append(dict(gene=g, comparison='Brain eQTL', n_tests=len(e),
                         best_test=best['study'] + ' ' + best['tissue'],
                         plasma_max_nlp=round(-math.log10(float(best['min_p_plasma'])), 2),
                         other_max_nlp=round(-math.log10(min(float(r['min_p_brain']) for r in e)), 2),
                         **{'PP.H%d' % k: float(best['PP.H%d' % k]) for k in range(5)}))
    else:
        rows.append(dict(gene=g, comparison='Brain eQTL', n_tests=0, best_test='', plasma_max_nlp='',
                         other_max_nlp='', **{'PP.H%d' % k: '' for k in range(5)}))
write('coloc_best_per_gene.tsv', rows)
for r in rows:
    if r['n_tests']:
        print('%-7s %-10s tests %2s  H3 %.2f  H4 %.2f  plasma -log10P %s' %
              (r['gene'], r['comparison'], r['n_tests'], r['PP.H3'], r['PP.H4'], r['plasma_max_nlp']))
    else:
        print('%-7s %-10s not tested' % (r['gene'], r['comparison']))

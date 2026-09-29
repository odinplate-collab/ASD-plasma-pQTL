# Same module tests as gandal2022_primary.py on the cortex microarray meta-analyses of Gandal et al. 2018 (ASD, SCZ, BD, MDD; IBD bowel as non-brain reference).
# Inputs (relative to ASD_ROOT):
#   brain/data/replication/gandal2018a/TableS1_Microarray_MetaAnalyses.csv   Gandal et al. 2018 Science, Table S1
#   panel12/data/gene_modules.tsv                                           (module_gene_sets.py)
# Outputs (relative to ASD_ROOT): brain/data/replication/replication_modules.tsv
# Usage: python replication_modules.py
import csv
import os
import statistics
from collections import defaultdict

from scipy.stats import mannwhitneyu, binomtest

ROOT = os.environ.get("ASD_ROOT", ".")
D = os.path.join(ROOT, 'brain', 'data', 'replication')
MODULES = os.path.join(ROOT, 'panel12', 'data')
TAB, NL = chr(9), chr(10)
P12 = ['AHSG', 'ANPEP', 'BCHE', 'BTD', 'C1RL', 'C3', 'CLEC3B', 'IGFBP5',
       'MBL2', 'POSTN', 'PTGDS', 'QSOX1']

MOD = defaultdict(list)
for r in csv.DictReader(open(MODULES + '/gene_modules.tsv', encoding='utf-8'), delimiter=TAB):
    MOD[r['module']].append(r['gene'])
CMP = sorted({g for m in MOD if m.startswith('CMP_') for g in MOD[m]})
MOD['COMPLEMENT_all_layers'] = CMP


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def run(label, lfc, fdr):
    """lfc, fdr: dict gene -> value. Returns list of module rows."""
    genes = [g for g in lfc if lfc[g] is not None]
    allv = {g: lfc[g] for g in genes}
    bg = list(allv.values())
    print(NL + '#' * 72)
    print('%s  |  %d genes  |  background median %+.4f, %.1f%% > 0'
          % (label, len(bg), statistics.median(bg), 100 * sum(v > 0 for v in bg) / len(bg)))
    print('#' * 72)
    print('%-8s %8s %10s' % ('gene', 'logFC', 'FDR'))
    for g in P12:
        if g in allv:
            print('%-8s %+8.3f %10.2e%s' % (g, allv[g], fdr.get(g) or 1,
                                          ' *' if (fdr.get(g) or 1) < 0.05 else ''))
        else:
            print('%-8s  not measured' % g)
    out = []
    print(NL + '%-30s %4s %9s %6s %9s %10s %10s'
          % ('module', 'n', 'med lFC', '% up', 'sig u/d', 'binom P', 'compet P'))
    order = ['COMPLEMENT_all_layers'] + sorted(m for m in MOD if m != 'COMPLEMENT_all_layers')
    for m in order:
        gi = [g for g in MOD[m] if g in allv]
        if len(gi) < 3:
            continue
        v = [allv[g] for g in gi]
        s = set(gi)
        rest = [allv[g] for g in allv if g not in s]
        up = sum(x > 0 for x in v)
        su = sum(1 for g in gi if allv[g] > 0 and (fdr.get(g) or 1) < 0.05)
        sd = sum(1 for g in gi if allv[g] < 0 and (fdr.get(g) or 1) < 0.05)
        pb = binomtest(up, len(v), 0.5).pvalue
        pc = mannwhitneyu(v, rest, alternative='two-sided', method='asymptotic').pvalue
        out.append(dict(dataset=label, module=m, n=len(v),
                        median_logFC=round(statistics.median(v), 4),
                        frac_up=round(up / len(v), 3), sig_up=su, sig_down=sd,
                        binom_p=pb, competitive_p=pc))
        print('%-30s %4d %+9.4f %5.0f%% %4d/%-4d %10.2e %10.2e'
              % (m[:30], len(v), statistics.median(v), 100 * up / len(v), su, sd, pb, pc))
    return out


rows = []
T = list(csv.DictReader(open(D + '/gandal2018a/TableS1_Microarray_MetaAnalyses.csv',
                             encoding='utf-8')))
for dis, lab in [('ASD', 'Gandal2018a microarray ASD'),
                 ('SCZ', 'Gandal2018a microarray SCZ'),
                 ('BD', 'Gandal2018a microarray BD'),
                 ('MDD', 'Gandal2018a microarray MDD'),
                 ('IBD', 'Gandal2018a IBD (bowel, non-brain control)')]:
    lfc, fdr = {}, {}
    for r in T:
        g = r['hgnc_symbol'] or r['external_gene_id']
        if not g or g in lfc:        # first row per gene symbol
            continue
        lfc[g] = f(r['%s.beta_log2FC' % dis])
        fdr[g] = f(r['%s.FDR' % dis])
    rows += run(lab, lfc, fdr)

with open(D + '/replication_modules.tsv', 'w', newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh, list(rows[0].keys()), delimiter=TAB, lineterminator=NL)
    w.writeheader()
    w.writerows(rows)
print(NL + 'written replication_modules.tsv')

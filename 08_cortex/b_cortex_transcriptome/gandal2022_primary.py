# ASD cortex RNA-seq (Gandal et al. 2022): the 12 panel genes, pathway modules and complement cascade layers (competitive tests).
# Inputs (relative to ASD_ROOT):
#   brain/data/asd_cortex/Gandal2022_MOESM5.xlsx   Gandal et al. 2022 Nature, Supplementary Data 3
#     (sheet DEGene_Statistics: whole-cortex and 11 region-level log2FC and FDR; 49 ASD, 54 control)
#   panel12/data/gene_modules.tsv                  (module_gene_sets.py)
# Outputs (relative to ASD_ROOT):
#   brain/data/asd_cortex/gandal2022_panel12.tsv            12 panel genes, whole cortex and regions up
#   brain/data/asd_cortex/gandal2022_modules.tsv            module tests (incl. COMPLEMENT_all_layers)
#   brain/data/asd_cortex/gandal2022_complement_genes.tsv   complement genes with their layer
# Usage: python gandal2022_primary.py
# Module tests: fraction of member genes with log2FC > 0 (binomial test vs 0.5) and a competitive
# two-sided Wilcoxon rank-sum test of member log2FC against all other genes.
import csv
import math
import os
import statistics
from collections import defaultdict

import openpyxl
from scipy.stats import mannwhitneyu, binomtest

ROOT = os.environ.get("ASD_ROOT", ".")
D = os.path.join(ROOT, 'brain', 'data', 'asd_cortex')
MODULES = os.path.join(ROOT, 'panel12', 'data')
TAB, NL = chr(9), chr(10)
P12 = ['AHSG', 'ANPEP', 'BCHE', 'BTD', 'C1RL', 'C3', 'CLEC3B', 'IGFBP5',
       'MBL2', 'POSTN', 'PTGDS', 'QSOX1']
REGIONS = ['BA9', 'BA44_45', 'BA24', 'BA4_6', 'BA38', 'BA20_37', 'BA41_42_22',
           'BA3_1_2_5', 'BA7', 'BA39_40', 'BA17']


def num(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


wb = openpyxl.load_workbook(D + '/Gandal2022_MOESM5.xlsx', read_only=True)
ws = wb['DEGene_Statistics']
it = ws.iter_rows(values_only=True)
hdr = list(next(it))
ix = {h: i for i, h in enumerate(hdr)}
G = {}
for r in it:
    sym = r[ix['external_gene_name']]
    if not sym or sym in G:          # first row per gene symbol
        continue
    rec = {'wc_lfc': num(r[ix['WholeCortex_ASD_logFC']]),
           'wc_fdr': num(r[ix['WholeCortex_ASD_FDR']])}
    for reg in REGIONS:
        rec[reg + '_lfc'] = num(r[ix['ASD_%s_logFC' % reg]])
        rec[reg + '_fdr'] = num(r[ix['ASD_%s_FDR' % reg]])
    G[sym] = rec
wb.close()
print('genes: %d' % len(G), flush=True)

MOD = defaultdict(list)
for r in csv.DictReader(open(MODULES + '/gene_modules.tsv', encoding='utf-8'), delimiter=TAB):
    MOD[r['module']].append(r['gene'])

# ---- 1. the 12 panel genes ----------------------------------------------------
out1 = []
print(NL + '=== 1. the 12 panel genes, ASD vs control, whole cortex ===', flush=True)
print('%-8s %8s %10s %14s' % ('gene', 'logFC', 'FDR', 'regions up (FDR<.05 up)'))
for g in P12:
    r = G.get(g)
    if not r or r['wc_lfc'] is None:
        print('%-8s  not in dataset' % g)
        out1.append(dict(gene=g, wc_logFC='', wc_FDR='', n_regions_up='',
                         n_regions_sig_up=''))
        continue
    up = sum(1 for reg in REGIONS if (r[reg + '_lfc'] or 0) > 0)
    sig = sum(1 for reg in REGIONS if (r[reg + '_lfc'] or 0) > 0
              and (r[reg + '_fdr'] if r[reg + '_fdr'] is not None else 1) < 0.05)
    print('%-8s %+8.3f %10.2e %8d/11  (%d)' % (g, r['wc_lfc'], r['wc_fdr'], up, sig))
    out1.append(dict(gene=g, wc_logFC=r['wc_lfc'], wc_FDR=r['wc_fdr'],
                     n_regions_up=up, n_regions_sig_up=sig))
with open(D + '/gandal2022_panel12.tsv', 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, list(out1[0].keys()), delimiter=TAB, lineterminator=NL)
    w.writeheader()
    w.writerows(out1)

# ---- 2. modules -------------------------------------------------------------
ALL = {g: r['wc_lfc'] for g, r in G.items() if r['wc_lfc'] is not None}
bg_all = list(ALL.values())
print(NL + 'background: %d genes, median logFC %+.4f, %.1f%% > 0'
      % (len(bg_all), statistics.median(bg_all), 100 * sum(v > 0 for v in bg_all) / len(bg_all)))

rows = []
print(NL + '=== 2. modules, whole cortex ===')
print('%-30s %4s %9s %7s %9s %11s %9s'
      % ('module', 'n', 'med lFC', '% up', 'sig up/dn', 'binom P', 'compet P'))
for m in sorted(MOD):
    v = [ALL[g] for g in MOD[m] if g in ALL]
    genes_in = [g for g in MOD[m] if g in ALL]
    if len(v) < 3:
        continue
    rest = [ALL[g] for g in ALL if g not in set(genes_in)]
    up = sum(x > 0 for x in v)
    su = sum(1 for g in genes_in if ALL[g] > 0 and (G[g]['wc_fdr'] or 1) < 0.05)
    sd = sum(1 for g in genes_in if ALL[g] < 0 and (G[g]['wc_fdr'] or 1) < 0.05)
    pb = binomtest(up, len(v), 0.5).pvalue
    pc = mannwhitneyu(v, rest, alternative='two-sided', method='asymptotic').pvalue
    rows.append(dict(module=m, n=len(v), median_logFC=round(statistics.median(v), 4),
                     mean_logFC=round(statistics.mean(v), 4),
                     frac_up=round(up / len(v), 3), sig_up=su, sig_down=sd,
                     binom_p=pb, competitive_p=pc))
    print('%-30s %4d %+9.4f %6.0f%% %5d/%-3d %11.2e %9.2e'
          % (m, len(v), statistics.median(v), 100 * up / len(v), su, sd, pb, pc))

# all complement layers together
cmp_genes = sorted({g for m in MOD if m.startswith('CMP_') for g in MOD[m] if g in ALL})
v = [ALL[g] for g in cmp_genes]
rest = [ALL[g] for g in ALL if g not in set(cmp_genes)]
up = sum(x > 0 for x in v)
su = sum(1 for g in cmp_genes if ALL[g] > 0 and (G[g]['wc_fdr'] or 1) < 0.05)
sd = sum(1 for g in cmp_genes if ALL[g] < 0 and (G[g]['wc_fdr'] or 1) < 0.05)
pb = binomtest(up, len(v), 0.5).pvalue
pc = mannwhitneyu(v, rest, alternative='two-sided', method='asymptotic').pvalue
rows.append(dict(module='COMPLEMENT_all_layers', n=len(v),
                 median_logFC=round(statistics.median(v), 4),
                 mean_logFC=round(statistics.mean(v), 4),
                 frac_up=round(up / len(v), 3), sig_up=su, sig_down=sd,
                 binom_p=pb, competitive_p=pc))
print('%-30s %4d %+9.4f %6.0f%% %5d/%-3d %11.2e %9.2e'
      % ('COMPLEMENT (all layers)', len(v), statistics.median(v), 100 * up / len(v), su, sd, pb, pc))
with open(D + '/gandal2022_modules.tsv', 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, list(rows[0].keys()), delimiter=TAB, lineterminator=NL)
    w.writeheader()
    w.writerows(rows)

# ---- 3. complement genes, gene by gene --------------------------------------
print(NL + '=== 3. complement genes, whole cortex, sorted by logFC ===')
cg = sorted(cmp_genes, key=lambda g: -ALL[g])
with open(D + '/gandal2022_complement_genes.tsv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, delimiter=TAB, lineterminator=NL)
    w.writerow(['gene', 'layer', 'wc_logFC', 'wc_FDR', 'regions_up', 'regions_sig_up'])
    for g in cg:
        layer = [m for m in MOD if m.startswith('CMP_') and g in MOD[m]][0]
        r = G[g]
        ru = sum(1 for reg in REGIONS if (r[reg + '_lfc'] or 0) > 0)
        rs = sum(1 for reg in REGIONS if (r[reg + '_lfc'] or 0) > 0
                 and (r[reg + '_fdr'] if r[reg + '_fdr'] is not None else 1) < 0.05)
        w.writerow([g, layer, r['wc_lfc'], r['wc_fdr'], ru, rs])
        flag = ' *' if (r['wc_fdr'] or 1) < 0.05 else ''
        print('  %-9s %-28s %+7.3f  FDR %-9.2e  up in %2d/11 regions%s'
              % (g, layer[6:], r['wc_lfc'], r['wc_fdr'], ru, flag))
print(NL + 'done', flush=True)

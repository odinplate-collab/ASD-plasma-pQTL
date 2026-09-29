# Genotype inputs for the diagnostic models of the 12-protein panel:
#   - the lead cis variant of each protein (protein-level, diagnosis-adjusted model, n = 90)
#   - dosages of every variant within each panel gene's span +/- 1 Mb (GRCh37)
# Inputs (relative to ASD_ROOT): pQTL/pQTL_re/ASD_pQTL_imputated_sample90_rm29903_n_alt.tsv,
#   derived/12_cis_summary_all_proteins/protein_level_cis_combined_n90_dxadjusted_ALL.csv.gz,
#   derived/11_protein_maps/measured_proteins_regions_1Mb.tsv
# Outputs (relative to ASD_ROOT): panel12/data/panel12_lead_instruments.tsv, genotype_12leads.tsv.gz,
#   genotype_12regions.tsv.gz
# Usage: python extract_12_regions.py
import csv, gzip, os

ROOT = os.environ.get('ASD_ROOT', '.')
RE   = ROOT + '/pQTL/pQTL_re'
PK   = ROOT + '/derived'
D    = ROOT + '/panel12/data'
GENO = RE + '/ASD_pQTL_imputated_sample90_rm29903_n_alt.tsv'
CIS  = PK + '/12_cis_summary_all_proteins/protein_level_cis_combined_n90_dxadjusted_ALL.csv.gz'
REG  = PK + '/11_protein_maps/measured_proteins_regions_1Mb.tsv'

P12 = ['AHSG','ANPEP','BCHE','BTD','C1RL','C3','CLEC3B','IGFBP5','MBL2','POSTN','PTGDS','QSOX1']

regions = {}
for r in csv.DictReader(open(REG, encoding='utf-8'), delimiter='\t'):
    if r['protein'] in P12:
        regions[r['protein']] = (r['chr'], int(r['lo']), int(r['hi']))
print('regions:', len(regions), flush=True)
assert len(regions) == 12, sorted(set(P12) - set(regions))

# lead cis variant per protein: smallest P in the protein-level combined model (n = 90, diagnosis-adjusted)
best = {}
with gzip.open(CIS, 'rt') as f:
    h = f.readline().rstrip('\n').split(',')
    gi, pi, si = h.index('protein'), h.index('p'), h.index('SNP')
    bi, ei, mi = h.index('beta'), h.index('se'), h.index('maf')
    for l in f:
        a = l.rstrip('\n').split(',')
        if a[gi] in regions:
            p = float(a[pi])
            if a[gi] not in best or p < best[a[gi]][0]:
                best[a[gi]] = (p, a[si], float(a[bi]), float(a[ei]), float(a[mi]))
with open(D + '/panel12_lead_instruments.tsv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['protein','variant_id','beta','se','p','maf','F'])
    for g in P12:
        p, s, b, e, m = best[g]
        w.writerow([g, s, '%.6g' % b, '%.6g' % e, '%.4g' % p, '%.4f' % m, '%.1f' % ((b/e)**2)])
        print('  lead %-9s %-22s F=%.1f MAF=%.3f' % (g, s, (b/e)**2, m), flush=True)

wanted_lead = {best[g][1]: g for g in P12}

with open(GENO) as f:
    hdr = f.readline().rstrip('\n').split('\t')
    cols = [c.replace('-', '_') for c in hdr[1:]]
    out = gzip.open(D + '/genotype_12regions.tsv.gz', 'wt', newline='')
    w = csv.writer(out, delimiter='\t', lineterminator='\n')
    w.writerow(['variant_id'] + cols)
    kept = 0
    lead_rows = {}
    for n, line in enumerate(f, 1):
        vid = line[:line.index('\t')]
        try:
            c, pos = vid.split(':')[0], int(vid.split(':')[1])
        except (ValueError, IndexError):
            continue
        for g, (ch, lo, hi) in regions.items():
            if c == ch and lo <= pos <= hi:
                row = line.rstrip('\n').split('\t')
                w.writerow(row); kept += 1
                if vid in wanted_lead:
                    lead_rows[vid] = row[1:]
                break
        if n % 1000000 == 0:
            print('  scanned %d, kept %d' % (n, kept), flush=True)
    out.close()
print('kept %d variants across the 12 regions' % kept, flush=True)
print('lead variants recovered: %d / 12' % len(lead_rows), flush=True)

with gzip.open(D + '/genotype_12leads.tsv.gz', 'wt', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['variant_id'] + cols)
    for g in P12:
        v = best[g][1]
        if v in lead_rows:
            w.writerow([v] + lead_rows[v])
print('written to', D, flush=True)

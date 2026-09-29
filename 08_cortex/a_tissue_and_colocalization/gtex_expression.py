# Median expression (TPM) of the 12 panel genes across GTEx v8 tissues, via the GTEx Portal API.
# Inputs (relative to ASD_ROOT): none (queries https://gtexportal.org/api/v2)
# Outputs (relative to ASD_ROOT): brain/results/gtex_tpm.tsv
# Usage: python gtex_expression.py
import csv
import json
import os
import urllib.request

ROOT = os.environ.get("ASD_ROOT", ".")
OUT = os.path.join(ROOT, 'brain', 'results')
TAB, NL = chr(9), chr(10)
P12 = ['C3', 'AHSG', 'MBL2', 'C1RL', 'POSTN', 'BCHE', 'BTD', 'ANPEP', 'CLEC3B', 'IGFBP5',
       'PTGDS', 'QSOX1']
API = 'https://gtexportal.org/api/v2'


def write(name, rows):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, list(rows[0].keys()), delimiter=TAB, lineterminator=NL)
        w.writeheader()
        w.writerows(rows)
    print('wrote %-32s %4d rows' % (name, len(rows)))


rows = []
for g in P12:
    # GENCODE v26 gene ID for the symbol, then GTEx v8 median TPM per tissue
    ref = json.load(urllib.request.urlopen(
        API + '/reference/gene?geneId=%s&gencodeVersion=v26&genomeBuild=GRCh38%%2Fhg38' % g, timeout=120))
    gid = [d for d in ref['data'] if d['geneSymbol'] == g][0]['gencodeId']
    med = json.load(urllib.request.urlopen(
        API + '/expression/medianGeneExpression?gencodeId=%s&datasetId=gtex_v8' % gid, timeout=120))
    for d in med['data']:
        rows.append(dict(gene=g, gencode=gid, tissue=d['tissueSiteDetailId'], tpm=d['median']))
write('gtex_tpm.tsv', rows)

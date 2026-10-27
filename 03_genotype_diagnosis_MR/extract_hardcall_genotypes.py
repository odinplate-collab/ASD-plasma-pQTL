# Hard-call genotypes, for every individual with whole-genome sequencing (probands, TD individuals,
# parents and siblings), at the lead cis instrument of Niu et al. for each measured protein
# (the strongest of the instruments selected in niu_mr_colocalization.py), so that variant choice
# is independent of this cohort.
# Genotypes are counts of each alternative allele (0, 1, 2; missing calls left empty), one column
# per variant (chr:pos:ref:alt, GRCh37).
# Inputs (relative to ASD_ROOT):
#   derived/proteome_wide_MR/pw_instruments.csv     (niu_mr_colocalization.py)
#   Variant/GrCh37_nochr/chr<1-22>.vcf.gz           (GRCh37 genotype calls)
# Outputs:
#   argument 1: per-individual genotype table (identifiable data; keep it outside any shared location)
#   derived/genotype_x_diagnosis/niu_lead_instruments.tsv (one instrument per protein; no individual data)
# Usage: python extract_hardcall_genotypes.py <genotype_table.tsv>
import gzip, os, sys, collections
import pandas as pd

ROOT = os.environ.get("ASD_ROOT", ".")
OUT = ROOT + "/derived/genotype_x_diagnosis"; os.makedirs(OUT, exist_ok=True)
DEST = sys.argv[1]

ins = pd.read_csv(ROOT + "/derived/proteome_wide_MR/pw_instruments.csv")
niu = ins[ins.source == "Niu"].sort_values("p", kind="mergesort").drop_duplicates("protein")      # strongest instrument per protein
niu = niu.rename(columns={"protein": "Protein", "key": "chrpos", "rs": "rsID", "beta": "beta_niu", "se": "SE", "p": "P"})
niu[["Protein", "chrpos", "rsID", "ea", "oa", "beta_niu", "SE", "P", "F"]].to_csv(OUT + "/niu_lead_instruments.tsv", sep="\t", index=False)

want = collections.defaultdict(set)              # chromosome -> positions
for k in niu.chrpos:
    c, p = str(k).split(":"); want[c].add(int(p))

geno, samples = {}, None
for chrom in sorted(want, key=lambda x: (len(x), x)):
    pos_set = want[chrom]; last = max(pos_set)
    with gzip.open(f"{ROOT}/Variant/GrCh37_nochr/chr{chrom}.vcf.gz", "rt") as f:
        for line in f:
            if line[0] == "#":
                if line.startswith("#CHROM"): samples = line.rstrip("\n").split("\t")[9:]
                continue
            t1 = line.index("\t"); t2 = line.index("\t", t1 + 1); pos = int(line[t1 + 1:t2])
            if pos > last: break
            if pos not in pos_set: continue
            v = line.rstrip("\n").split("\t")
            for k, alt in enumerate(v[4].split(","), start=1):
                cnt = []
                for x in v[9:]:
                    g = x.split(":")[0].replace("|", "/").split("/")
                    cnt.append("" if "." in g else str(sum(a == str(k) for a in g)))
                geno[f"{chrom}:{pos}:{v[3]}:{alt}"] = cnt
    print(f"chr{chrom}: done", flush=True)
G = pd.DataFrame(geno, index=samples).rename_axis("wgs_id"); G.to_csv(DEST, sep="\t")
found = {":".join(c.split(":")[:2]) for c in G.columns}
print("individuals", len(samples), "| instruments", len(niu), "| positions found", len(found & set(niu.chrpos)))

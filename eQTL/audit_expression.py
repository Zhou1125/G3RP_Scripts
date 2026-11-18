#!/usr/bin/env python3

from pathlib import Path
from collections import Counter, defaultdict
import csv
import re

P = Path(".")
RNA = Path("inputs/rna")
OUT = P / "01_audit"
OUT.mkdir(parents=True, exist_ok=True)

vcf_file = OUT / "vcf_samples.txt"
if not vcf_file.exists():
    raise SystemExit("ERROR: vcf_samples.txt missing")

vcf_samples = [
    x.strip() for x in vcf_file.read_text().splitlines()
    if x.strip()
]
vcf_set = set(vcf_samples)

print("===== 01 VCF =====")
print("VCF samples:", len(vcf_samples))
print("Unique:", len(vcf_set))

files = []
for f in RNA.rglob("*"):
    if not f.is_file():
        continue

    if f.parent.name not in ("05featurecounts", "06featurecounts"):
        continue

    if f.name.endswith(".out_counts.txt"):
        typ = "gene_candidate"
    elif f.name.endswith(".out_counts_exon.txt"):
        typ = "exon_candidate"
    else:
        continue

    m = re.match(r"^(SRR\d+)\.", f.name)
    if not m:
        continue

    files.append({
        "sample": m.group(1),
        "type": typ,
        "path": str(f),
        "group": f.parent.parent.name,
        "size": f.stat().st_size
    })

print("\n===== 02 COUNT FILE INVENTORY =====")

for typ in ("gene_candidate", "exon_candidate"):
    subset = [x for x in files if x["type"] == typ]
    ids = {x["sample"] for x in subset}
    counts = Counter(x["sample"] for x in subset)

    print("\nTYPE:", typ)
    print("FILES:", len(subset))
    print("UNIQUE_SAMPLES:", len(ids))
    print("MATCHED_VCF:", len(ids & vcf_set))
    print("VCF_ONLY:", len(vcf_set - ids))
    print("RNA_ONLY:", len(ids - vcf_set))
    print("DUPLICATE_SAMPLE_FILES:",
          sum(v > 1 for v in counts.values()))

    print("GROUPS:", dict(Counter(x["group"] for x in subset)))

    (OUT / f"{typ}_vcf_only.txt").write_text(
        "\n".join(sorted(vcf_set - ids)) + "\n"
    )
    (OUT / f"{typ}_rna_only.txt").write_text(
        "\n".join(sorted(ids - vcf_set)) + "\n"
    )

with open(OUT / "count_file_manifest.tsv", "w") as handle:
    w = csv.writer(handle, delimiter="\t")
    w.writerow(["SampleID", "Type", "Group", "Size", "Path"])
    for x in sorted(files, key=lambda a: (a["type"], a["sample"])):
        w.writerow([
            x["sample"], x["type"], x["group"],
            x["size"], x["path"]
        ])

print("\n===== 03 COUNT FILE HEADERS =====")

for typ in ("gene_candidate", "exon_candidate"):
    candidates = [x for x in files if x["type"] == typ]
    if not candidates:
        continue

    f = Path(candidates[0]["path"])
    print("\nTYPE:", typ)
    print("FILE:", f)

    with open(f) as handle:
        for _ in range(5):
            line = handle.readline()
            if not line:
                break
            print(line.rstrip()[:500])

print("\n===== 04 FEATURECOUNTS COMMANDS =====")

scripts = sorted(RNA.rglob("01featurecounts.sh"))

for f in scripts:
    if f.parent.name not in ("05featurecounts", "06featurecounts"):
        continue

    print("\nSCRIPT:", f)
    with open(f, errors="replace") as handle:
        for line in handle:
            if (
                "featureCounts" in line
                or "-f " in line
                or "-t " in line
                or "-g " in line
                or "-s " in line
                or "-p " in line
                or "--countReadPairs" in line
            ):
                print(line.rstrip()[:500])

print("\n===== 05 SAMPLE MATCH =====")

gene_ids = {
    x["sample"] for x in files
    if x["type"] == "gene_candidate"
}

exon_ids = {
    x["sample"] for x in files
    if x["type"] == "exon_candidate"
}

with open(OUT / "sample_match.tsv", "w") as handle:
    w = csv.writer(handle, delimiter="\t")
    w.writerow([
        "SampleID", "In_VCF",
        "GeneCount", "ExonCount"
    ])

    for s in sorted(vcf_set | gene_ids | exon_ids):
        w.writerow([
            s,
            int(s in vcf_set),
            int(s in gene_ids),
            int(s in exon_ids)
        ])

print("VCF:", len(vcf_set))
print("GENE:", len(gene_ids))
print("EXON:", len(exon_ids))
print("ALL_THREE:", len(vcf_set & gene_ids & exon_ids))

print("\n===== AUDIT COMPLETE =====")
print("OUTPUT:", OUT)

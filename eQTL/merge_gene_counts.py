#!/usr/bin/env python3

from pathlib import Path
import csv
import numpy as np
import pandas as pd

P = Path(".")
AUDIT = P / "01_audit"
OUT = P / "03_expression"
OUT.mkdir(parents=True, exist_ok=True)

print("===== STEP04: MERGE GENE COUNTS =====", flush=True)

# VCF sample order
samples = [
    x.strip()
    for x in (AUDIT / "vcf_samples.txt").read_text().splitlines()
    if x.strip()
]

assert len(samples) == 58
assert len(set(samples)) == 58

# Read manifest
manifest = pd.read_csv(
    AUDIT / "count_file_manifest.tsv",
    sep="\t",
    dtype=str
)

manifest = manifest[
    manifest["Type"] == "gene_candidate"
].copy()

assert not manifest["SampleID"].duplicated().any()

file_map = dict(
    zip(manifest["SampleID"], manifest["Path"])
)

assert set(file_map) == set(samples), (
    "VCF and RNA samples do not match"
)

# FeatureCounts columns
meta_cols = [
    "Geneid", "Chr", "Start",
    "End", "Strand", "Length"
]

count_series = []
reference_meta = None
sample_qc = []

for index, sample in enumerate(samples, 1):

    path = file_map[sample]

    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        low_memory=False
    )

    if df.shape[1] != 7:
        raise ValueError(
            f"{sample}: expected 7 columns, got {df.shape[1]}"
        )

    if list(df.columns[:6]) != meta_cols:
        raise ValueError(
            f"{sample}: unexpected featureCounts columns"
        )

    if df["Geneid"].isna().any():
        raise ValueError(f"{sample}: missing Geneid")

    if df["Geneid"].duplicated().any():
        raise ValueError(f"{sample}: duplicate Geneid")

    meta = df[meta_cols].copy()

    if reference_meta is None:
        reference_meta = meta.reset_index(drop=True)
    else:
        if not meta.reset_index(drop=True).equals(reference_meta):
            raise ValueError(
                f"{sample}: annotation/order differs from reference"
            )

    counts = pd.to_numeric(
        df.iloc[:, 6],
        errors="raise"
    )

    if counts.isna().any():
        raise ValueError(f"{sample}: missing counts")

    if (counts < 0).any():
        raise ValueError(f"{sample}: negative counts")

    if not np.isfinite(counts.to_numpy()).all():
        raise ValueError(f"{sample}: nonfinite counts")

    if not (counts % 1 == 0).all():
        raise ValueError(f"{sample}: noninteger counts")

    counts = counts.astype(np.int64)
    counts.name = sample

    count_series.append(counts.reset_index(drop=True))

    group = manifest.loc[
        manifest["SampleID"] == sample, "Group"
    ].iloc[0]

    sample_qc.append({
        "SampleID": sample,
        "SourceGroup": group,
        "TotalAssignedCounts": int(counts.sum()),
        "NonzeroGenes": int((counts > 0).sum()),
        "ZeroGenes": int((counts == 0).sum())
    })

    print(
        f"[{index:02d}/{len(samples)}] "
        f"{sample} group={group} "
        f"genes={len(counts)} "
        f"total={int(counts.sum())}",
        flush=True
    )

# Merge matrix
counts_matrix = pd.concat(
    count_series, axis=1
)

counts_matrix.insert(
    0, "Geneid",
    reference_meta["Geneid"]
)

# Construct annotation
annotation = []

for row in reference_meta.itertuples(index=False):

    chromosomes = str(row.Chr).split(";")
    starts = str(row.Start).split(";")
    ends = str(row.End).split(";")
    strands = str(row.Strand).split(";")

    if not (
        len(chromosomes)
        == len(starts)
        == len(ends)
        == len(strands)
    ):
        raise ValueError(
            f"{row.Geneid}: inconsistent exon annotation fields"
        )

    if len(set(chromosomes)) != 1:
        raise ValueError(
            f"{row.Geneid}: exons on multiple chromosomes"
        )

    if len(set(strands)) != 1:
        raise ValueError(
            f"{row.Geneid}: inconsistent strands"
        )

    chrom = chromosomes[0]
    strand = strands[0]

    if strand not in ("+", "-"):
        raise ValueError(
            f"{row.Geneid}: invalid strand {strand}"
        )

    start = min(map(int, starts))
    end = max(map(int, ends))

    if start < 1 or end < start:
        raise ValueError(
            f"{row.Geneid}: invalid coordinates"
        )

    tss = start if strand == "+" else end

    annotation.append({
        "gene_id": row.Geneid,
        "chr": chrom,
        "start": start,
        "end": end,
        "strand": strand,
        "tss": tss,
        "length": int(row.Length)
    })

annotation = pd.DataFrame(annotation)

# Save output
counts_path = OUT / "gene_counts_58.tsv.gz"
annotation_path = OUT / "gene_annotation.tsv"
qc_path = OUT / "sample_count_qc.tsv"
sample_path = OUT / "samples_58.tsv"

counts_matrix.to_csv(
    counts_path,
    sep="\t",
    index=False,
    compression="gzip"
)

annotation.to_csv(
    annotation_path,
    sep="\t",
    index=False
)

qc = pd.DataFrame(sample_qc)
qc.to_csv(qc_path, sep="\t", index=False)

qc[["SampleID", "SourceGroup"]].to_csv(
    sample_path,
    sep="\t",
    index=False
)

# Basic QC
values = counts_matrix.iloc[:, 1:].to_numpy()
nonzero_any = int(np.any(values > 0, axis=1).sum())
expressed_half = int(
    ((values > 0).sum(axis=1) >= 29).sum()
)

print()
print("===== FINAL SUMMARY =====")
print("SAMPLES =", len(samples))
print("GENES =", len(annotation))
print("GENES_NONZERO_ANY =", nonzero_any)
print("GENES_NONZERO_AT_LEAST_29_SAMPLES =", expressed_half)
print("CHROMOSOMES =", annotation["chr"].nunique())
print("MATRIX_SHAPE =", counts_matrix.shape)
print("SOURCE_GROUPS:")
print(qc["SourceGroup"].value_counts().sort_index().to_string())

print()
print("===== OUTPUT FILES =====")
for path in [
    counts_path,
    annotation_path,
    qc_path,
    sample_path
]:
    print(path, path.stat().st_size, "bytes")

print("STEP04_SUCCESS")

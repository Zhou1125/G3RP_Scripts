#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import re

P = Path(".")

RNA = Path(
    "inputs/rna"
)

GTF = Path(
    "inputs/reference/goat.gtf"
)

OUT = P / "03_expression/16_gtex_style"
OUT.mkdir(parents=True, exist_ok=True)

print("===== STEP16A GTEx EXPRESSION FILTER =====", flush=True)

# --------------------------------------------------
# 01. Read counts and sample information
# --------------------------------------------------

samples = pd.read_csv(
    P / "03_expression/samples_58.tsv",
    sep="\t", dtype=str
)

ids = samples["SampleID"].tolist()

counts = pd.read_csv(
    P / "03_expression/gene_counts_58.tsv.gz",
    sep="\t", dtype={"Geneid":str}
)

ann = pd.read_csv(
    P / "03_expression/gene_annotation.tsv",
    sep="\t",
    dtype={"gene_id":str, "chr":str}
)

assert len(ids) == 58
assert len(set(ids)) == 58
assert counts.columns[1:].tolist() == ids
assert counts["Geneid"].tolist() == ann["gene_id"].tolist()
assert counts["Geneid"].is_unique

genes = counts["Geneid"].tolist()
gene_set = set(genes)

C = counts[ids].to_numpy(dtype=np.float64)

assert np.isfinite(C).all()
assert (C >= 0).all()

print("SAMPLES =", len(ids))
print("RAW_GENES =", len(genes))

# --------------------------------------------------
# 02. Merge StringTie TPM
# Audit duplicated Gene IDs instead of summing them
# --------------------------------------------------

T = np.full(C.shape, np.nan)

qc_rows = []
duplicate_rows = []
global_duplicate_ids = set()

required = [
    "Gene ID", "Gene Name", "Reference",
    "Strand", "Start", "End",
    "Coverage", "FPKM", "TPM"
]

for j, sid in enumerate(ids):

    paths = list(
        RNA.glob(
            f"*/04stringtie/{sid}/{sid}.gene_abund.tab"
        )
    )

    if len(paths) != 1:
        raise RuntimeError(
            f"{sid}: expected one TPM file, found {len(paths)}"
        )

    f = paths[0]

    df = pd.read_csv(
        f,
        sep="\t",
        dtype={
            "Gene ID":str,
            "Reference":str,
            "Strand":str
        },
        low_memory=False
    )

    if not set(required).issubset(df.columns):
        raise RuntimeError(
            f"{sid}: missing required columns"
        )

    if df["Gene ID"].isna().any():
        raise RuntimeError(
            f"{sid}: missing Gene ID"
        )

    df["Gene ID"] = df["Gene ID"].str.strip()

    if (df["Gene ID"] == "").any():
        raise RuntimeError(
            f"{sid}: empty Gene ID"
        )

    df["TPM"] = pd.to_numeric(
        df["TPM"], errors="coerce"
    )

    if (
        df["TPM"].isna().any()
        or not np.isfinite(df["TPM"]).all()
        or (df["TPM"] < 0).any()
    ):
        raise RuntimeError(
            f"{sid}: invalid TPM values"
        )

    original_rows = len(df)
    original_tpm_sum = float(df["TPM"].sum())

    # Identify all rows with repeated Gene ID
    dup_mask = df["Gene ID"].duplicated(keep=False)

    dup_df = df.loc[dup_mask, required].copy()
    dup_ids = set(dup_df["Gene ID"])

    if dup_ids:

        global_duplicate_ids.update(dup_ids)

        dup_df.insert(0, "SampleID", sid)
        dup_df.insert(1, "SourceFile", str(f))

        duplicate_rows.append(dup_df)

        # Do not retain arbitrary first row
        # Do not sum ambiguous TPM values
        df = df.loc[~dup_mask].copy()

    assert df["Gene ID"].is_unique

    mapping = pd.Series(
        df["TPM"].to_numpy(dtype=float),
        index=df["Gene ID"]
    )

    T[:,j] = counts["Geneid"].map(mapping).to_numpy(float)

    matched = int(np.isfinite(T[:,j]).sum())

    missing_positive = int(
        (np.isnan(T[:,j]) & (C[:,j] > 0)).sum()
    )

    qc_rows.append({
        "SampleID":sid,
        "SourceGroup":samples.iloc[j]["SourceGroup"],
        "StringTieRows":original_rows,
        "DuplicateGeneIDs":len(dup_ids),
        "DuplicateRows":int(dup_mask.sum()),
        "MatchedCountGenes":matched,
        "MissingTPM":len(genes)-matched,
        "MissingTPM_WithPositiveCounts":missing_positive,
        "OriginalTPMSum":original_tpm_sum
    })

    print(
        f"[{j+1:02d}/58] {sid} "
        f"matched={matched} "
        f"duplicate_gene_ids={len(dup_ids)} "
        f"duplicate_rows={int(dup_mask.sum())}",
        flush=True
    )

# --------------------------------------------------
# 03. Save duplicate-ID audits
# --------------------------------------------------

pd.DataFrame(qc_rows).to_csv(
    OUT / "TPM_sample_audit.tsv",
    sep="\t",
    index=False
)

dup_columns = ["SampleID", "SourceFile"] + required

if duplicate_rows:
    all_dup = pd.concat(
        duplicate_rows,
        ignore_index=True
    )
else:
    all_dup = pd.DataFrame(columns=dup_columns)

all_dup.to_csv(
    OUT / "StringTie_duplicate_gene_rows.tsv",
    sep="\t",
    index=False
)

if len(all_dup):

    summary = (
        all_dup.groupby("Gene ID")
        .agg(
            DuplicateRows=("SampleID","size"),
            AffectedSamples=("SampleID","nunique"),
            DistinctChromosomes=("Reference","nunique"),
            DistinctStrands=("Strand","nunique"),
            MaxTPM=("TPM","max")
        )
        .reset_index()
    )

else:

    summary = pd.DataFrame(columns=[
        "Gene ID", "DuplicateRows", "AffectedSamples",
        "DistinctChromosomes", "DistinctStrands",
        "MaxTPM"
    ])

summary.to_csv(
    OUT / "StringTie_duplicate_gene_summary.tsv",
    sep="\t",
    index=False
)

duplicate_count_ids = global_duplicate_ids & gene_set

print("\n===== STRINGTIE DUPLICATE AUDIT =====")
print("TOTAL_DUPLICATE_IDS =", len(global_duplicate_ids))
print("DUPLICATE_IDS_IN_COUNTS =", len(duplicate_count_ids))
print("DUPLICATE_ROWS =", len(all_dup))

# Safety stop if ambiguity is unexpectedly widespread
if len(duplicate_count_ids) > len(genes) * 0.10:
    raise RuntimeError(
        "More than 10% of count genes have ambiguous "
        "StringTie Gene IDs. Review duplicate audit."
    )

# Save TPM without silently filling missing values with zero
pd.DataFrame(
    T, index=genes, columns=ids
).to_csv(
    OUT / "TPM_58.tsv.gz",
    sep="\t",
    index_label="gene_id",
    compression="gzip",
    na_rep="NA"
)

# --------------------------------------------------
# 04. Extract gene TSS from ARS1.2 GTF
# --------------------------------------------------

gene_pattern = re.compile(r'\bgene_id\s+"([^"]+)"')
bio_pattern = re.compile(r'\bgene_biotype\s+"([^"]+)"')

gtf_records = {}
gtf_duplicates = set()

with GTF.open() as handle:

    for line in handle:

        if line.startswith("#"):
            continue

        z = line.rstrip("\n").split("\t")

        if len(z) < 9 or z[2] != "gene":
            continue

        m = gene_pattern.search(z[8])

        if not m:
            continue

        gid = m.group(1)

        if gid not in gene_set:
            continue

        if gid in gtf_records:
            gtf_duplicates.add(gid)
            continue

        chrom = z[0].removeprefix("chr")
        strand = z[6]

        if strand not in ("+","-"):
            continue

        start = int(z[3])
        end = int(z[4])

        bm = bio_pattern.search(z[8])

        gtf_records[gid] = {
            "gene_id":gid,
            "gtf_chr":chrom,
            "gtf_strand":strand,
            "gtf_tss":start if strand=="+" else end,
            "gene_biotype":(
                bm.group(1) if bm else "unknown"
            )
        }

for gid in gtf_duplicates:
    gtf_records.pop(gid, None)

meta = pd.DataFrame(
    list(gtf_records.values()),
    columns=[
        "gene_id", "gtf_chr", "gtf_strand",
        "gtf_tss", "gene_biotype"
    ]
)

audit = ann.merge(
    meta,
    on="gene_id",
    how="left",
    validate="one_to_one"
)

# --------------------------------------------------
# 05. GTEx V8 expression thresholds
# --------------------------------------------------

min_samples = int(np.ceil(len(ids) * 0.20))

assert min_samples == 12

audit["N_Count_ge6"] = (C >= 6).sum(axis=1)
audit["N_TPM_ge0.1"] = (T >= 0.1).sum(axis=1)
audit["N_TPM_Reported"] = np.isfinite(T).sum(axis=1)

audit["ExpressionPass"] = (
    (audit["N_Count_ge6"] >= min_samples)
    & (audit["N_TPM_ge0.1"] >= min_samples)
)

audit["Autosome"] = audit["chr"].isin(
    {str(i) for i in range(1,30)}
)

audit["AnnotationMatch"] = (
    audit["gtf_chr"].notna()
    & (audit["chr"] == audit["gtf_chr"])
    & (audit["strand"] == audit["gtf_strand"])
    & (audit["tss"] == audit["gtf_tss"])
)

# Globally exclude ambiguous StringTie IDs
audit["StringTieDuplicateID"] = (
    audit["gene_id"].isin(global_duplicate_ids)
)

audit["FinalSelected"] = (
    audit["ExpressionPass"]
    & audit["Autosome"]
    & audit["AnnotationMatch"]
    & ~audit["StringTieDuplicateID"]
)

audit.to_csv(
    OUT / "GTEx_filter_audit.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip"
)

selected = audit.loc[
    audit["FinalSelected"],
    ["gene_id","gtf_chr","gtf_tss","gene_biotype"]
].copy()

selected["chr_num"] = selected["gtf_chr"].astype(int)

selected = selected.sort_values(
    ["chr_num","gtf_tss","gene_id"]
).drop(columns="chr_num")

if selected.empty:
    raise RuntimeError(
        "No genes passed GTEx filtering"
    )

selected.to_csv(
    OUT / "GTEx_selected_genes.tsv",
    sep="\t",
    index=False
)

print("\n===== FINAL SUMMARY =====")
print("SAMPLES =",len(ids))
print("RAW_GENES =",len(genes))
print("MIN_SAMPLES =",min_samples)
print("EXPRESSION_PASS =",int(audit["ExpressionPass"].sum()))
print("AUTOSOMAL_EXPRESSION_PASS =",int(
    (audit["ExpressionPass"] & audit["Autosome"]).sum()
))
print("DUPLICATE_IDS_IN_COUNTS =",len(duplicate_count_ids))
print("DUPLICATE_IDS_EXCLUDED_FROM_ELIGIBLE =",int(
    (
        audit["ExpressionPass"]
        & audit["Autosome"]
        & audit["AnnotationMatch"]
        & audit["StringTieDuplicateID"]
    ).sum()
))
print("FINAL_SELECTED_GENES =",len(selected))
print("STEP16A_SUCCESS",flush=True)

#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd

P = Path(".")
D = P / "05_fastqtl/22_fdr"
OUT = P / "05_fastqtl/23_nominal"
OUT.mkdir(parents=True, exist_ok=True)

models = {}

for name in ["M1", "M1_PEER5"]:
    f = D / name / "all_genes_FDR.tsv.gz"
    df = pd.read_csv(f, sep="\t")

    assert len(df) == 20093
    assert df["gene_id"].is_unique

    models[name] = df.set_index("gene_id")

a = models["M1"]
b = models["M1_PEER5"]

assert set(a.index) == set(b.index)

sig_a = a["qval"].le(0.05)
sig_b = b["qval"].le(0.05)

A = set(a.index[sig_a])
B = set(b.index[sig_b])

summary = {
    "M1_EGENES": len(A),
    "M1_PEER5_EGENES": len(B),
    "SHARED_EGENES": len(A & B),
    "M1_ONLY": len(A - B),
    "PEER5_ONLY": len(B - A),
    "UNION_EGENES": len(A | B)
}

print("===== STEP23A eGENE COMPARISON =====")

for k, v in summary.items():
    print(k, "=", v)

comparison = pd.DataFrame({
    "gene_id": sorted(A | B)
})

for name, df in models.items():
    for col in [
        "qval", "pval_beta", "variant_id",
        "pval_nominal", "slope",
        "pval_nominal_threshold"
    ]:
        comparison[f"{col}_{name}"] = (
            comparison["gene_id"].map(df[col])
        )

comparison["Significant_M1"] = comparison["gene_id"].isin(A)
comparison["Significant_PEER5"] = comparison["gene_id"].isin(B)

comparison["Shared_eGene"] = (
    comparison["Significant_M1"]
    & comparison["Significant_PEER5"]
)

comparison["Same_lead_variant"] = (
    comparison["Shared_eGene"]
    & (
        comparison["variant_id_M1"]
        == comparison["variant_id_M1_PEER5"]
    )
)

comparison.to_csv(
    OUT / "eGene_model_comparison.tsv",
    sep="\t", index=False
)

for name, genes in [
    ("shared", A & B),
    ("M1_only", A - B),
    ("PEER5_only", B - A),
    ("union", A | B)
]:
    (OUT / f"{name}_egenes.txt").write_text(
        "\n".join(sorted(genes)) + "\n"
    )

print("SHARED_SAME_LEAD =", int(
    comparison["Same_lead_variant"].sum()
))

# Data-driven threshold: must exceed each significant
# eGene's gene-specific nominal cutoff.
rows = []

for name, df in models.items():
    selected = df[df["qval"] <= 0.05]

    t = pd.to_numeric(
        selected["pval_nominal_threshold"],
        errors="coerce"
    )

    if t.isna().any() or not np.isfinite(t).all():
        raise RuntimeError(f"{name}: invalid nominal thresholds")

    if ((t < 0) | (t > 1)).any():
        raise RuntimeError(f"{name}: thresholds outside [0,1]")

    maximum = float(t.max())

    # Retain nominal associations beyond gene-specific cutoffs
    planned = 1e-4

    if not np.isfinite(maximum) or maximum >= planned:
        raise RuntimeError(
            f"{name}: nominal cap {planned} does not exceed "
            f"all gene-specific thresholds ({maximum}); review before running"
        )

    rows.append({
        "Model": name,
        "Max_eGene_nominal_threshold": maximum,
        "Planned_nominal_threshold": planned
    })

    print(name, "MAX_GENE_THRESHOLD =", maximum)
    print(name, "NOMINAL_OUTPUT_THRESHOLD =", planned)

pd.DataFrame(rows).to_csv(
    OUT / "nominal_thresholds.tsv",
    sep="\t", index=False
)

print("STEP23A_SUCCESS")

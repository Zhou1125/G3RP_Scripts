#!/usr/bin/env python3

from pathlib import Path
import gzip
import numpy as np
import pandas as pd

P = Path(
    "."
)

OUT = P / "04_covariates/17_GTEx"
OUT.mkdir(parents=True, exist_ok=True)

BED = P / (
    "03_expression/16_gtex_style/"
    "GTEx_style_TMM_INT_58.bed.gz"
)

PCA = P / (
    "02_genotype/13_ld_pca/"
    "skin58_pca.eigenvec"
)

PSAM = P / (
    "02_genotype/14_dosage/"
    "skin58_DSfix.psam"
)

print("===== STEP17 GTEx COVARIATES =====", flush=True)

# 01. VCF-matched sample information
samples = pd.read_csv(
    P / "03_expression/samples_58.tsv",
    sep="\t",
    dtype=str
)

ids = samples["SampleID"].tolist()

assert len(ids) == 58
assert len(set(ids)) == 58

# 02. Verify official phenotype BED
with gzip.open(BED, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
    n_genes = sum(1 for line in f if line.strip())

assert header[:4] == [
    "#chr", "start", "end", "phenotype_id"
]
assert header[4:] == ids
assert n_genes == 20093

print("BED_GENES =", n_genes)
print("BED_SAMPLES =", len(ids))

# 03. Verify dosage PGEN sample set
psam = pd.read_csv(
    PSAM,
    sep=r"\s+",
    dtype=str
)

iid = "#IID" if "#IID" in psam.columns else "IID"

assert iid in psam.columns
assert len(psam) == 58
assert set(psam[iid]) == set(ids)

print("DOSAGE_PGEN_SAMPLE_MATCH = TRUE")

# 04. Load LD-pruned genotype PCs
pca = pd.read_csv(
    PCA,
    sep=r"\s+",
    dtype=str
)

assert "IID" in pca.columns
assert pca["IID"].is_unique
assert set(pca["IID"]) == set(ids)

pca = pca.set_index("IID").loc[ids]

pc_cols = [f"PC{i}" for i in range(1, 6)]

pcs = pca[pc_cols].apply(
    pd.to_numeric,
    errors="raise"
)

assert np.isfinite(pcs.to_numpy()).all()

# Standardize PCs
pc_sd = pcs.std(axis=0, ddof=0)

assert (pc_sd > 0).all()

pcs = (pcs - pcs.mean()) / pc_sd

pcs.columns = [
    f"Genotype_PC{i}" for i in range(1, 6)
]

# 05. SourceGroup dummy coding
groups = samples["SourceGroup"].astype(str)

levels = sorted(groups.unique(), key=int)

print("SOURCE_GROUPS =", levels)
print("SOURCE_GROUP_COUNTS:")
print(groups.value_counts().sort_index().to_string())

cat = pd.Categorical(
    groups,
    categories=levels
)

dummy = pd.get_dummies(
    cat,
    drop_first=True,
    dtype=float
)

dummy.index = ids
dummy.columns = [
    f"SourceGroup_{x}" for x in levels[1:]
]

pcs.index = ids

assert dummy.shape[1] == 6

# 06. Models
models = {
    "M0_GenoPC3": pcs.iloc[:, :3],
    "M0_GenoPC5": pcs.iloc[:, :5],
    "M1_Source_GenoPC3": pd.concat(
        [dummy, pcs.iloc[:, :3]], axis=1
    ),
    "M1_Source_GenoPC5": pd.concat(
        [dummy, pcs.iloc[:, :5]], axis=1
    )
}

summary = []

for name, df in models.items():

    assert df.index.tolist() == ids

    X = df.to_numpy(dtype=float)

    assert np.isfinite(X).all()

    design = np.column_stack(
        [np.ones(58), X]
    )

    rank = np.linalg.matrix_rank(design)
    condition = np.linalg.cond(design)

    if rank != design.shape[1]:
        raise RuntimeError(
            f"{name}: rank deficient"
        )

    if condition > 1e6:
        raise RuntimeError(
            f"{name}: excessive condition number"
        )

    residual_df = 58 - rank - 1

    if residual_df <= 0:
        raise RuntimeError(
            f"{name}: insufficient regression DF"
        )

    # FastQTL / GTEx format:
    # covariates x samples
    path = OUT / f"{name}.covariates.tsv"

    df.T.to_csv(
        path,
        sep="\t",
        index=True,
        index_label="id",
        float_format="%.10g"
    )

    summary.append({
        "Model": name,
        "N_samples": 58,
        "N_covariates": X.shape[1],
        "Design_rank": rank,
        "Condition_number": condition,
        "Residual_DF_with_SNP": residual_df,
        "Path": str(path)
    })

    print(
        name,
        "COV=", X.shape[1],
        "RANK=", rank,
        "CONDITION=", round(condition, 4),
        "REGRESSION_DF=", residual_df,
        flush=True
    )

# 07. Prepare observed covariates for PEER
# Rows are samples, columns are known covariates.
known = models["M1_Source_GenoPC3"].copy()

known.index.name = "SampleID"

known.to_csv(
    OUT / "PEER_observed_covariates.tsv",
    sep="\t",
    float_format="%.10g"
)

# 08. Source group association with genotype PCs
source = groups.to_numpy()

for col in pcs.columns:

    y = pcs[col].to_numpy()
    mean = y.mean()

    ss_total = ((y - mean)**2).sum()
    ss_between = 0.0

    for g in levels:
        v = y[source == g]
        ss_between += len(v) * (v.mean()-mean)**2

    eta2 = ss_between / ss_total

    print(
        col,
        "SOURCE_GROUP_ETA2=",
        round(eta2, 5)
    )

pd.DataFrame(summary).to_csv(
    OUT / "covariate_model_audit.tsv",
    sep="\t",
    index=False
)

print("\n===== FINAL SUMMARY =====")
print("SAMPLES =", len(ids))
print("EXPRESSION_GENES =", n_genes)
print("GENOTYPE_PC_SOURCE = LD_PRUNED_227239")
print("SOURCE_GROUPS =", len(levels))
print("MODELS =", len(models))

print("STEP17_SUCCESS", flush=True)

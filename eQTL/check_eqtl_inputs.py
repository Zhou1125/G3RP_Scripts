#!/usr/bin/env python3

from pathlib import Path
import gzip
import numpy as np
import pandas as pd

P = Path(".")

BED = P / "03_expression/16_gtex_style/GTEx_style_TMM_INT_58.bed.gz"
COVDIR = P / "04_covariates/17_GTEx"
PEERDIR = P / "04_covariates/18_PEER/PEER5"

samples = pd.read_csv(
    P / "03_expression/samples_58.tsv",
    sep="\t", dtype=str
)

ids = samples["SampleID"].tolist()

assert len(ids) == 58
assert len(set(ids)) == 58

print("===== STEP19 INPUT AUDIT =====", flush=True)

# Phenotype BED
with gzip.open(BED, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
    genes = sum(1 for line in f if line.strip())

assert header[:4] == [
    "#chr", "start", "end", "phenotype_id"
]
assert header[4:] == ids
assert genes == 20093
assert Path(str(BED)+".tbi").stat().st_size > 0

print("BED_GENES =", genes)
print("BED_SAMPLES =", len(ids))

files = {
    "M1": COVDIR / "M1_Source_GenoPC3.covariates.tsv",
    "M1_PEER5": PEERDIR / "M1_PEER5.covariates.tsv"
}

for name, path in files.items():
    cov = pd.read_csv(path, sep="\t", index_col=0)

    assert cov.columns.tolist() == ids
    assert cov.index.is_unique

    X = cov.to_numpy(dtype=float)
    assert np.isfinite(X).all()

    D = np.column_stack([np.ones(58), X.T])
    rank = np.linalg.matrix_rank(D)
    condition = np.linalg.cond(D)
    df = 58-rank-1

    print(
        name,
        "COVARIATES =", X.shape[0],
        "RANK =", rank,
        "CONDITION =", round(condition,4),
        "SNP_REGRESSION_DF =", df
    )

    assert rank == D.shape[1]
    assert condition < 1e6
    assert df >= 20

peer_file = PEERDIR / "PEER5_factors.tsv"
peer = pd.read_csv(peer_file, sep="\t", index_col=0)

assert peer.shape == (5,58)
assert peer.columns.tolist() == ids
assert np.isfinite(peer.to_numpy(dtype=float)).all()

print("PEER5_FACTORS = 5")
print("PEER5_SAMPLE_ORDER = MATCHED")
print("STEP19_INPUT_SUCCESS")

#!/usr/bin/env python3

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import subprocess
from pathlib import Path

P = Path(".")

BIN = P / "tools/bin/fastQTL"

VCF = Path(
    "inputs/genotypes/skin58.vcf.gz"
)

BED = P / "03_expression/16_gtex_style/GTEx_style_TMM_INT_58.bed.gz"

COVS = {
    "M1": P / "04_covariates/17_GTEx/M1_Source_GenoPC3.covariates.tsv",
    "M1_PEER5": P / "04_covariates/18_PEER/PEER5/M1_PEER5.covariates.tsv"
}

HEADERS = {
    "permute": (
        "gene_id num_var beta_shape1 beta_shape2 true_df pval_true_df "
        "variant_id tss_distance ma_samples ma_count maf ref_factor "
        "pval_nominal slope slope_se pval_perm pval_beta"
    ).split(),
    "nominal": (
        "gene_id variant_id tss_distance ma_samples ma_count "
        "maf pval_nominal slope slope_se"
    ).split()
}

parser = argparse.ArgumentParser()
parser.add_argument("--model", choices=sorted(COVS), required=True)
parser.add_argument("--mode", choices=["permute", "nominal"], required=True)
parser.add_argument("--chunks", type=int, default=100)
parser.add_argument("--workers", type=int, default=4)
parser.add_argument("--nominal-threshold", type=float, default=1e-4)
args = parser.parse_args()

if args.chunks < 29 or args.workers < 1:
    parser.error("chunks must be >=29 and workers >=1")
if not 0 < args.nominal_threshold <= 1:
    parser.error("nominal threshold must be in (0, 1]")

OUT = P / "05_fastqtl/21_full" / args.model / args.mode
OUT.mkdir(parents=True, exist_ok=True)

for f in [
    BIN, VCF, Path(str(VCF)+".tbi"),
    BED, Path(str(BED)+".tbi"), COVS[args.model]
]:
    if not f.is_file() or f.stat().st_size == 0:
        raise RuntimeError("Missing input: " + str(f))

expected_cols = len(HEADERS[args.mode])

def fingerprint(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

# Reuse chunks only when inputs and analysis settings are unchanged.
signature = {
    'model': args.model, 'mode': args.mode, 'chunks': args.chunks,
    'nominal_threshold': args.nominal_threshold if args.mode == 'nominal' else None,
    'window': 1000000, 'maf': 0.05, 'ma_samples': 3,
    'permutations': [1000, 10000], 'seed_base': 20261031,
    'inputs': {str(f): fingerprint(f) for f in
               [BIN, VCF, Path(str(VCF)+'.tbi'), BED,
                Path(str(BED)+'.tbi'), COVS[args.model]]}
}
manifest = OUT / 'run_signature.json'
if manifest.exists():
    if json.loads(manifest.read_text()) != signature:
        raise RuntimeError('Inputs/settings changed; use a fresh output directory.')
else:
    if list(OUT.glob('chunk_*')):
        raise RuntimeError('Existing chunks have no input signature; use a fresh output directory.')
    manifest.write_text(json.dumps(signature, indent=2) + '\n')

def verify_gzip(path):
    if not path.is_file() or path.stat().st_size == 0:
        return False

    cp = subprocess.run(
        ["gzip", "-t", str(path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )

    if cp.returncode:
        return False

    with gzip.open(path, "rt") as fh:
        first = fh.readline().rstrip("\n")

    return (args.mode == "nominal" and first == "") or (bool(first) and len(first.split("\t")) == expected_cols)

def run_chunk(i):
    prefix = OUT / f"chunk_{i:03d}"

    result = Path(str(prefix)+".txt.gz")
    ok = Path(str(prefix)+".ok")
    log = Path(str(prefix)+".fastqtl.log")
    captured = Path(str(prefix)+".runner.log")

    if ok.is_file() and verify_gzip(result):
        return i, "REUSED"

    ok.unlink(missing_ok=True)

    cmd = [
        str(BIN),
        "--vcf", str(VCF),
        "--bed", str(BED),
        "--cov", str(COVS[args.model]),
        "--window", "1000000",
        "--maf-threshold", "0.05",
        "--ma-sample-threshold", "3",
        "--chunk", str(i), str(args.chunks),
        "--seed", str(20261031+i),
        "--log", str(log),
        "--out", str(result)
    ]

    if args.mode == "permute":
        cmd += ["--permute", "1000", "10000"]
    else:
        cmd += ["--threshold", str(args.nominal_threshold)]

    with captured.open("w") as fh:
        cp = subprocess.run(
            cmd, stdout=fh, stderr=subprocess.STDOUT
        )

    if cp.returncode != 0 or not verify_gzip(result):
        raise RuntimeError(
            f"Chunk {i} failed. "
            f"Exit={cp.returncode}. "
            f"See {captured} and {log}"
        )

    ok.write_text("OK\n")
    return i, "FINISHED"

print(
    f"MODEL={args.model} MODE={args.mode} "
    f"CHUNKS={args.chunks} WORKERS={args.workers}",
    flush=True
)

failures = []

with concurrent.futures.ThreadPoolExecutor(
    max_workers=args.workers
) as executor:

    futures = [
        executor.submit(run_chunk, i)
        for i in range(1, args.chunks+1)
    ]

    for future in concurrent.futures.as_completed(futures):
        try:
            i, status = future.result()
            print(f"CHUNK={i:03d} {status}", flush=True)
        except Exception as exc:
            failures.append(str(exc))
            print("ERROR:", exc, flush=True)

if failures:
    raise RuntimeError(
        f"{len(failures)} chunks failed. "
        "Review individual chunk logs before resubmitting."
    )

# Merge outputs without loading full results into memory.
merged = OUT / f"{args.model}.{args.mode}.tsv.gz"
temporary = OUT / f"{args.model}.{args.mode}.tsv.gz.partial"

rows = 0
seen_genes = set()

with gzip.open(temporary, "wt", compresslevel=1) as dst:
    dst.write("\t".join(HEADERS[args.mode]) + "\n")

    for i in range(1, args.chunks+1):
        f = OUT / f"chunk_{i:03d}.txt.gz"

        with gzip.open(f, "rt") as src:
            for line in src:
                if not line.strip():
                    continue

                if len(line.rstrip("\n").split("\t")) != expected_cols:
                    raise RuntimeError(f"Invalid columns in chunk {i}")

                if args.mode == "permute":
                    gid = line.split("\t", 1)[0]

                    if gid in seen_genes:
                        raise RuntimeError(
                            f"Duplicated permutation gene: {gid}"
                        )

                    seen_genes.add(gid)

                dst.write(line)
                rows += 1

if rows < 1 or not verify_gzip(temporary):
    raise RuntimeError("Merged output invalid")

temporary.replace(merged)

print("MERGED_ROWS =", rows, flush=True)
print("OUTPUT =", merged, flush=True)
print(
    f"STEP21_{args.model}_{args.mode.upper()}_SUCCESS",
    flush=True
)

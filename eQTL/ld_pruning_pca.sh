#!/usr/bin/env bash
set -euo pipefail

P=.
OUT="$P/02_genotype/13_ld_pca"
mkdir -p "$OUT"

PLINK2=./tools/bin/plink2

VCF=inputs/genotypes/skin58.vcf.gz

echo "===== SOFTWARE ====="
"$PLINK2" --version

echo "===== STEP13A: IMPORT AUTOSOMAL VCF ====="

"$PLINK2" \
  --vcf "$VCF" \
  --chr-set 29 no-xy \
  --allow-extra-chr \
  --chr 1-29 \
  --double-id \
  --make-pgen \
  --threads 4 \
  --out "$OUT/skin58_autosomes"

echo "===== STEP13B: GENOTYPE QC ====="

"$PLINK2" \
  --pfile "$OUT/skin58_autosomes" \
  --chr-set 29 no-xy \
  --snps-only just-acgt \
  --max-alleles 2 \
  --maf 0.05 \
  --geno 0.05 \
  --make-pgen \
  --threads 4 \
  --out "$OUT/skin58_qc"

echo "===== STEP13C: LD PRUNING ====="

"$PLINK2" \
  --pfile "$OUT/skin58_qc" \
  --chr-set 29 no-xy \
  --indep-pairwise 200 50 0.2 \
  --threads 4 \
  --out "$OUT/skin58_ld"

echo "===== STEP13D: PCA ====="

"$PLINK2" \
  --pfile "$OUT/skin58_qc" \
  --chr-set 29 no-xy \
  --extract "$OUT/skin58_ld.prune.in" \
  --pca 10 \
  --threads 4 \
  --out "$OUT/skin58_pca"

echo "===== STEP13 SUMMARY ====="

echo "LD SNPs:"
wc -l "$OUT/skin58_ld.prune.in"

echo "PCA samples:"
awk 'NR>1{n++} END{print n+0}' "$OUT/skin58_pca.eigenvec"

echo "Eigenvalues:"
head -10 "$OUT/skin58_pca.eigenval"

echo "STEP13_SUCCESS"

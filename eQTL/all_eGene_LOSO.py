#!/usr/bin/env python3
import csv
import gzip
import math
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t

P = Path('.')

VCF = Path(
    'inputs/genotypes/skin58.vcf.gz'
)

BED = P / '03_expression/16_gtex_style/GTEx_style_TMM_INT_58.bed.gz'

TABIX = Path(
    './tools/bin/tabix'
)

OUT = P / '05_fastqtl/25B_LOSO'
OUT.mkdir(parents=True, exist_ok=True)

COV_FILES = {
    'M1': P / '04_covariates/17_GTEx/M1_Source_GenoPC3.covariates.tsv',
    'M1_PEER5': P / '04_covariates/18_PEER/PEER5/M1_PEER5.covariates.tsv'
}

def read_tsv(path):
    with path.open(newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def covariates(path, destination_samples):
    with path.open(newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        original_samples = header[1:]
        rows = list(reader)

    if len(original_samples) != 58:
        raise RuntimeError('Covariate sample count != 58')

    if set(original_samples) != set(destination_samples):
        raise RuntimeError('Covariate sample mismatch')

    if any(len(r) != len(header) for r in rows):
        raise RuntimeError('Malformed covariate matrix')

    names = [r[0] for r in rows]

    matrix = np.array(
        [[float(v) for v in r[1:]] for r in rows],
        dtype=float
    ).T

    order = [
        original_samples.index(s)
        for s in destination_samples
    ]

    matrix = matrix[order, :]

    if not np.isfinite(matrix).all():
        raise RuntimeError('Nonfinite covariates')

    return names, matrix

def fit(y, x, cov):
    good = (
        np.isfinite(x) &
        np.isfinite(y) &
        np.isfinite(cov).all(axis=1)
    )

    x = x[good]
    y = y[good]
    cov = cov[good]

    n = len(x)

    result = dict(
        n=n, rank=0, df=math.nan,
        slope=math.nan, se=math.nan,
        p=math.nan,
        ci95_low=math.nan,
        ci95_high=math.nan,
        status='NOT_ESTIMABLE'
    )

    if n < 15:
        result['status'] = 'N_LT_15'
        return result

    C = np.column_stack((np.ones(n), cov))
    rank = int(np.linalg.matrix_rank(C))
    df = n - rank - 1

    result.update(rank=rank, df=df)

    if df < 5:
        result['status'] = 'DF_LT_5'
        return result

    # Frisch-Waugh-Lovell residualization
    xc = x - C @ np.linalg.lstsq(C, x, rcond=None)[0]
    yc = y - C @ np.linalg.lstsq(C, y, rcond=None)[0]

    ss = float(xc @ xc)

    if ss <= 1e-9:
        result['status'] = 'NO_RESIDUAL_GENOTYPE_VARIATION'
        return result

    slope = float((xc @ yc) / ss)
    residual = yc - slope * xc
    rss = float(residual @ residual)

    se = math.sqrt(max(0.0, rss / df / ss))

    if se > 0:
        p = float(2 * student_t.sf(abs(slope / se), df))
    else:
        p = 0.0 if slope != 0 else 1.0

    ci = float(student_t.ppf(0.975, df)) * se

    result.update(
        slope=slope,
        se=se,
        p=p,
        ci95_low=slope-ci,
        ci95_high=slope+ci,
        status='OK'
    )

    return result

def write_table(name, rows):
    target = OUT / name
    temporary = Path(str(target) + '.partial')

    with temporary.open('w', newline='') as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=list(rows[0]),
            delimiter='\t',
            lineterminator='\n'
        )
        writer.writeheader()
        writer.writerows(rows)

    temporary.replace(target)
    print('OUTPUT =', target, flush=True)

print('===== STEP25B ALL eGENE LEAD LOSO =====', flush=True)

for f in [VCF, Path(str(VCF)+'.tbi'), BED, TABIX] + list(COV_FILES.values()):
    if not f.is_file() or not f.stat().st_size:
        raise RuntimeError(f'MISSING_INPUT: {f}')

# Load all significant eGenes and their original FastQTL results.
models = {}

for model, expected in [('M1', 175), ('M1_PEER5', 272)]:
    summary = read_tsv(
        P / f'05_fastqtl/24_significant/{model}.egene_pair_summary.tsv'
    )

    if len(summary) != expected:
        raise RuntimeError(f'{model}: unexpected eGene count')

    if len({r['gene_id'] for r in summary}) != expected:
        raise RuntimeError(f'{model}: duplicate genes')

    permutation = {}

    with gzip.open(
        P / f'05_fastqtl/21_full/{model}/permute/{model}.permute.tsv.gz',
        'rt'
    ) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            permutation[row['gene_id']] = row

    for row in summary:
        gene = row['gene_id']

        if gene not in permutation:
            raise RuntimeError(f'Missing permutation: {gene}')

        if permutation[gene]['variant_id'] != row['permutation_lead_variant']:
            raise RuntimeError(f'Lead SNP mismatch: {gene}')

        if float(row['qval']) > 0.05:
            raise RuntimeError(f'Non-significant gene: {gene}')

    models[model] = (summary, permutation)

genes = {
    r['gene_id']
    for summary, _ in models.values()
    for r in summary
}

variants = {
    r['permutation_lead_variant']
    for summary, _ in models.values()
    for r in summary
}

if len(genes) != 326:
    raise RuntimeError('Expected 326 distinct eGenes')

print('MODEL_GENE_LEADS = 447')
print('UNIQUE_EGENES =', len(genes))
print('UNIQUE_LEAD_SNPS =', len(variants))

# Read selected expression phenotypes once.
phenotypes = {}

with gzip.open(BED, 'rt') as fh:
    header = next(fh).rstrip('\n').split('\t')

    if len(header) != 62:
        raise RuntimeError('Expected BED with 58 samples')

    samples = header[4:]

    if len(set(samples)) != 58:
        raise RuntimeError('Duplicate BED sample IDs')

    for line in fh:
        fields = line.rstrip('\n').split('\t')

        if len(fields) != 62:
            raise RuntimeError('Malformed BED row')

        gene = fields[3]

        if gene in genes:
            if gene in phenotypes:
                raise RuntimeError(f'Duplicate phenotype: {gene}')

            phenotypes[gene] = np.array(
                [float(v) for v in fields[4:]],
                dtype=float
            )

if set(phenotypes) != genes:
    raise RuntimeError('Missing eGene expression phenotypes')

print('BED_EXPRESSION_MATCH = 326', flush=True)

# Read and align both covariate matrices.
cov = {}

for model, file in COV_FILES.items():
    names, matrix = covariates(file, samples)
    cov[model] = matrix

    rank = np.linalg.matrix_rank(
        np.column_stack((np.ones(58), matrix))
    )

    print(
        f'{model} N_COV={len(names)} DESIGN_RANK={rank}',
        flush=True
    )

# Reconstruct source groups using the six binary indicators in M1.
names, _ = covariates(COV_FILES['M1'], samples)
matrix = cov['M1']

indicators = []

for i, name in enumerate(names):
    x = matrix[:, i]

    binary = np.all(
        (np.abs(x) < 1e-8) |
        (np.abs(x-1) < 1e-8)
    )

    if binary and 0 < x.sum() < 58:
        indicators.append((name, x))

if len(indicators) != 6:
    raise RuntimeError('Expected six source group indicators')

S = np.column_stack([x for _, x in indicators])

if np.any(S.sum(axis=1) > 1 + 1e-7):
    raise RuntimeError('Overlapping source groups')

groups = np.array([
    next(
        (indicators[j][0] for j in range(6) if S[i,j] > 0.5),
        'Source_baseline'
    )
    for i in range(58)
], dtype=object)

group_sizes = Counter(groups)

if sorted(group_sizes.values()) != [6,6,6,6,7,11,16]:
    raise RuntimeError(f'Source group mismatch: {dict(group_sizes)}')

print('SOURCE_GROUP_SIZES =', dict(sorted(group_sizes.items())))

# Match the original VCF sample order.
vcf_header = subprocess.run(
    [str(TABIX), '-H', str(VCF)],
    capture_output=True,
    text=True,
    check=True
).stdout

chrom_headers = [
    line for line in vcf_header.splitlines()
    if line.startswith('#CHROM\t')
]

if len(chrom_headers) != 1:
    raise RuntimeError('VCF header invalid')

vcf_samples = chrom_headers[0].split('\t')[9:]

if len(vcf_samples) != 58 or set(vcf_samples) != set(samples):
    raise RuntimeError('VCF/BED sample mismatch')

indices = [vcf_samples.index(s) for s in samples]

print('VCF_BED_COV_SAMPLES_MATCH = 58')

# Extract all unique lead SNP positions.
positions = set()

for variant in variants:
    parts = variant.split('_', 3)

    if len(parts) != 4 or not parts[1].isdigit():
        raise RuntimeError(f'Invalid variant ID: {variant}')

    positions.add((parts[0], int(parts[1])))

region_file = OUT / 'lead_positions_1based.tsv'

with region_file.open('w') as fh:
    for chrom, pos in sorted(
        positions,
        key=lambda z: (
            int(z[0]) if z[0].isdigit() else 999,
            z[0], z[1]
        )
    ):
        fh.write(f'{chrom}\t{pos}\n')

vcf_lines = subprocess.run(
    [str(TABIX), '-R', str(region_file), str(VCF)],
    capture_output=True,
    text=True,
    check=True
).stdout

dosages = {}

for line in vcf_lines.splitlines():
    fields = line.split('\t')

    variant = '_'.join([
        fields[0], fields[1], fields[3], fields[4]
    ])

    if variant not in variants:
        continue

    if variant in dosages:
        raise RuntimeError(f'Duplicate VCF variant: {variant}')

    fmt = fields[8].split(':')

    if 'DS' not in fmt:
        raise RuntimeError(f'Missing DS: {variant}')

    ds_index = fmt.index('DS')
    values = []

    for i in indices:
        parts = fields[9+i].split(':')
        raw = parts[ds_index] if ds_index < len(parts) else '.'

        if ',' in raw:
            raise RuntimeError(f'Multiallelic DS: {variant}')

        values.append(
            float(raw) if raw not in ('.','NA','') else math.nan
        )

    x = np.asarray(values, dtype=float)

    if np.any(
        (x[np.isfinite(x)] < -1e-4) |
        (x[np.isfinite(x)] > 2.0001)
    ):
        raise RuntimeError(f'Invalid DS: {variant}')

    dosages[variant] = x

if set(dosages) != variants:
    raise RuntimeError('Not all lead SNPs matched the VCF')

print('EXACT_DS_VARIANTS_MATCH =', len(dosages), flush=True)

# Full-sample reproducibility and seven leave-one-source-out fits.
calibration = []
loso_rows = []
gene_summaries = []

for model in ('M1', 'M1_PEER5'):

    summary, permutation = models[model]
    C = cov[model]

    for row in summary:
        gene = row['gene_id']
        variant = row['permutation_lead_variant']

        y = phenotypes[gene]
        x = dosages[variant]

        full = fit(y, x, C)

        original_p = float(permutation[gene]['pval_nominal'])
        original_slope = float(permutation[gene]['slope'])

        if full['status'] == 'OK':
            p_diff = abs(
                math.log10(max(full['p'],1e-300)) -
                math.log10(max(original_p,1e-300))
            )
            slope_diff = (
                abs(full['slope']-original_slope) /
                max(0.01,abs(original_slope))
            )
        else:
            p_diff = math.nan
            slope_diff = math.nan

        calibration.append({
            'model': model,
            'gene_id': gene,
            'lead_variant': variant,
            'FastQTL_slope': original_slope,
            'OLS_slope': full['slope'],
            'FastQTL_p': original_p,
            'OLS_p': full['p'],
            'abs_log10_p_diff': p_diff,
            'relative_slope_diff': slope_diff,
            'n_full': full['n'],
            'full_status': full['status']
        })

        current = []

        for group in sorted(group_sizes):
            keep = groups != group

            result = fit(
                y[keep],
                x[keep],
                C[keep, :]
            )

            item = {
                'model': model,
                'gene_id': gene,
                'lead_variant': variant,
                'removed_group': group,
                'removed_n': group_sizes[group],
                'n': result['n'],
                'df': result['df'],
                'cov_rank': result['rank'],
                'slope': result['slope'],
                'slope_se': result['se'],
                'p': result['p'],
                'ci95_low': result['ci95_low'],
                'ci95_high': result['ci95_high'],
                'status': result['status'],
                'full_slope': full['slope'],
                'full_p': full['p']
            }

            loso_rows.append(item)
            current.append(item)

        usable = [
            r for r in current
            if r['status'] == 'OK'
        ]

        if full['status'] == 'OK' and full['slope'] != 0:
            ratios = [
                abs(r['slope']/full['slope'])
                for r in usable
            ]
            same_direction = sum(
                r['slope'] * full['slope'] > 0
                for r in usable
            )
        else:
            ratios = []
            same_direction = 0

        if ratios:
            stable = (
                len(usable)==7 and
                same_direction==7 and
                min(ratios)>=0.5 and
                max(ratios)<=2.0
            )
            label = (
                'DIRECTION_AND_MAGNITUDE_STABLE'
                if stable else 'REVIEW'
            )
        else:
            label = 'UNESTIMABLE'

        gene_summaries.append({
            'model': model,
            'gene_id': gene,
            'lead_variant': variant,
            'qval': row['qval'],
            'full_slope': full['slope'],
            'full_p': full['p'],
            'full_n': full['n'],
            'full_status': full['status'],
            'n_loso_estimable': len(usable),
            'n_loso_same_direction': same_direction,
            'n_loso_p_lt_0_05': sum(
                r['p']<0.05 for r in usable
            ),
            'min_abs_slope_ratio':
                min(ratios) if ratios else math.nan,
            'max_abs_slope_ratio':
                max(ratios) if ratios else math.nan,
            'sensitivity_label': label
        })

    print(
        f'MODEL={model} COMPLETED_GENES={len(summary)}',
        flush=True
    )

write_table(
    'full_model_FastQTL_reproducibility.tsv',
    calibration
)

write_table(
    'all_447_leads_LOSO_7groups.tsv',
    loso_rows
)

write_table(
    'all_447_leads_LOSO_summary.tsv',
    gene_summaries
)

# Do not accept LOSO results unless full-data OLS
# reproduces the original FastQTL associations adequately.
valid = [
    r for r in calibration
    if r['full_status']=='OK' and r['n_full']==58
]

passed = [
    r for r in valid
    if r['abs_log10_p_diff']<=0.1
    and r['relative_slope_diff']<=0.05
]

print('REPRODUCIBILITY_FULL_DS_N =', len(valid))
print('AGREEMENT_WITHIN_TOLERANCE =', len(passed))

if len(valid)<400 or len(passed)/max(1,len(valid))<0.95:
    raise RuntimeError(
        'FASTQTL_REPRODUCIBILITY_QC_FAILED: '
        'inspect full_model_FastQTL_reproducibility.tsv'
    )

for model in ('M1','M1_PEER5'):
    rows = [
        r for r in gene_summaries
        if r['model']==model
    ]

    same7 = sum(
        r['n_loso_estimable']==7 and
        r['n_loso_same_direction']==7
        for r in rows
    )

    stable = sum(
        r['sensitivity_label']=='DIRECTION_AND_MAGNITUDE_STABLE'
        for r in rows
    )

    print(
        f'{model} TOTAL={len(rows)} '
        f'LOSO_SAME_DIRECTION_ALL7={same7} '
        f'HEURISTIC_STABLE={stable}',
        flush=True
    )

print('STEP25B_LOSO_SUCCESS', flush=True)

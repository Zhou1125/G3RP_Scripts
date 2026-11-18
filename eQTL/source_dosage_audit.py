#!/usr/bin/env python3
import csv
import math
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np

P = Path('.')

TABIX = Path(
    './tools/bin/tabix'
)

VCF = Path(
    'inputs/genotypes/skin58.vcf.gz'
)

COV = (
    P / '04_covariates/17_GTEx/'
    'M1_Source_GenoPC3.covariates.tsv'
)

OUT = P / '05_fastqtl/25_source_dosage_audit'
OUT.mkdir(parents=True, exist_ok=True)

for f in [TABIX, VCF, Path(str(VCF) + '.tbi'), COV]:
    if not f.is_file() or f.stat().st_size == 0:
        raise RuntimeError(f'INPUT_MISSING: {f}')

print('===== STEP25A SOURCE DOSAGE AUDIT =====', flush=True)

# Identify the six source-group dummy covariates.
with COV.open(newline='') as fh:
    reader = csv.reader(fh, delimiter='\t')
    header = next(reader)
    samples = header[1:]

    if len(samples) != 58 or len(set(samples)) != 58:
        raise RuntimeError('Expected 58 unique covariate samples')

    indicators = []
    all_names = []

    for row in reader:
        if not row:
            continue

        if len(row) != len(header):
            raise RuntimeError('Malformed covariate row')

        all_names.append(row[0])
        x = np.asarray([float(v) for v in row[1:]])

        if not np.isfinite(x).all():
            raise RuntimeError('Nonfinite covariate values')

        binary = np.all(
            (np.abs(x) < 1e-8) |
            (np.abs(x - 1) < 1e-8)
        )

        if binary and 0 < x.sum() < len(x):
            indicators.append((row[0], x))

if len(indicators) != 6:
    raise RuntimeError(
        'Expected six source-group indicator rows; '
        f'found {len(indicators)}: {all_names}'
    )

B = np.vstack([x for _, x in indicators])

if np.any(B.sum(axis=0) > 1 + 1e-8):
    raise RuntimeError('Source indicators are not one-hot')

groups = []

for col in B.T:
    idx = np.flatnonzero(col > 0.5)
    group = (
        indicators[int(idx[0])][0]
        if len(idx) else 'Source_baseline'
    )
    groups.append(group)

group_counts = Counter(groups)

expected_sizes = sorted([11, 7, 16, 6, 6, 6, 6])

if (
    len(group_counts) != 7 or
    sorted(group_counts.values()) != expected_sizes
):
    raise RuntimeError(
        f'Source group counts mismatch: {dict(group_counts)}'
    )

print('SOURCE_GROUPS =', dict(sorted(group_counts.items())))

# Read exact sample order from VCF.
proc = subprocess.run(
    [str(TABIX), '-H', str(VCF)],
    capture_output=True,
    text=True
)

if proc.returncode:
    raise RuntimeError(proc.stderr[:1000])

chrom_headers = [
    line for line in proc.stdout.splitlines()
    if line.startswith('#CHROM\t')
]

if len(chrom_headers) != 1:
    raise RuntimeError('Invalid VCF header')

vcf_samples = chrom_headers[0].split('\t')[9:]

if (
    len(vcf_samples) != 58 or
    set(vcf_samples) != set(samples)
):
    raise RuntimeError('VCF/covariate sample mismatch')

indices = [vcf_samples.index(s) for s in samples]

print('VCF_SAMPLE_MATCH = 58')

# Load the 447 model-gene lead records.
records = []

for model, expected_n in [('M1', 175), ('M1_PEER5', 272)]:

    f = (
        P / '05_fastqtl/24_significant' /
        f'{model}.egene_pair_summary.tsv'
    )

    with f.open(newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')

        required = {
            'gene_id', 'qval',
            'permutation_lead_variant',
            'num_significant_pairs'
        }

        if not required.issubset(reader.fieldnames or []):
            raise RuntimeError(f'Missing columns: {f}')

        rows = list(reader)

    if (
        len(rows) != expected_n or
        len({r['gene_id'] for r in rows}) != expected_n
    ):
        raise RuntimeError(f'{model} eGene count mismatch')

    for row in rows:
        if (
            float(row['qval']) > 0.05 or
            int(row['num_significant_pairs']) < 1
        ):
            raise RuntimeError('Invalid significant eGene')

        records.append((model, row))

unique_leads = sorted({
    row['permutation_lead_variant']
    for _, row in records
})

print('MODEL_GENE_LEAD_ROWS =', len(records))
print('UNIQUE_LEAD_VARIANTS =', len(unique_leads))

# Create a two-column 1-based region file.
targets = {}

for variant in unique_leads:

    tokens = variant.split('_', 3)

    if len(tokens) != 4 or not tokens[1].isdigit():
        raise RuntimeError(f'Invalid variant ID: {variant}')

    chrom, pos, ref, alt = tokens
    targets[(chrom, pos, ref, alt)] = variant

positions = sorted(
    {(chrom, int(pos)) for chrom, pos, _, _ in targets},
    key=lambda x: (
        int(x[0]) if x[0].isdigit() else 999,
        x[0], x[1]
    )
)

regions = OUT / 'unique_lead_positions_1based.tsv'

with regions.open('w') as fh:
    for chrom, pos in positions:
        fh.write(f'{chrom}\t{pos}\n')

# Query every unique position in one indexed VCF operation.
proc = subprocess.run(
    [str(TABIX), '-R', str(regions), str(VCF)],
    capture_output=True,
    text=True
)

if proc.returncode:
    raise RuntimeError('tabix -R: ' + proc.stderr[:1000])

matched = {}

for line in proc.stdout.splitlines():

    fields = line.split('\t')

    if len(fields) < 9 + len(vcf_samples):
        raise RuntimeError('Malformed VCF row')

    key = (fields[0], fields[1], fields[3], fields[4])

    if key in targets:
        if key in matched:
            raise RuntimeError('Duplicate VCF allele')
        matched[key] = fields

if len(matched) != len(targets):

    missing = [
        v for k, v in targets.items()
        if k not in matched
    ]

    raise RuntimeError(
        f'Exact VCF allele mismatches: {missing[:10]}'
    )

print('VCF_EXACT_LEAD_ALLELES =', len(matched))

cache = {}

def get_dosage(variant):

    if variant in cache:
        return cache[variant]

    key = tuple(variant.split('_', 3))
    fields = matched[key]

    fmt = fields[8].split(':')

    if 'DS' not in fmt:
        raise RuntimeError(f'DS missing: {variant}')

    ds_index = fmt.index('DS')
    gp_index = fmt.index('GP') if 'GP' in fmt else None

    ds = []
    gp = []

    for i in indices:

        parts = fields[9 + i].split(':')

        raw = (
            parts[ds_index]
            if ds_index < len(parts) else '.'
        )

        if raw in ('.', '', 'NA'):
            ds.append(np.nan)
        else:
            if ',' in raw:
                raise RuntimeError('Ambiguous multi-allelic DS')

            value = float(raw)

            if not -1e-4 <= value <= 2 + 1e-4:
                raise RuntimeError('Dosage outside [0,2]')

            ds.append(value)

        if (
            gp_index is None or
            gp_index >= len(parts) or
            parts[gp_index] in ('.', '')
        ):
            gp.append(np.nan)
        else:
            values = parts[gp_index].split(',')

            gp.append(
                max(map(float, values))
                if len(values) == 3 and '.' not in values
                else np.nan
            )

    result = (
        np.asarray(ds, dtype=float),
        np.asarray(gp, dtype=float)
    )

    cache[variant] = result
    return result

summaries = []
group_details = []

for index, (model, row) in enumerate(records, 1):

    variant = row['permutation_lead_variant']
    dosage, gp = get_dosage(variant)

    valid = np.isfinite(dosage)
    x = dosage[valid]

    if len(x) < 2:
        raise RuntimeError(f'Insufficient dosage: {variant}')

    total_ss = float(np.square(x - x.mean()).sum())
    within_ss = 0.0
    informative = 0

    for group in sorted(group_counts):

        mask = np.asarray(
            [g == group for g in groups]
        ) & valid

        values = dosage[mask]
        n = len(values)

        mac = (
            min(float(values.sum()), float(2*n-values.sum()))
            if n else 0.0
        )

        variance = (
            float(values.var()) if n else math.nan
        )

        # Exploratory thresholds, not formal association tests.
        variable = bool(
            n >= 5 and mac >= 2 and variance >= 0.02
        )

        informative += int(variable)

        if n:
            within_ss += float(
                np.square(values-values.mean()).sum()
            )

        group_details.append({
            'model': model,
            'gene_id': row['gene_id'],
            'lead_variant': variant,
            'source_group': group,
            'n_group': group_counts[group],
            'n_valid_DS': n,
            'alt_dosage_frequency':
                values.mean()/2 if n else math.nan,
            'dosage_minor_allele_count': mac,
            'dosage_variance': variance,
            'n_DS_above_0_5': int((values > 0.5).sum()),
            'within_group_informative': int(variable)
        })

    r2 = (
        1 - within_ss / total_ss
        if total_ss > 1e-12 else math.nan
    )

    gp_valid = np.isfinite(gp)

    gp_fraction = (
        float((gp[gp_valid] >= 0.9).mean())
        if gp_valid.any() else math.nan
    )

    r2 = (
        min(1.0, max(0.0, r2))
        if math.isfinite(r2) else math.nan
    )

    review = (
        informative < 2 or
        (math.isfinite(r2) and r2 >= 0.80)
    )

    summaries.append({
        'model': model,
        'gene_id': row['gene_id'],
        'qval': row['qval'],
        'lead_variant': variant,
        'num_significant_pairs': row['num_significant_pairs'],
        'n_valid_DS': int(valid.sum()),
        'n_missing_DS': int((~valid).sum()),
        'global_DS_MAF': min(
            float(x.mean()/2),
            float(1-x.mean()/2)
        ),
        'R2_dosage_explained_by_source': r2,
        'within_source_variance_fraction':
            within_ss / total_ss
            if total_ss > 1e-12 else math.nan,
        'n_informative_source_groups': informative,
        'GP_observed_n': int(gp_valid.sum()),
        'fraction_GP_max_ge_0_9': gp_fraction,
        'review_flag':
            'REVIEW' if review else 'MULTIGROUP_VARIATION'
    })

    if index % 50 == 0:
        print(
            f'PROGRESS = {index}/{len(records)}',
            flush=True
        )

def write_table(path, rows):

    tmp = Path(str(path) + '.partial')

    with tmp.open('w', newline='') as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=list(rows[0]),
            delimiter='\t',
            lineterminator='\n'
        )

        writer.writeheader()
        writer.writerows(rows)

    tmp.replace(path)

write_table(
    OUT / 'lead_source_dosage_summary.tsv',
    summaries
)

write_table(
    OUT / 'lead_by_source_group.tsv',
    group_details
)

for model in ('M1', 'M1_PEER5'):

    rows = [
        r for r in summaries
        if r['model'] == model
    ]

    n_variable = sum(
        r['n_informative_source_groups'] >= 2
        for r in rows
    )

    n_high_r2 = sum(
        math.isfinite(r['R2_dosage_explained_by_source'])
        and r['R2_dosage_explained_by_source'] >= 0.8
        for r in rows
    )

    n_review = sum(
        r['review_flag'] == 'REVIEW'
        for r in rows
    )

    print(
        f'{model} LEADS={len(rows)} '
        f'GROUPS_VAR_GE2={n_variable} '
        f'SOURCE_R2_GE0.8={n_high_r2} '
        f'REVIEW={n_review}',
        flush=True
    )

print('STEP25A_SOURCE_DOSAGE_AUDIT_SUCCESS', flush=True)

#!/usr/bin/env python3
import csv
import gzip
import math
from collections import Counter
from pathlib import Path

P = Path('.')
FDR = 0.05
NOMINAL_CAP = 1e-4
MODELS = ('M1', 'M1_PEER5')

OUT = P / '05_fastqtl/24_significant'
OUT.mkdir(parents=True, exist_ok=True)

NOMINAL_COLS = [
    'gene_id', 'variant_id', 'tss_distance',
    'ma_samples', 'ma_count', 'maf',
    'pval_nominal', 'slope', 'slope_se'
]

EXTRA_COLS = [
    'egene_qval',
    'egene_pval_beta',
    'gene_nominal_threshold',
    'permutation_lead_variant',
    'permutation_lead_nominal_p'
]


def load_egenes(model):

    f = (
        P / '05_fastqtl/22_fdr'
        / model / 'egenes_FDR005.tsv.gz'
    )

    if not f.is_file():
        raise FileNotFoundError(f)

    genes = {}

    with gzip.open(f, 'rt', newline='') as handle:

        reader = csv.DictReader(handle, delimiter='\t')

        mandatory = {
            'gene_id', 'qval', 'pval_beta',
            'variant_id', 'pval_nominal',
            'pval_nominal_threshold'
        }

        if (
            reader.fieldnames is None
            or not mandatory.issubset(reader.fieldnames)
        ):
            raise RuntimeError(
                f'{model}: missing eGene columns'
            )

        for r in reader:

            gene = r['gene_id']

            if gene in genes:
                raise RuntimeError(
                    f'{model}: duplicate eGene {gene}'
                )

            qval = float(r['qval'])
            threshold = float(
                r['pval_nominal_threshold']
            )

            if not (
                math.isfinite(qval)
                and 0 <= qval <= FDR
            ):
                raise RuntimeError(
                    f'{model}: invalid eGene q-value {gene}'
                )

            if not (
                math.isfinite(threshold)
                and 0 < threshold <= NOMINAL_CAP
            ):
                raise RuntimeError(
                    f'{model}: invalid threshold '
                    f'{gene} {threshold}'
                )

            genes[gene] = r

    if not genes:
        raise RuntimeError(f'{model}: zero eGenes')

    maximum = max(
        float(x['pval_nominal_threshold'])
        for x in genes.values()
    )

    print(
        f'{model} EGENES={len(genes)} '
        f'MAX_GENE_CUTOFF={maximum}',
        flush=True
    )

    return genes


def process_model(model, genes):

    f = (
        P / '05_fastqtl/21_full'
        / model / 'nominal'
        / (model + '.nominal.tsv.gz')
    )

    if not f.is_file():
        raise FileNotFoundError(f)

    out = OUT / (
        model + '.significant_pairs.tsv.gz'
    )

    counts = Counter()
    pairs = {}
    best = {}
    seen_nominal = set()

    total = 0
    selected_genes = 0

    with (
        gzip.open(f, 'rt', newline='') as source,
        gzip.open(out, 'wt', newline='') as dest
    ):

        reader = csv.DictReader(
            source, delimiter='\t'
        )

        if reader.fieldnames != NOMINAL_COLS:
            raise RuntimeError(
                f'{model}: unexpected nominal columns: '
                f'{reader.fieldnames}'
            )

        writer = csv.DictWriter(
            dest,
            fieldnames=NOMINAL_COLS + EXTRA_COLS,
            delimiter='\t',
            lineterminator='\n'
        )

        writer.writeheader()

        for r in reader:

            total += 1

            gene = r['gene_id']
            variant = r['variant_id']

            pair = (gene, variant)

            if pair in seen_nominal:
                raise RuntimeError(
                    f'{model}: duplicate nominal pair {pair}'
                )

            seen_nominal.add(pair)

            pv = float(r['pval_nominal'])

            if not (
                math.isfinite(pv)
                and 0 <= pv <= NOMINAL_CAP + 1e-10
            ):
                raise RuntimeError(
                    f'{model}: invalid nominal p-value '
                    f'{pair} {pv}'
                )

            if gene not in genes:
                continue

            selected_genes += 1

            g = genes[gene]

            threshold = float(
                g['pval_nominal_threshold']
            )

            # GTEx-style gene-specific significance
            if pv >= threshold:
                continue

            slope = float(r['slope'])

            if not math.isfinite(slope):
                raise RuntimeError(
                    f'{model}: invalid slope {pair}'
                )

            counts[gene] += 1

            pairs[pair] = (pv, slope)

            if (
                gene not in best
                or pv < best[gene][0]
            ):
                best[gene] = (pv, variant)

            r.update({
                'egene_qval':
                    g['qval'],
                'egene_pval_beta':
                    g['pval_beta'],
                'gene_nominal_threshold':
                    g['pval_nominal_threshold'],
                'permutation_lead_variant':
                    g['variant_id'],
                'permutation_lead_nominal_p':
                    g['pval_nominal']
            })

            writer.writerow(r)

    missing = set(genes) - set(counts)

    missing_lead = [
        gene
        for gene, g in genes.items()
        if (gene, g['variant_id']) not in pairs
    ]

    # Gene-level summary, retaining all significant eGenes
    summary = OUT / (
        model + '.egene_pair_summary.tsv'
    )

    fields = [
        'gene_id',
        'qval',
        'pval_beta',
        'pval_nominal_threshold',
        'permutation_lead_variant',
        'permutation_lead_nominal_p',
        'num_significant_pairs',
        'min_pair_nominal_p',
        'min_pair_variant_id',
        'permutation_lead_present',
        'min_variant_matches_perm_lead'
    ]

    with summary.open('w', newline='') as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter='\t',
            lineterminator='\n'
        )

        writer.writeheader()

        for gene in sorted(genes):

            g = genes[gene]
            b = best.get(gene)

            writer.writerow({
                'gene_id': gene,
                'qval': g['qval'],
                'pval_beta': g['pval_beta'],
                'pval_nominal_threshold':
                    g['pval_nominal_threshold'],
                'permutation_lead_variant':
                    g['variant_id'],
                'permutation_lead_nominal_p':
                    g['pval_nominal'],
                'num_significant_pairs':
                    counts[gene],
                'min_pair_nominal_p':
                    b[0] if b else 'NA',
                'min_pair_variant_id':
                    b[1] if b else 'NA',
                'permutation_lead_present':
                    (gene, g['variant_id']) in pairs,
                'min_variant_matches_perm_lead':
                    (
                        b is not None
                        and b[1] == g['variant_id']
                    )
            })

    print(
        f'{model} NOMINAL_ROWS={total} '
        f'NOMINAL_IN_EGENES={selected_genes} '
        f'SIGNIFICANT_PAIRS={len(pairs)} '
        f'EGENES_WITH_PAIRS={len(counts)} '
        f'LEAD_MISSING={len(missing_lead)}',
        flush=True
    )

    if missing or missing_lead:
        raise RuntimeError(
            f'{model}: missing eGenes='
            f'{sorted(missing)[:8]}, '
            f'permutation leads={missing_lead[:8]}. '
            'Check inputs.'
        )

    return pairs


print(
    '===== STEP24 FASTQTL SIGNIFICANT PAIRS =====',
    flush=True
)

egenes = {
    name: load_egenes(name)
    for name in MODELS
}

pairsets = {
    name: process_model(name, egenes[name])
    for name in MODELS
}

A = pairsets['M1']
B = pairsets['M1_PEER5']

shared_pairs = set(A) & set(B)

shared_genes = (
    set(egenes['M1'])
    & set(egenes['M1_PEER5'])
)

print('SHARED_EGENES =', len(shared_genes))
print('SHARED_EXACT_PAIRS =', len(shared_pairs))
print('M1_ONLY_PAIRS =', len(A) - len(shared_pairs))
print('PEER5_ONLY_PAIRS =', len(B) - len(shared_pairs))

# Compare identical gene-variant pairs across models
shared_out = OUT / (
    'shared_exact_significant_pairs.tsv.gz'
)

shared_fields = [
    'gene_id',
    'variant_id',
    'M1_pval_nominal',
    'M1_slope',
    'M1_PEER5_pval_nominal',
    'M1_PEER5_slope',
    'same_effect_direction'
]

with gzip.open(
    shared_out, 'wt', newline=''
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=shared_fields,
        delimiter='\t',
        lineterminator='\n'
    )

    writer.writeheader()

    for gene, variant in sorted(shared_pairs):

        p1, s1 = A[(gene, variant)]
        p2, s2 = B[(gene, variant)]

        writer.writerow({
            'gene_id': gene,
            'variant_id': variant,
            'M1_pval_nominal': p1,
            'M1_slope': s1,
            'M1_PEER5_pval_nominal': p2,
            'M1_PEER5_slope': s2,
            'same_effect_direction': s1 * s2 > 0
        })


same_leads = sum(
    egenes['M1'][gene]['variant_id']
    == egenes['M1_PEER5'][gene]['variant_id']
    for gene in shared_genes
)

summary_values = {
    'M1_eGenes': len(egenes['M1']),
    'M1_PEER5_eGenes': len(egenes['M1_PEER5']),
    'shared_eGenes': len(shared_genes),
    'shared_same_lead': same_leads,
    'M1_pairs': len(A),
    'M1_PEER5_pairs': len(B),
    'shared_exact_pairs': len(shared_pairs),
    'M1_only_pairs': len(A) - len(shared_pairs),
    'PEER5_only_pairs': len(B) - len(shared_pairs),
    'same_slope_sign_shared_pairs': sum(
        A[x][1] * B[x][1] > 0
        for x in shared_pairs
    )
}

summary_file = OUT / (
    'Step24_model_comparison_summary.tsv'
)

with summary_file.open(
    'w', newline=''
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=list(summary_values),
        delimiter='\t',
        lineterminator='\n'
    )

    writer.writeheader()
    writer.writerow(summary_values)

# Confirm compatibility with independently verified Step23A
if (
    len(egenes['M1']) != 175
    or len(egenes['M1_PEER5']) != 272
    or len(shared_genes) != 121
    or same_leads != 71
):
    raise RuntimeError(
        'Step24 counts disagree with Step23A'
    )

for key, value in summary_values.items():
    print(f'{key}={value}', flush=True)

print('STEP24_ALL_SUCCESS', flush=True)

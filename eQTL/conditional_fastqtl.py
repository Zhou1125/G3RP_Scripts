#!/usr/bin/env python3
import csv
import gzip
import hashlib
import json
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

P = Path('.')
OUT = P / '05_fastqtl/26_conditional'

BED = P / '03_expression/16_gtex_style/GTEx_style_TMM_INT_58.bed.gz'
BIN = P / 'tools/bin/fastQTL'

VCF = Path(
    'inputs/genotypes/skin58.vcf.gz'
)

TABIX = Path('./tools/bin/tabix')
BGZIP = Path('./tools/bin/bgzip')

COVS = {
    'M1': P / '04_covariates/17_GTEx/M1_Source_GenoPC3.covariates.tsv',
    'M1_PEER5': P / '04_covariates/18_PEER/PEER5/M1_PEER5.covariates.tsv',
}

NGENES = {'M1': 175, 'M1_PEER5': 272}
CHUNKS = 32
WORKERS = 4

FIELDS = (
    'gene_id', 'rank', 'variant_id', 'tss_distance',
    'pval_nominal', 'pval_conditional', 'selected_as_lead'
)

def sha256(data):
    return hashlib.sha256(data).hexdigest()

def require(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f'Missing input: {path}')

def check_chunk(file, allowed):
    if not file.is_file() or file.stat().st_size == 0:
        return False

    try:
        with gzip.open(file, 'rt') as f:
            for line in f:
                if not line.strip():
                    continue

                v = line.rstrip('\n').split()

                if len(v) != 7 or v[0] not in allowed:
                    return False

                if int(v[1]) < 0 or v[6] not in ('0', '1'):
                    return False

                for s in (v[4], v[5]):
                    if not math.isfinite(float(s)):
                        return False

        return True

    except (OSError, EOFError, UnicodeError, ValueError, IndexError):
        return False

def prepare(model):
    root = OUT / model
    root.mkdir(parents=True, exist_ok=True)

    egfile = P / f'05_fastqtl/22_fdr/{model}/egenes_FDR005.tsv.gz'
    require(egfile)

    genes = {}

    with gzip.open(egfile, 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            gene = r['gene_id']
            threshold = float(r['pval_nominal_threshold'])

            if gene in genes or not (0.0 < threshold < 0.0001):
                raise RuntimeError(
                    f'Invalid gene threshold: {model} {gene} {threshold}'
                )

            genes[gene] = (threshold, r['qval'], r['variant_id'])

    if len(genes) != NGENES[model]:
        raise RuntimeError(f'{model}: eGene count mismatch')

    selected = []

    with gzip.open(BED, 'rt') as fh:
        header = next(fh)

        if not header.startswith('#') or len(header.rstrip('\n').split('\t')) != 62:
            raise RuntimeError('Unexpected BED header')

        selected.append(header)
        used = set()

        for line in fh:
            v = line.split('\t', 5)

            if len(v) < 5:
                raise RuntimeError('Malformed BED line')

            gene = v[3]

            if gene in genes:
                if gene in used:
                    raise RuntimeError(f'Duplicate BED gene: {gene}')

                selected.append(line)
                used.add(gene)

    if used != set(genes):
        raise RuntimeError(f'{model}: selected eGenes missing from BED')

    bed_bytes = ''.join(selected).encode()

    thresh_bytes = ''.join(
        f'{g}\t{genes[g][0]:.17g}\n'
        for g in sorted(genes)
    ).encode()

    signature = {
        'model': model,
        'bed_sha256': sha256(bed_bytes),
        'threshold_sha256': sha256(thresh_bytes),
        'vcf': str(VCF),
        'cov': str(COVS[model]),
        'bin': str(BIN),
        'chunks': CHUNKS,
        'maf': 0.05,
        'ma_samples': 3,
        'window': 1000000,
        'map_full': True
    }

    manifest = root / 'analysis_settings.json'

    if manifest.is_file():
        if json.loads(manifest.read_text()) != signature:
            raise RuntimeError(
                f'{model}: prior analysis settings differ'
            )
    else:
        if list((root / 'chunks').glob('chunk_*.ok')):
            raise RuntimeError('Existing chunks lack parameter manifest')

        manifest.write_text(
            json.dumps(signature, indent=2) + '\n'
        )

    bed_out = root / 'egenes_only.bed.gz'
    threshold_out = root / 'egenes_thresholds.txt'

    if (
        not bed_out.is_file()
        or not Path(str(bed_out) + '.tbi').is_file()
    ):
        with bed_out.open('wb') as dest:
            cp = subprocess.run(
                [str(BGZIP), '-c'],
                input=bed_bytes,
                stdout=dest,
                stderr=subprocess.PIPE
            )

        if cp.returncode != 0:
            raise RuntimeError(
                'BGZIP_FAILED: ' + cp.stderr.decode(errors='replace')
            )

        subprocess.run(
            [str(TABIX), '-f', '-p', 'bed', str(bed_out)],
            check=True
        )

    else:
        with gzip.open(bed_out, 'rb') as f:
            if sha256(f.read()) != signature['bed_sha256']:
                raise RuntimeError('Stored BED differs from current input')

    threshold_out.write_bytes(thresh_bytes)

    print(f'{model} PREPARED_GENES={len(genes)}', flush=True)

    return root, genes, bed_out, threshold_out

def run(model):
    root, genes, bed, threshold = prepare(model)

    directory = root / 'chunks'
    directory.mkdir(exist_ok=True)

    def task(i):
        prefix = directory / f'chunk_{i:03d}'

        result = Path(str(prefix) + '.txt.gz')
        ok = Path(str(prefix) + '.ok')

        if ok.is_file() and check_chunk(result, genes):
            return i, 'REUSED'

        ok.unlink(missing_ok=True)

        cmd = [
            str(BIN),
            '--vcf', str(VCF),
            '--bed', str(bed),
            '--cov', str(COVS[model]),
            '--map', str(threshold),
            '--map-full',
            '--window', '1000000',
            '--maf-threshold', '0.05',
            '--ma-sample-threshold', '3',
            '--chunk', str(i), str(CHUNKS),
            '--seed', str(20261101 + i),
            '--log', str(prefix) + '.fastqtl.log',
            '--out', str(result)
        ]

        with Path(str(prefix) + '.runner.log').open('w') as f:
            cp = subprocess.run(
                cmd, stdout=f, stderr=subprocess.STDOUT
            )

        if cp.returncode != 0 or not check_chunk(result, genes):
            raise RuntimeError(
                f'{model} chunk {i} failed; see {prefix}.runner.log'
            )

        ok.write_text('OK\n')
        return i, 'FINISHED'

    failures = []

    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futures = [
            pool.submit(task, i)
            for i in range(1, CHUNKS + 1)
        ]

        for future in as_completed(futures):
            try:
                i, status = future.result()
                print(
                    f'{model} CHUNK={i:03d} {status}',
                    flush=True
                )
            except Exception as e:
                print(f'{model} ERROR={e}', flush=True)
                failures.append(str(e))

    if failures:
        raise RuntimeError(
            f'{model}: {len(failures)} mapping chunks failed'
        )

    allrows = []
    candidate_leads = []

    for i in range(1, CHUNKS + 1):
        f = directory / f'chunk_{i:03d}.txt.gz'

        with gzip.open(f, 'rt') as fh:
            for line in fh:
                if not line.strip():
                    continue

                v = line.split()
                r = dict(zip(FIELDS, v))

                gene = r['gene_id']
                r['model'] = model
                r['gene_qval'] = genes[gene][1]
                r['gene_nominal_threshold'] = genes[gene][0]
                r['original_permutation_lead'] = genes[gene][2]

                allrows.append(r)

                if r['selected_as_lead'] == '1':
                    candidate_leads.append(r)

    seen = set()
    leads_by_gene = {}

    for r in candidate_leads:
        key = (r['gene_id'], int(r['rank']))

        if key in seen:
            raise RuntimeError(f'Duplicate gene-rank: {key}')

        seen.add(key)
        leads_by_gene.setdefault(r['gene_id'], []).append(r)

    by_gene = []

    for gene, (threshold, qval, original_lead) in sorted(genes.items()):
        rows = leads_by_gene.get(gene, [])

        sig = sum(
            float(r['pval_conditional']) <= threshold
            for r in rows
        )

        first = next(
            (r for r in rows if r['rank'] == '0'),
            None
        )

        by_gene.append({
            'gene_id': gene,
            'gene_qval': qval,
            'original_permutation_lead': original_lead,
            'n_forward_selected_signals': len(rows),
            'n_backward_threshold_pass': sig,
            'rank0_matches_permutation_lead':
                'NA' if first is None
                else int(first['variant_id'] == original_lead)
        })

    fields = list(FIELDS) + [
        'model', 'gene_qval', 'gene_nominal_threshold',
        'original_permutation_lead'
    ]

    def write_table(path, rows, columns):
        tmp = Path(str(path) + '.partial')

        with tmp.open('w', newline='') as fh:
            w = csv.DictWriter(
                fh, fieldnames=columns,
                delimiter='\t', lineterminator='\n'
            )
            w.writeheader()
            w.writerows(rows)

        os.replace(tmp, path)

    write_table(
        root / 'conditional_candidate_leads.tsv',
        candidate_leads, fields
    )

    write_table(
        root / 'per_egene_conditional_summary.tsv',
        by_gene, list(by_gene[0])
    )

    path = root / 'conditional_all_candidates.tsv.gz'
    tmp = Path(str(path) + '.partial')

    with gzip.open(tmp, 'wt') as fh:
        w = csv.DictWriter(
            fh, fieldnames=fields,
            delimiter='\t', lineterminator='\n'
        )
        w.writeheader()
        w.writerows(allrows)

    os.replace(tmp, path)

    backward_pass = sum(
        r['n_backward_threshold_pass']
        for r in by_gene
    )

    multi = sum(
        r['n_forward_selected_signals'] >= 2
        for r in by_gene
    )

    print(
        f'{model} EGENES={len(genes)} '
        f'GENES_WITH_SIGNALS={len(leads_by_gene)} '
        f'FORWARD_SIGNALS={len(candidate_leads)} '
        f'BACKWARD_THRESHOLD_PASS={backward_pass} '
        f'MULTISIGNAL_GENES={multi} '
        f'NO_SIGNAL_GENES={len(genes)-len(leads_by_gene)}',
        flush=True
    )

    return {
        'model': model,
        'egenes': len(genes),
        'genes_with_candidate_signals': len(leads_by_gene),
        'forward_candidate_signals': len(candidate_leads),
        'backward_threshold_pass': backward_pass,
        'genes_with_2plus_signals': multi
    }

if __name__ == '__main__':
    for file in [
        BED, BIN, VCF, Path(str(VCF) + '.tbi'),
        TABIX, BGZIP, *COVS.values()
    ]:
        require(file)

    OUT.mkdir(parents=True, exist_ok=True)

    stats = [
        run(model)
        for model in ('M1', 'M1_PEER5')
    ]

    dest = OUT / 'Step26_conditional_model_summary.tsv'

    with dest.open('w', newline='') as fh:
        w = csv.DictWriter(
            fh, fieldnames=list(stats[0]),
            delimiter='\t', lineterminator='\n'
        )
        w.writeheader()
        w.writerows(stats)

    print('STEP26_CONDITIONAL_SUCCESS', flush=True)

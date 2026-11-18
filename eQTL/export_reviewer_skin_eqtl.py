#!/usr/bin/env python3
import csv
import gzip
import math
import subprocess
from pathlib import Path

P = Path('.')
OUT = P / '05_fastqtl/30_reviewer_skin_eqtl'
OUT.mkdir(parents=True, exist_ok=True)

BED = P / '03_expression/16_gtex_style/GTEx_style_TMM_INT_58.bed.gz'
COV = P / '04_covariates/17_GTEx/M1_Source_GenoPC3.covariates.tsv'
VCF = Path(
    'inputs/genotypes/skin58.vcf.gz'
)
TABIX = './tools/bin/tabix'

with gzip.open(BED, 'rt') as fh:
    header = next(fh).rstrip('\n').split('\t')
    samples = header[4:]
    if len(samples) != 58 or len(set(samples)) != 58:
        raise RuntimeError('BED sample count error')
    asip_rows = [
        v for line in fh
        if (v := line.rstrip('\n').split('\t'))[3] == 'ASIP'
    ]

if len(asip_rows) != 1:
    raise RuntimeError('ASIP must occur once in BED')

expression = dict(zip(samples, map(float, asip_rows[0][4:])))

with COV.open(newline='') as fh:
    reader = csv.reader(fh, delimiter='\t')
    cov_samples = next(reader)[1:]
    cov_rows = list(reader)

if set(cov_samples) != set(samples):
    raise RuntimeError('Covariate sample mismatch')

indicators = []
for row in cov_rows:
    x = [float(z) for z in row[1:]]
    if len(x) != 58:
        raise RuntimeError('Invalid covariate row')
    if all(abs(v)<1e-8 or abs(v-1)<1e-8 for v in x) \
       and 0<sum(x)<58:
        indicators.append((row[0], dict(zip(cov_samples, x))))

if len(indicators) != 6:
    raise RuntimeError('Expected six source-group indicators')

groups = {}
for s in samples:
    chosen = [name for name, vals in indicators if vals[s]>0.5]
    if len(chosen)>1:
        raise RuntimeError('Overlapping source groups')
    groups[s] = chosen[0] if chosen else 'Source_baseline'

vcf_header = subprocess.run(
    [TABIX, '-H', str(VCF)],
    capture_output=True, text=True, check=True
).stdout

h = [
    x for x in vcf_header.splitlines()
    if x.startswith('#CHROM\t')
]
if len(h)!=1:
    raise RuntimeError('VCF header error')

vcf_samples = h[0].split('\t')[9:]
if len(vcf_samples)!=58 or set(vcf_samples)!=set(samples):
    raise RuntimeError('VCF sample mismatch')

def read_ds(chrom, pos, ref, alt):
    lines = subprocess.run(
        [TABIX, str(VCF), f'{chrom}:{pos}-{pos}'],
        capture_output=True, text=True, check=True
    ).stdout.splitlines()

    records = [x.split('\t') for x in lines]
    records = [
        r for r in records
        if (r[0],r[1],r[3],r[4])==(chrom,pos,ref,alt)
    ]
    if len(records)!=1:
        raise RuntimeError('Exact allele not found')

    r = records[0]
    fmt = r[8].split(':')
    if 'DS' not in fmt:
        raise RuntimeError('DS not found')

    j = fmt.index('DS')
    result = {}
    for s, v in zip(vcf_samples,r[9:]):
        tokens = v.split(':')
        raw = tokens[j] if j<len(tokens) else '.'
        if ',' in raw:
            raise RuntimeError('Multiallelic DS')
        result[s] = (
            float(raw) if raw not in ('.','','NA')
            else math.nan
        )
    return result

lead = read_ds('13','63241845','C','T')
g128v = read_ds('13','63249342','G','T')

# Optional: derive uncorrected log2(CPM+1) from raw featureCounts.
logcpm = {}
counts_file = P/'03_expression/gene_counts_58.tsv.gz'

if counts_file.is_file():
    with gzip.open(counts_file,'rt') as fh:
        ch = next(fh).rstrip('\n').split('\t')
        if all(s in ch for s in samples):
            indices = [ch.index(s) for s in samples]
            totals = [0.0]*58
            asip_counts = None

            for line in fh:
                a = line.rstrip('\n').split('\t')
                if len(a)<=max(indices):
                    continue
                v = [float(a[i]) for i in indices]
                totals = [x+y for x,y in zip(totals,v)]
                if a[0]=='ASIP':
                    asip_counts = v

            if asip_counts is not None:
                for s,c,total in zip(samples,asip_counts,totals):
                    logcpm[s] = (
                        math.log2(c/total*1e6+1)
                        if total>0 else math.nan
                    )

data = OUT/'skin58_ASIP_existing_data.tsv'
with data.open('w',newline='') as fh:
    w = csv.writer(fh,delimiter='\t')
    w.writerow([
        'sample_id','source_group','ASIP_TMM_INT',
        'ASIP_log2CPM_untmm','G128V_T_DS','ASIP_lead_T_DS'
    ])
    for s in samples:
        w.writerow([
            s,groups[s],expression[s],logcpm.get(s,''),
            g128v[s],lead[s]
        ])

metadata = OUT/'skin58_metadata_TO_COMPLETE.tsv'

if not metadata.exists():
    with metadata.open('w',newline='') as fh:
        w = csv.writer(fh,delimiter='\t')
        w.writerow([
            'sample_id','source_group','coat_color',
            'breed_population','sex','age',
            'skin_anatomical_site','study_accession',
            'RNA_accession','RNA_library_batch',
            'genotype_source','genotype_validation'
        ])
        for s in samples:
            w.writerow([s,groups[s]]+['']*10)
else:
    print('Existing metadata template preserved')

# Export ASIP-specific robustness tables.
sources = [
    (
        P/'05_fastqtl/25B_LOSO/all_447_leads_LOSO_7groups.tsv',
        OUT/'ASIP_PEER5_LOSO.tsv'
    ),
    (
        P/'05_fastqtl/25_source_dosage_audit/lead_by_source_group.tsv',
        OUT/'ASIP_PEER5_source_dosage.tsv'
    )
]
for src,dst in sources:
    if not src.is_file():
        continue
    with src.open(newline='') as fi, dst.open('w',newline='') as fo:
        reader = csv.DictReader(fi,delimiter='\t')
        writer = csv.DictWriter(
            fo,fieldnames=reader.fieldnames,
            delimiter='\t',lineterminator='\n'
        )
        writer.writeheader()
        for row in reader:
            if row['model']=='M1_PEER5' and row['gene_id']=='ASIP':
                writer.writerow(row)

print('SAMPLES =',len(samples))
print('SOURCE_GROUPS =', {
    k:list(groups.values()).count(k)
    for k in sorted(set(groups.values()))
})
print('ASIP_RAW_LOGCPM_AVAILABLE =',len(logcpm))
print('DATA =',data)
print('METADATA =',metadata)
print('REVIEWER_EQTL_EXPORT_SUCCESS')

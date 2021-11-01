import argparse
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser(description='')
parser.add_argument('indir', type=Path, help='Path to the RatGTEx pipeline base directory')
parser.add_argument('outdir', type=Path, help='Output directory path')
parser.add_argument('tissues', nargs='+', type=str, help='Tissues to include')
args = parser.parse_args()

(args.outdir / 'singleTissueEqtl').mkdir(exist_ok=False)

sig = []
for tissue in args.tissues:
    print(tissue)
    d = pd.read_csv(args.indir / tissue / f'{tissue}.cis_qtl_signif.txt.gz', sep='\t')
    d = d[['phenotype_id', 'variant_id', 'pval_nominal', 'slope']]
    d = d.rename(
        columns={
            'variant_id': 'variantId',
            'pval_nominal': 'pValue',
            'slope': 'nes',
        }
    )
    d['chromosome'] = d['variantId'].str.extract(r'chr(\d+):\d+$')
    d['pos'] = d['variantId'].str.extract(r'chr\d+:(\d+)$')
    d['tissueSiteDetailId'] = tissue
    sig.append(d)

sig = pd.concat(sig)

for gene, d in sig.groupby('phenotype_id'):
    d = d.drop(columns=['phenotype_id'])
    d.to_csv(args.outdir / 'singleTissueEqtl' / f'{gene}.txt', sep='\t', index=False, float_format='%g')

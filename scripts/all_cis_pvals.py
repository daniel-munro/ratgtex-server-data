"""Save all cis p-values in zip archive of per-gene files"""
import argparse
import gzip
from pathlib import Path
import zipfile

parser = argparse.ArgumentParser(description='Save all cis p-values in zip archive of per-gene files')
parser.add_argument('indir', type=Path, help='Path to the RatGTEx pipeline base directory')
parser.add_argument("version", type=str, help="rn6 or rn7")
parser.add_argument('outdir', type=Path, help='Output directory path')
parser.add_argument('tissue', help='Tissue whose cis p-values are to be processed')
args = parser.parse_args()

infile = args.indir / args.version / args.tissue / f'{args.tissue}.cis_qtl_all_pvals.txt.gz'
outfile = args.outdir / 'cis_pvals' / f'{args.tissue}.{args.version}.zip'

seen = set()
last_gene = None
with gzip.open(infile, 'rt') as f:
    with zipfile.ZipFile(outfile, 'w', zipfile.ZIP_DEFLATED) as out:
        next(f)
        for line in f:
            (gene, var, pval) = line.split('\t')
            if gene != last_gene:
                assert gene not in seen # Entries for a gene should be consecutive
                seen.add(gene)
                if last_gene is not None:
                    out.writestr(f'{args.tissue}.{args.version}/{last_gene}.txt', pvals)
                pvals = 'variant_id\tpval_nominal\n'
                last_gene = gene
            pvals += f'{var}\t{pval}' # pval still has newline
        out.writestr(f'{args.tissue}.{args.version}/{gene}.txt', pvals)

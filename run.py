import argparse
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(description='Process RatGTEx pipeline results into data files for the site')
parser.add_argument('indir', type=Path, help='Path to the RatGTEx pipeline base directory')
parser.add_argument('outdir', type=Path, help='Output directory path')
parser.add_argument('tissues', nargs='+', type=str, help='Tissues to include')
parser.add_argument('--keep_existing', action='store_true', default=False, help='Do not regenerate output files that already exist')
args = parser.parse_args()

args.outdir.mkdir(exist_ok=True)
(args.outdir / 'eqtl').mkdir(exist_ok=True)
(args.outdir / 'cis_pvals').mkdir(exist_ok=True)

files = ['top_assoc.txt', 'eqtls_indep.txt']
files += [f'{tissue}.cis_qtl_signif.txt.gz' for tissue in args.tissues]
files += [f'{tissue}.trans_qtl_pairs.txt.gz' for tissue in args.tissues]
files = [args.outdir / 'eqtl' / x for x in files]
if not (args.keep_existing and all([file.exists() for file in files])):
    print('Assembling eQTL files')
    script = Path(__file__).parent / 'eqtl_files.R'
    subprocess.run(['Rscript', script, str(args.indir), str(args.outdir), *args.tissues])

if not (args.keep_existing and (args.outdir / 'tissueInfo.txt').exists()):
    print('Making tissue info table')
    script = Path(__file__).parent / 'tissueInfo.R'
    subprocess.run(['Rscript', script, str(args.indir), str(args.outdir), *args.tissues])

files = [args.outdir / x for x in ['gene.txt', 'autocomplete.json']]
if not (args.keep_existing and all([file.exists() for file in files])):
    print('Assembling gene info')
    script = Path(__file__).parent / 'gene.R'
    subprocess.run(['Rscript', script, str(args.indir), str(args.outdir), *args.tissues])

files = [args.outdir / x for x in ['medianGeneExpression.txt.gz', 'topExpressedGene.txt']]
if not (args.keep_existing and all([file.exists() for file in files])):
    print('Summarizing gene expression')
    script = Path(__file__).parent / 'medianGeneExpression.R'
    subprocess.run(['Rscript', script, str(args.indir), str(args.outdir), *args.tissues])

infile = args.indir / 'ref' / 'Rattus_norvegicus.Rnor_6.0.99.genes.gtf'
if not (args.keep_existing and (args.outdir / 'exon.txt').exists()):
    print('Making exon table')
    script = Path(__file__).parent / 'exon.py'
    subprocess.run(['python3', script, str(infile), str(args.outdir)])

if not (args.keep_existing and (args.outdir / 'singleTissueEqtl.zip').exists()):
    print('Assembling all significant cis associations')
    script = Path(__file__).parent / 'singleTissueEqtl.py'
    dirname = str(args.outdir / 'singleTissueEqtl')
    subprocess.run(f'mkdir -p {dirname} && rm -r {dirname}', shell=True)
    subprocess.run(['python3', script, str(args.indir), str(args.outdir), *args.tissues])
    subprocess.run(f'cd {args.outdir} && zip -rq singleTissueEqtl.zip singleTissueEqtl && rm -r singleTissueEqtl', shell=True)

files = [args.outdir / 'cis_pvals' / f'{tissue}.zip' for tissue in args.tissues]
if not (args.keep_existing and all([file.exists() for file in files])):
    print('Extracting all cis p-values')
    script = Path(__file__).parent / 'all_cis_pvals.R'
    for tissue in args.tissues:
        print(tissue)
        dirname = str(args.outdir / 'cis_pvals' / tissue)
        subprocess.run(f'mkdir -p {dirname} && rm -r {dirname} && mkdir {dirname}', shell=True)
        subprocess.run(['Rscript', script, str(args.indir), str(args.outdir), tissue])
        subprocess.run(f'cd {args.outdir / "cis_pvals"} && zip -rq {tissue}.zip {tissue} && rm -r {tissue}', shell=True)

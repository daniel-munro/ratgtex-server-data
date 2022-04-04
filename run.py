"""Call scripts to process RatGTEx pipeline results into data files for the RatGTEx site"""

import argparse
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(
    description="Process RatGTEx pipeline results into data files for the site"
)
parser.add_argument(
    "indir", type=Path, help="Path to the RatGTEx pipeline base directory"
)
parser.add_argument("outdir", type=Path, help="Output directory path")
parser.add_argument("tissues", nargs="+", type=str, help="Tissues to include")
parser.add_argument(
    "--keep_existing",
    action="store_true",
    default=False,
    help="Do not rerun a script if all its output files already exist",
)
args = parser.parse_args()

args.outdir.mkdir(exist_ok=True)
(args.outdir / "cis_pvals").mkdir(exist_ok=True)
(args.outdir / "covar").mkdir(exist_ok=True)
(args.outdir / "eqtl").mkdir(exist_ok=True)
(args.outdir / "expr").mkdir(exist_ok=True)
(args.outdir / "ref").mkdir(exist_ok=True)

files = ["top_assoc.txt", "eqtls_indep.txt"]
files += [f"{tissue}.cis_qtl_signif.txt.gz" for tissue in args.tissues]
files += [f"{tissue}.trans_qtl_pairs.txt.gz" for tissue in args.tissues]
files = [args.outdir / "eqtl" / x for x in files]
if not (args.keep_existing and all([file.exists() for file in files])):
    print("Assembling eQTL files")
    script = Path(__file__).parent / "eqtl_files.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), str(args.outdir), *args.tissues],
        check=True,
    )

if not (args.keep_existing and (args.outdir / "tissueInfo.txt").exists()):
    print("Making tissue info table")
    script = Path(__file__).parent / "tissueInfo.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), str(args.outdir), *args.tissues],
        check=True,
    )

files = [args.outdir / x for x in ["gene.txt", "autocomplete.json"]]
if not (args.keep_existing and all([file.exists() for file in files])):
    print("Assembling gene info")
    script = Path(__file__).parent / "gene.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), str(args.outdir), *args.tissues],
        check=True,
    )

files = [
    args.outdir / x for x in ["medianGeneExpression.txt.gz", "topExpressedGene.txt"]
]
if not (args.keep_existing and all([file.exists() for file in files])):
    print("Summarizing gene expression")
    script = Path(__file__).parent / "medianGeneExpression.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), str(args.outdir), *args.tissues],
        check=True,
    )

infile = args.indir / "ref" / "Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
if not (args.keep_existing and (args.outdir / "exon.txt").exists()):
    print("Making exon table")
    script = Path(__file__).parent / "exon.py"
    subprocess.run(["python3", script, str(infile), str(args.outdir)], check=True)

if not (args.keep_existing and (args.outdir / "singleTissueEqtl.zip").exists()):
    print("Assembling all significant cis associations")
    script = Path(__file__).parent / "singleTissueEqtl.py"
    subprocess.run(
        ["python3", script, str(args.indir), str(args.outdir), *args.tissues],
        check=True,
    )

files = [args.outdir / "cis_pvals" / f"{tissue}.zip" for tissue in args.tissues]
if not (args.keep_existing and all([file.exists() for file in files])):
    print("Extracting all cis p-values")
    script = Path(__file__).parent / "all_cis_pvals.py"
    for tissue in args.tissues:
        print(tissue)
        subprocess.run(
            ["python3", script, str(args.indir), str(args.outdir), tissue], check=True
        )

files = ["RatGTEx_rats.tsv", "RatGTEx_samples.tsv", "rats.html", "samples.html"]
files = [args.outdir / "ref" / x for x in files]
if not (args.keep_existing and all([file.exists() for file in files])):
    print("Making sample and rat info tables")
    script = Path(__file__).parent / "sample_table.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), str(args.outdir), *args.tissues],
        check=True,
    )

## Copy existing expression files
infiles = []
outfiles = []
for tissue in args.tissues:
    for units in {"iqn", "log2", "tpm"}:
        infiles.append(args.indir / tissue / f"{tissue}.expr.{units}.bed.gz")
        outfiles.append(args.outdir / "expr" / f"{tissue}.expr.{units}.bed.gz")
if not (args.keep_existing and all([file.exists() for file in outfiles])):
    print("Copying expression files")
    for infile in infiles:
        subprocess.run(["cp", str(infile), args.outdir / "expr"], check=True)

## Copy existing covariate files
outfiles = []
for tissue in args.tissues:
    outfiles.append(args.outdir / "covar" / f"{tissue}.covar.txt")
if not (args.keep_existing and all([file.exists() for file in outfiles])):
    print("Copying covariate files")
    for tissue in args.tissues:
        infile = args.indir / tissue / "covar.txt"
        outfile = args.outdir / "covar" / f"{tissue}.covar.txt"
        subprocess.run(["cp", str(infile), str(outfile)], check=True)

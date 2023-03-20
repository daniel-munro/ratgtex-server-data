"""Call scripts to process RatGTEx pipeline results into data files for the RatGTEx site"""

import argparse
from pathlib import Path
import subprocess
import yaml

parser = argparse.ArgumentParser(
    description="Process RatGTEx pipeline results into data files for the site"
)
parser.add_argument(
    "-i", "--indir", type=Path, required=True, help="Path to the RatGTEx pipeline base directory"
)
parser.add_argument(
    "-v", "--version", type=str, required=True, help="rn6 or rn7. All tissues in the config file that have directories inside '{version}/' will be processed."
)
parser.add_argument("-o", "--outdir", type=Path, required=True, help="Output directory path")
parser.add_argument(
    "--keep_existing",
    action="store_true",
    default=False,
    help="Do not rerun a script if all its output files already exist",
)
args = parser.parse_args()

config = yaml.safe_load(open(args.indir / "config.yaml"))
rn = args.version
tissues = [tissue for tissue in config if (args.indir / rn / tissue).exists()]

args.outdir.mkdir(exist_ok=True)
(args.outdir / "cis_pvals").mkdir(exist_ok=True)
(args.outdir / "covar").mkdir(exist_ok=True)
(args.outdir / "eqtl").mkdir(exist_ok=True)
(args.outdir / "expr").mkdir(exist_ok=True)
(args.outdir / "fastq_map").mkdir(exist_ok=True)
(args.outdir / "geno").mkdir(exist_ok=True)
(args.outdir / "rat_ids").mkdir(exist_ok=True)
(args.outdir / "ref").mkdir(exist_ok=True)
(args.outdir / "splice").mkdir(exist_ok=True)

## eQTL files
outfiles = [f"{rn}.top_assoc.txt", f"{rn}.eqtls_indep.txt"]
outfiles += [f"{tissue}.{rn}.cis_qtl_signif.txt.gz" for tissue in tissues]
outfiles += [f"{tissue}.{rn}.trans_qtl_pairs.txt.gz" for tissue in tissues]
outfiles = [args.outdir / "eqtl" / x for x in outfiles]
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Assembling eQTL files")
    script = Path(__file__).parent / "scripts" / "eqtl_files.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), rn, str(args.outdir), *tissues],
        check=True,
    )

## Splice files
outfiles = [f"{rn}.top_assoc_splice.txt", f"{rn}.sqtls_indep.txt"]
outfiles += [f"{tissue}.{rn}.splice.trans_qtl_pairs.txt.gz" for tissue in tissues]
outfiles = [args.outdir / "splice" / x for x in outfiles]
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Assembling sQTL files")
    script = Path(__file__).parent / "scripts" / "sqtl_files.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), rn, str(args.outdir), *tissues],
        check=True,
    )
outfiles = [args.outdir / "splice" / f"{tissue}.{rn}.leafcutter.bed.gz" for tissue in tissues]
outfiles += [args.outdir / "splice" / f"{tissue}.{rn}.covar_splice.txt" for tissue in tissues]
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Copying splicing phenotype and covariate files")
    for tissue in tissues:
        for suffix in ["leafcutter.bed.gz", "covar_splice.txt"]:
            infile = args.indir / rn / tissue / "splice" / f"{tissue}.{suffix}"
            outfile = args.outdir / "splice" / f"{tissue}.{rn}.{suffix}"
            subprocess.run(["cp", str(infile), str(outfile)], check=True)

## Tissue info table
if not (args.keep_existing and (args.outdir / f"{rn}.tissueInfo.txt").exists()):
    print("Making tissue info table")
    script = Path(__file__).parent / "scripts" / "tissueInfo.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), rn, str(args.outdir), *tissues],
        check=True,
    )

## Gene info table
outfiles = [args.outdir / x for x in [f"{rn}.gene.txt", f"{rn}.autocomplete.json"]]
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Assembling gene info")
    script = Path(__file__).parent / "scripts" / "gene.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), rn, str(args.outdir), *tissues],
        check=True,
    )

## Gene expression files
outfiles = [
    args.outdir / x for x in [f"{rn}.medianGeneExpression.txt.gz", f"{rn}.topExpressedGene.txt"]
]
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Summarizing gene expression")
    script = Path(__file__).parent / "scripts" / "medianGeneExpression.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), rn, str(args.outdir), *tissues],
        check=True,
    )

## Exon table
infile = {
    "rn6": args.indir / "ref_rn6" / "Rattus_norvegicus.Rnor_6.0.99.genes.gtf",
    "rn7": args.indir / "ref_rn7" / "Rattus_norvegicus.mRatBN7.2.108.genes.gtf",
}[rn]
outfile = args.outdir / "{rn}.exon.txt"
if not (args.keep_existing and outfile.exists()):
    print("Making exon table")
    script = Path(__file__).parent / "scripts" / "exon.py"
    subprocess.run(["python3", script, str(infile), str(outfile)], check=True)

## All significant cis-eQTL SNPs
if not (args.keep_existing and (args.outdir / f"{rn}.singleTissueEqtl.zip").exists()):
    print("Assembling all significant cis-eQTL associations")
    script = Path(__file__).parent / "scripts" / "singleTissueEqtl.py"
    subprocess.run(
        ["python3", script, str(args.indir), rn, str(args.outdir), *tissues],
        check=True,
    )

## All cis-window p-values
outfiles = [args.outdir / "cis_pvals" / f"{tissue}.{rn}.zip" for tissue in tissues]
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Extracting all cis p-values")
    script = Path(__file__).parent / "scripts" / "all_cis_pvals.py"
    for tissue in tissues:
        print(tissue)
        subprocess.run(
            ["python3", script, str(args.indir), rn, str(args.outdir), tissue], check=True
        )

## Rat and sample info tables
outfiles = [f"{rn}.RatGTEx_rats.tsv", f"{rn}.RatGTEx_samples.tsv", f"{rn}.rats.html", f"{rn}.samples.html"]
outfiles = [args.outdir / "ref" / x for x in outfiles]
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Making sample and rat info tables")
    script = Path(__file__).parent / "scripts" / "sample_table.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), rn, str(args.outdir), *tissues],
        check=True,
    )

## Copy existing expression files
infiles = []
outfiles = []
for tissue in tissues:
    for units in ["iqn", "log2", "tpm"]:
        infiles.append(args.indir / rn / tissue / f"{tissue}.expr.{units}.bed.gz")
        outfiles.append(args.outdir / "expr" / f"{tissue}.{rn}.expr.{units}.bed.gz")
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Copying expression files")
    for infile in infiles:
        subprocess.run(["cp", str(infile), args.outdir / "expr"], check=True)

## Copy existing covariate, fastq_map, and rat_ids files
outfiles = []
for tissue in tissues:
    outfiles.append(args.outdir / "covar" / f"{tissue}.{rn}.covar.txt")
    outfiles.append(args.outdir / "fastq_map" / f"{tissue}.{rn}.fastq_map.txt")
    outfiles.append(args.outdir / "rat_ids" / f"{tissue}.{rn}.rat_ids.txt")
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Copying covariate, fastq_map, and rat_ids files")
    for tissue in tissues:
        for ftype in ['covar', 'fastq_map', 'rat_ids']:
            infile = args.indir / rn / tissue / f"{ftype}.txt"
            outfile = args.outdir / ftype / f"{tissue}.{rn}.{ftype}.txt"
            subprocess.run(["cp", str(infile), str(outfile)], check=True)

## Copy genotype files
datasets = list(set(config[tissue]['geno_dataset'] for tissue in tissues))
outfiles = []
for dataset in datasets:
    for ext in ["vcf.gz", "vcf.gz.tbi"]:
        outfiles.append(args.outdir / "geno" / f"{dataset}.{rn}.{ext}")
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Copying genotype files")
    for dataset in datasets:
        for ext in ["vcf.gz", "vcf.gz.tbi"]:
            infile = args.indir / f"geno_{rn}" / f"{dataset}.{rn}.{ext}"
            outfile = args.outdir / "geno" / f"{dataset}.{rn}.{ext}"
            subprocess.run(["cp", str(infile), str(outfile)], check=True)

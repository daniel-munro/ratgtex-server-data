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
parser.add_argument("-o", "--outdir", type=Path, required=True, help="Output directory path")
parser.add_argument(
    "--keep_existing",
    action="store_true",
    default=False,
    help="Do not rerun a script if all its output files already exist",
)
args = parser.parse_args()

config = yaml.safe_load(open(args.indir / "config.yaml"))
version = "v3"
rn = "rn7"
v = "v3_rn7"
tissues = [tissue for tissue in config["tissues"] if (args.indir / version / tissue).exists()]

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
outfiles = [f"top_assoc.{v}.txt", f"eqtls_indep.{v}.txt"]
outfiles += [f"cis_qtl_signif.{tissue}.{v}.txt.gz" for tissue in tissues]
outfiles += [f"trans_qtl_pairs.{tissue}.{v}.txt.gz" for tissue in tissues]
outfiles = [args.outdir / "eqtl" / x for x in outfiles]
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Assembling eQTL files", flush=True)
    script = Path(__file__).parent / "scripts" / "eqtl_files.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), rn, str(args.outdir), *tissues],
        check=True,
    )

## Splice files
outfiles = [f"top_assoc_splice.{v}.txt", f"sqtls_indep.{v}.txt"]
outfiles += [f"trans_qtl_pairs.{tissue}.{v}.txt.gz" for tissue in tissues]
outfiles = [args.outdir / "splice" / x for x in outfiles]
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Assembling sQTL files", flush=True)
    script = Path(__file__).parent / "scripts" / "sqtl_files.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), rn, str(args.outdir), *tissues],
        check=True,
    )
outfiles = [args.outdir / "splice" / f"leafcutter.{tissue}.{v}.bed.gz" for tissue in tissues]
outfiles += [args.outdir / "splice" / f"covar_splice.{tissue}.{v}.txt" for tissue in tissues]
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Copying splicing phenotype and covariate files", flush=True)
    for tissue in tissues:
        for suffix in ["leafcutter.bed.gz", "covar_splice.txt"]:
            infile = args.indir / rn / tissue / "splice" / f"{tissue}.{suffix}"
            outfile = args.outdir / "splice" / f"{tissue}.{v}.{suffix}"
            subprocess.run(["cp", str(infile), str(outfile)], check=True)

## Tissue info table
if not (args.keep_existing and (args.outdir / f"tissueInfo.{v}.txt").exists()):
    print("Making tissue info table", flush=True)
    script = Path(__file__).parent / "scripts" / "tissueInfo.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), rn, str(args.outdir), *tissues],
        check=True,
    )

## Gene info table
outfiles = [args.outdir / x for x in [f"gene.{v}.txt", f"autocomplete.{v}.json"]]
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Assembling gene info", flush=True)
    script = Path(__file__).parent / "scripts" / "gene.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), rn, str(args.outdir), *tissues],
        check=True,
    )

## Gene expression files
outfiles = [
    args.outdir / x for x in [f"medianGeneExpression.{v}.txt.gz", f"topExpressedGene.{v}.txt"]
]
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Summarizing gene expression", flush=True)
    script = Path(__file__).parent / "scripts" / "medianGeneExpression.R"
    subprocess.run(
        ["Rscript", script, str(args.indir), rn, str(args.outdir), *tissues],
        check=True,
    )

## Exon table
infile = args.indir / "ref" / "GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf"
outfile = args.outdir / f"exon.{v}.txt"
if not (args.keep_existing and outfile.exists()):
    print("Making exon table", flush=True)
    script = Path(__file__).parent / "scripts" / "exon.py"
    subprocess.run(["python3", script, str(infile), str(outfile)], check=True)

## All significant cis-eQTL SNPs
if not (args.keep_existing and (args.outdir / f"singleTissueEqtl.{v}.zip").exists()):
    print("Assembling all significant cis-eQTL associations", flush=True)
    script = Path(__file__).parent / "scripts" / "singleTissueEqtl.py"
    subprocess.run(
        ["python3", script, str(args.indir), rn, str(args.outdir), *tissues],
        check=True,
    )

## All cis-window p-values
outfiles = [args.outdir / "cis_pvals" / f"{tissue}.{v}.zip" for tissue in tissues]
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Extracting all cis p-values", flush=True)
    script = Path(__file__).parent / "scripts" / "all_cis_pvals.py"
    for tissue in tissues:
        print(tissue, flush=True)
        subprocess.run(
            ["python3", script, str(args.indir), rn, str(args.outdir), tissue], check=True
        )

## Rat and sample info tables
outfiles = [f"RatGTEx_rats.{v}.tsv", f"RatGTEx_samples.{v}.tsv", f"rats.{v}.html", f"samples.{v}.html"]
outfiles = [args.outdir / "ref" / x for x in outfiles]
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Making sample and rat info tables", flush=True)
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
        outfiles.append(args.outdir / "expr" / f"expr.{units}.{tissue}.{v}.bed.gz")
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Copying expression files", flush=True)
    for tissue in tissues:
        for units in ["iqn", "log2", "tpm"]:
            infile = args.indir / rn / tissue / f"{tissue}.expr.{units}.bed.gz"
            outfile = args.outdir / "expr" / f"expr.{units}.{tissue}.{v}.bed.gz"
            subprocess.run(["cp", str(infile), str(outfile)], check=True)

## Copy existing covariate, fastq_map, and rat_ids files
outfiles = []
for tissue in tissues:
    outfiles.append(args.outdir / "covar" / f"covar.{tissue}.{v}.txt")
    outfiles.append(args.outdir / "fastq_map" / f"fastq_map.{tissue}.{v}.txt")
    outfiles.append(args.outdir / "rat_ids" / f"rat_ids.{tissue}.{v}.txt")
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Copying covariate, fastq_map, and rat_ids files", flush=True)
    for tissue in tissues:
        for ftype in ['covar', 'fastq_map', 'rat_ids']:
            infile = args.indir / rn / tissue / f"{ftype}.txt"
            outfile = args.outdir / ftype / f"{ftype}.{tissue}.{v}.txt"
            subprocess.run(["cp", str(infile), str(outfile)], check=True)

## Copy genotype files
datasets = list(set(config[tissue]['geno_dataset'] for tissue in tissues))
outfiles = []
for dataset in datasets:
    for ext in ["vcf.gz", "vcf.gz.tbi"]:
        outfiles.append(args.outdir / "geno" / f"{dataset}.{v}.{ext}")
if not (args.keep_existing and all(file.exists() for file in outfiles)):
    print("Copying genotype files", flush=True)
    for dataset in datasets:
        for ext in ["vcf.gz", "vcf.gz.tbi"]:
            infile = args.indir / f"geno_{rn}" / f"{dataset}.{ext}"
            outfile = args.outdir / "geno" / f"{dataset}.{v}.{ext}"
            subprocess.run(["cp", str(infile), str(outfile)], check=True)

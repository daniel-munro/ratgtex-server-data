"""Check if all required data files from the RatGTEx pipeline exist."""

import argparse
from pathlib import Path
import pandas as pd
import yaml

parser = argparse.ArgumentParser(
    description="Check if all required data files from the RatGTEx pipeline exist."
)
parser.add_argument(
    "indir", type=Path, help="Path to the RatGTEx pipeline base directory"
)
parser.add_argument(
    "version", type=str, help="rn6 or rn7. Only checks the files for this version."
)
args = parser.parse_args()

rn = args.version
config = yaml.safe_load(open(args.indir / "config.yaml"))
tissues = [tissue for tissue in config if (args.indir / rn / tissue).exists()]

files = []

anno = {
    "rn6": args.indir / "ref_rn6" / "Rattus_norvegicus.Rnor_6.0.99.genes.gtf",
    "rn7": args.indir / "ref_rn7" / "Rattus_norvegicus.mRatBN7.2.108.genes.gtf",
}[rn]
files.append(anno)

general = [
    args.indir / f"geno_{rn}" / "alleles.txt.gz",
    args.indir / f"geno_{rn}" / "genotyping_log.csv",
    args.indir / "tissue_info.txt",
    args.indir / f"ref_{rn}" / "GENES_RAT.txt",
]
files += general

accessions = ["GSE201236", "GSE173137", "GSE173138", "GSE173136", "GSE173140", "GSE173139"]
files += [args.indir / "samples" / f"{accession}_series_matrix.txt.gz" for accession in accessions]

for tissue in tissues:
    fnames = [
        "covar.txt",
        "fastq_map.txt",
        "rat_ids.txt",
        f"{tissue}.aFC.txt",
        f"{tissue}.cis_independent_qtl.txt.gz",
        f"{tissue}.cis_qtl_all_pvals.txt.gz",
        f"{tissue}.cis_qtl_signif.txt.gz",
        f"{tissue}.cis_qtl.txt.gz",
        f"{tissue}.expr.iqn.bed.gz",
        f"{tissue}.expr.log2.bed.gz",
        f"{tissue}.expr.tpm.bed.gz",
        f"{tissue}.trans_qtl_pairs.txt.gz",
    ]
    files += [args.indir / rn / tissue / fname for fname in fnames]
    fnames = [
        f"{tissue}.leafcutter.bed.gz",
        f"{tissue}.covar_splice.txt",
        f"{tissue}_splice.cis_qtl.txt.gz",
        f"{tissue}_splice.cis_independent_qtl.txt.gz",
        f"{tissue}_splice.trans_qtl_pairs.txt.gz",
    ]
    files += [args.indir / rn / tissue / "splice" / fname for fname in fnames]

datasets = list(set(config[tissue]['geno_dataset'] for tissue in tissues))
for ext in ["vcf.gz", "vcf.gz.tbi"]:
    files += [args.indir / f"geno_{rn}" / f"{dataset}.{ext}" for dataset in datasets]

print("Tissues that will be included:")
print(" ".join(tissues))

any_missing = False
for file in files:
    if not file.exists():
        print(f"{file} is missing")
        any_missing = True

if not any_missing:
    print("All necessary files are present.")
    print("Remember to update site HTML files as needed, including copying data from rats.html and samples.html into about/samples/index.html.")

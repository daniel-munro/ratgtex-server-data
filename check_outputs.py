"""Check if all required data files for the RatGTEx site have been produced or copied."""

import argparse
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser(
    description="Check if all required data files for the RatGTEx site have been produced or copied."
)
parser.add_argument("outdir", type=Path, help="Path of the site data directory")
parser.add_argument(
    "version", type=str, help="rn6 or rn7. Only checks the files for this version."
)
args = parser.parse_args()

rn = args.version
tissue_info = pd.read_csv(args.outdir / f"{rn}.tissueInfo.txt", sep="\t")
tissues = list(tissue_info["tissueSiteDetailId"])

files = []

## Gene files
gene = [f"{rn}.gene.txt", f"{rn}.autocomplete.json"]
gene += [f"{rn}.medianGeneExpression.txt.gz", f"{rn}.topExpressedGene.txt"]
gene += [f"{rn}.exon.txt"]
files += [args.outdir / x for x in gene]

## eQTL files
eqtl = [f"{rn}.top_assoc.txt", f"{rn}.eqtls_indep.txt"]
eqtl += [f"{tissue}.{rn}.cis_qtl_signif.txt.gz" for tissue in tissues]
eqtl += [f"{tissue}.{rn}.trans_qtl_pairs.txt.gz" for tissue in tissues]
files += [args.outdir / "eqtl" / x for x in eqtl]

files += [args.outdir / f"{rn}.singleTissueEqtl.zip"]
files += [args.outdir / "cis_pvals" / f"{tissue}.{rn}.zip" for tissue in tissues]

## Sample and rat tables
ref = [f"{rn}.RatGTEx_rats.tsv", f"{rn}.RatGTEx_samples.tsv"]
files += [args.outdir / "ref" / x for x in ref]

## Expression files
for tissue in tissues:
    for units in ["iqn", "log2", "tpm"]:
        files.append(args.outdir / "expr" / f"{tissue}.{rn}.expr.{units}.bed.gz")

## Covariate, fastq_map, and rat_ids files
for tissue in tissues:
    files.append(args.outdir / "covar" / f"{tissue}.{rn}.covar.txt")
    files.append(args.outdir / "fastq_map" / f"{tissue}.{rn}.fastq_map.txt")
    files.append(args.outdir / "rat_ids" / f"{tissue}.{rn}.rat_ids.txt")

## Genotype files
geno = []
for dset in tissue_info["dataset"].unique():
    geno += [f"{dset}.{rn}.vcf.gz", f"{dset}.{rn}.vcf.gz.tbi"]
files += [args.outdir / "geno" / x for x in geno]

## Splice phenotypes, covariates, and sQTL files
splice = [f"{rn}.top_assoc_splice.txt", f"{rn}.sqtls_indep.txt"]
splice += [f"{tissue}.{rn}.leafcutter.bed.gz" for tissue in tissues]
splice += [f"{tissue}.{rn}.covar_splice.txt" for tissue in tissues]
splice += [f"{tissue}.{rn}.splice.cis_qtl_signif.txt.gz" for tissue in tissues]
splice += [f"{tissue}.{rn}.splice.trans_qtl_pairs.txt.gz" for tissue in tissues]
files += [args.outdir / "splice" / x for x in splice]

any_missing = False
for file in files:
    if not file.exists():
        print(f"{file} is missing")
        any_missing = True

if not any_missing:
    print("All necessary files are present.")

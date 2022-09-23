"""Check if all required data files for the RatGTEx site have been produced or copied."""

import argparse
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser(
    description="Check if all required data files for the RatGTEx site have been produced or copied."
)
parser.add_argument("outdir", type=Path, help="Path of the site data directory")
args = parser.parse_args()

tissue_info = pd.read_csv(args.outdir / "tissueInfo.txt", sep="\t")
tissues = list(tissue_info["tissueSiteDetailId"])

files = []

## Gene files
gene = ["gene.txt", "autocomplete.json"]
gene += ["medianGeneExpression.txt.gz", "topExpressedGene.txt"]
gene += ["exon.txt"]
files += [args.outdir / x for x in gene]

## eQTL files
eqtl = ["top_assoc.txt", "eqtls_indep.txt"]
eqtl += [f"{tissue}.cis_qtl_signif.txt.gz" for tissue in tissues]
eqtl += [f"{tissue}.trans_qtl_pairs.txt.gz" for tissue in tissues]
files += [args.outdir / "eqtl" / x for x in eqtl]

files += [args.outdir / "singleTissueEqtl.zip"]
files += [args.outdir / "cis_pvals" / f"{tissue}.zip" for tissue in tissues]

## Sample and rat tables
ref = ["RatGTEx_rats.tsv", "RatGTEx_samples.tsv"]
files += [args.outdir / "ref" / x for x in ref]

## Expression files
for tissue in tissues:
    for units in ["iqn", "log2", "tpm"]:
        files.append(args.outdir / "expr" / f"{tissue}.expr.{units}.bed.gz")

## Covariate, fastq_map, and rat_ids files
for tissue in tissues:
    files.append(args.outdir / "covar" / f"{tissue}.covar.txt")
    files.append(args.outdir / "fastq_map" / f"{tissue}.fastq_map.txt")
    files.append(args.outdir / "rat_ids" / f"{tissue}.rat_ids.txt")

## Genotype files
geno = []
for dset in tissue_info["dataset"].unique():
    geno += [f"{dset}.vcf.gz", f"{dset}.vcf.gz.tbi"]
files += [args.outdir / "geno" / x for x in geno]

## Splice phenotypes, covariates, and sQTL files
splice = ["top_assoc_splice.txt", "sqtls_indep.txt"]
splice += [f"{tissue}.leafcutter.bed.gz" for tissue in tissues]
splice += [f"{tissue}.covar_splice.txt" for tissue in tissues]
splice += [f"{tissue}_splice.trans_qtl_pairs.txt.gz" for tissue in tissues]
files += [args.outdir / "splice" / x for x in splice]

any_missing = False
for file in files:
    if not file.exists():
        print(f"{file} is missing")
        any_missing = True

if not any_missing:
    print("All necessary files are present.")
    print("Remember to update site HTML files as needed, including copying data from rats.html and samples.html into about/samples/index.html.")

"""Check if all required data files for the RatGTEx site have been produced or copied."""

import argparse
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser(
    description="Check if all required data files for the RatGTEx site have been produced or copied."
)
parser.add_argument("outdir", type=Path, help="Path of the site data directory")
args = parser.parse_args()

v = "v3_rn7"

tissue_info = pd.read_csv(args.outdir / f"tissueInfo.{v}.txt", sep="\t")
tissues = list(tissue_info["tissueSiteDetailId"])

files = []

## Gene files
gene = [f"gene.{v}.txt", f"autocomplete.{v}.json"]
gene += [f"medianGeneExpression.{v}.txt.gz", f"topExpressedGene.{v}.txt"]
gene += [f"exon.{v}.txt"]
files += [args.outdir / x for x in gene]

## eQTL files
eqtl = [f"top_assoc.{v}.txt", f"eqtls_indep.{v}.txt"]
eqtl += [f"cis_qtl_signif.{tissue}.{v}.txt.gz" for tissue in tissues]
eqtl += [f"trans_qtl_pairs.{tissue}.{v}.txt.gz" for tissue in tissues]
files += [args.outdir / "eqtl" / x for x in eqtl]

files += [args.outdir / f"singleTissueEqtl.{v}.zip"]
files += [args.outdir / "cis_pvals" / f"{tissue}.{v}.zip" for tissue in tissues]

## Sample and rat tables
ref = [f"RatGTEx_rats.{v}.tsv", f"RatGTEx_samples.{v}.tsv"]
files += [args.outdir / "ref" / x for x in ref]

## Expression files
for tissue in tissues:
    for units in ["iqn", "log2", "tpm"]:
        files.append(args.outdir / "expr" / f"expr.{units}.{tissue}.{v}.bed.gz")

## Covariate, fastq_map, and rat_ids files
for tissue in tissues:
    files.append(args.outdir / "covar" / f"covar.{tissue}.{v}.txt")
    files.append(args.outdir / "fastq_map" / f"fastq_map.{tissue}.{v}.txt")
    files.append(args.outdir / "rat_ids" / f"rat_ids.{tissue}.{v}.txt")

## Genotype files
geno = []
for dset in tissue_info["dataset"].unique():
    geno += [f"{dset}.{v}.vcf.gz", f"{dset}.{v}.vcf.gz.tbi"]
files += [args.outdir / "geno" / x for x in geno]

## Splice phenotypes, covariates, and sQTL files
splice = [f"top_assoc_splice.{v}.txt", f"sqtls_indep.{v}.txt"]
splice += [f"leafcutter.{tissue}.{v}.bed.gz" for tissue in tissues]
splice += [f"covar_splice.{tissue}.{v}.txt" for tissue in tissues]
splice += [f"splice.cis_qtl_signif.{tissue}.{v}.txt.gz" for tissue in tissues]
splice += [f"splice.trans_qtl_pairs.{tissue}.{v}.txt.gz" for tissue in tissues]
files += [args.outdir / "splice" / x for x in splice]

any_missing = False
for file in files:
    if not file.exists():
        print(f"{file} is missing")
        any_missing = True

if not any_missing:
    print("All necessary files are present.")

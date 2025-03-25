"""Assemble all significant cis associations in zip archive of per-gene files"""

import argparse
from pathlib import Path
import zipfile
import pandas as pd

parser = argparse.ArgumentParser(
    description="Assemble all significant cis associations in zip archive of per-gene files"
)
parser.add_argument(
    "indir", type=Path, help="Path to the RatGTEx pipeline base directory"
)
parser.add_argument("version", type=str, help="e.g. v3")
parser.add_argument("outdir", type=Path, help="Output directory path")
parser.add_argument("tissues", nargs="+", type=str, help="Tissues to include")
args = parser.parse_args()

sig = []
for tissue in args.tissues:
    print(tissue)
    d = pd.read_csv(args.indir / args.version / tissue / f"{tissue}.cis_qtl_signif.txt.gz", sep="\t")
    d = d[["phenotype_id", "variant_id", "pval_nominal", "slope"]]
    d = d.rename(
        columns={
            "variant_id": "variantId",
            "pval_nominal": "pValue",
            "slope": "nes",
        }
    )
    # d["chromosome"] = d["variantId"].str.extract(r"chr(\d+):\d+$")
    # d["pos"] = d["variantId"].str.extract(r"chr\d+:(\d+)$")
    d[["chromosome", "pos"]] = d["variantId"].str.split(":", expand=True)
    d["chromosome"] = d["chromosome"].str.replace("chr", "")
    d["tissueSiteDetailId"] = tissue
    sig.append(d)

sig = pd.concat(sig)

outfile = args.outdir / f"{args.version}.singleTissueEqtl.zip"
with zipfile.ZipFile(outfile, 'w', zipfile.ZIP_DEFLATED) as out:
    for gene, d in sig.groupby("phenotype_id"):
        d = d.drop(columns=["phenotype_id"])
        d_str = d.to_csv(
            sep="\t",
            index=False,
            float_format="%g",
        )
        out.writestr(f"{args.version}.singleTissueEqtl/{gene}.txt", d_str)

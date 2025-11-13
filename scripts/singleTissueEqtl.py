"""Assemble all significant cis associations in zip archive of per-gene files

Inputs:
    {indir}/v4/{tissue}/{tissue}.cis_qtl_signif.txt.gz

Outputs:
    {outdir}/singleTissueEqtl.v4.zip
"""

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
parser.add_argument("version", type=str, help="e.g. v4")
parser.add_argument("outdir", type=Path, help="Output directory path for all server data")
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
    d[["chromosome", "pos"]] = d["variantId"].str.split(":", expand=True)
    d["tissueSiteDetailId"] = tissue
    sig.append(d)

sig = pd.concat(sig)

outfile = args.outdir / "eqtl" / f"singleTissueEqtl.{args.version}.zip"
with zipfile.ZipFile(outfile, 'w', zipfile.ZIP_DEFLATED) as out:
    for gene, d in sig.groupby("phenotype_id"):
        d = d.drop(columns=["phenotype_id"])
        d_str = d.to_csv(
            sep="\t",
            index=False,
            float_format="%g",
        )
        out.writestr(f"singleTissueEqtl.{args.version}/{gene}.txt", d_str)

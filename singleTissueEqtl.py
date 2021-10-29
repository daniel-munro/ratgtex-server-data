import pandas as pd

tissues = ["Eye", "IL", "LHb", "NAcc", "OFC", "PL"]

# fname = {tis: f"/Users/dan/br/data/tensorqtl/{tis[0]}QCT.cis_qtl_signif.txt.gz" for tis in tissues[:5]}
# fname["Eye"] = "/Users/dan/eye/data/tensorqtl/main5.cis_qtl_signif.txt.gz"

sig = []
for tissue in tissues:
    print(tissue)
    # d = pd.read_csv(fname[tissue], sep="\t")
    d = pd.read_csv(f"/Users/dan/ratgtex/data/tensorqtl/{tissue}.cis_qtl_signif.txt.gz", sep="\t")
    d = d[["phenotype_id", "variant_id", "pval_nominal", "slope"]]
    d = d.rename(
        columns={
            # "phenotype_id": "geneId",
            "variant_id": "variantId",
            "pval_nominal": "pValue",
            "slope": "nes",
        }
    )
    d["chromosome"] = d["variantId"].str.extract(r"chr(\d+):\d+$")
    d["pos"] = d["variantId"].str.extract(r"chr\d+:(\d+)$")
    d["tissueSiteDetailId"] = tissue
    sig.append(d)
    # d.to_csv(f"/Users/dan/br/portal_server/data/singleTissueEqtl/{tissue}.txt.gz", sep="\t", index=False, float_format="%g")

sig = pd.concat(sig)
# print("saving...")
# sig.to_csv("/Users/dan/br/portal_server/data/singleTissueEqtl.txt.gz", sep="\t", index=False, float_format="%g")

for gene, d in sig.groupby("phenotype_id"):
    d = d.drop(columns=["phenotype_id"])
    d.to_csv(f"/Users/dan/code/portal_data/singleTissueEqtl/{gene}.txt", sep="\t", index=False, float_format="%g")

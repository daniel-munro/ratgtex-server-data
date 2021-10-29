from gtfparse import read_gtf

d = read_gtf("/Users/dan/br/data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf")
d = d.loc[d["feature"] == "exon", :]
# Seems to be a formatting problem with this GTF file:
d["exon_number"] = d["exon_number"].str.replace('"', '').astype(int)
d = d.rename(
    columns={
        "seqname": "chromosome",
        "exon_id": "exonId",
        "exon_number": "exonNumber",
        "feature": "featureType",
        "gene_id": "geneId",
        "gene_name": "geneSymbol",
        "transcript_id": "transcriptId",
    }
)
d = d[
    [
        "chromosome",
        "end",
        "exonId",
        "exonNumber",
        "featureType",
        "geneId",
        "geneSymbol",
        "start",
        "strand",
        "transcriptId",
    ]
]
d.to_csv("/Users/dan/br/portal_server/data/exon.txt", sep="\t", index=False)

import argparse
from gtfparse import read_gtf
from pathlib import Path

parser = argparse.ArgumentParser(description='Prepare exon annotations for server API')
parser.add_argument('input', type=str, help='GTF gene annotation file')
parser.add_argument('outdir', type=Path, help='Output directory path')
args = parser.parse_args()

d = read_gtf(args.input)
d = d.loc[d['feature'] == 'exon', :]
# Seems to be a formatting problem with this GTF file:
d['exon_number'] = d['exon_number'].str.replace('"', '').astype(int)
d = d.rename(
    columns={
        'seqname': 'chromosome',
        'exon_id': 'exonId',
        'exon_number': 'exonNumber',
        'feature': 'featureType',
        'gene_id': 'geneId',
        'gene_name': 'geneSymbol',
        'transcript_id': 'transcriptId',
    }
)
d = d[
    [
        'chromosome',
        'end',
        'exonId',
        'exonNumber',
        'featureType',
        'geneId',
        'geneSymbol',
        'start',
        'strand',
        'transcriptId',
    ]
]
d.to_csv(args.outdir / 'exon.txt', sep='\t', index=False)

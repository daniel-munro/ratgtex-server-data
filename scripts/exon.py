import argparse
from gtfparse import read_gtf
from pathlib import Path

parser = argparse.ArgumentParser(description='Prepare exon annotations for server API')
parser.add_argument('input', type=Path, help='GTF gene annotation file')
parser.add_argument('output', type=Path, help='Output file name (TXT/TSV)')
args = parser.parse_args()

d = read_gtf(args.input)
# Newer versions return a polars DF by default, but not all versions allow
# return type to be specified, so this handles older and newer versions:
if type(d).__module__ == 'polars.dataframe.frame':
    d = d.to_pandas()

d = d.loc[d['feature'] == 'exon', :]
# Seems to be a formatting problem with this GTF file:
d['exon_number'] = d['exon_number'].str.replace('"', '').str.replace("'", '').astype(int)
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
d.to_csv(args.output, sep='\t', index=False)

# ratgtex-server-data
Process RatGTEx pipeline results into data files for the RatGTEx site

Once the [RatGTEx pipeline](https://github.com/daniel-munro/ratgtex-pipeline) is run for all tissues, this code processes the results and assembles them into data files for [ratgtex.org](https://ratgtex.org). Most of them are put into a convenient format to be used by the [API](github.com/daniel-munro/ratgtex), which is used mainly for the data visualizations.

## Inputs (in specified input directory):

- `geno/alleles.txt.gz`: Three-column tab-delimited file with variant ID, REF, and ALT alleles for the union of SNPs extracted from the genotype VCF files.
- `tissue_info.txt`: Tab-delimited table of basic information about each tissue, e.g. full name and color to use in website visualizations.
- `ref/GENES_RAT.txt`: Gene information downloaded from the Rat Genome Database. Used to show descriptive gene names.
- `samples/GSE173141_series_matrix.txt.gz`: Downloaded GEO series metadata for brain regions that has been edited. More files will need to be added as inputs to `sample_table.R` as datasets are added to GEO.
- `ref/Rattus_norvegicus.Rnor_6.0.99.genes.bed`
- `ref/Rattus_norvegicus.Rnor_6.0.99.genes.gtf`
- `{tissue}/rat_ids.txt`
- `{tissue}/{tissue}.aFC.txt`
- `{tissue}/{tissue}.cis_independent_qtl.txt.gz`
- `{tissue}/{tissue}.cis_qtl_all_pvals.txt.gz`
- `{tissue}/{tissue}.cis_qtl_signif.txt.gz`
- `{tissue}/{tissue}.cis_qtl.txt.gz`
- `{tissue}/{tissue}.expr.tpm.bed.gz`
- `{tissue}/{tissue}.trans_qtl_pairs.txt.gz`

## Outputs (in specified output directory):

- `autocomplete.json`
- `cis_pvals/{tissue}.zip`
- `eqtl/eqtls_indep.txt`
- `eqtl/top_assoc.txt`
- `eqtl/{tissue}.cis_qtl_signif.txt.gz`
- `eqtl/{tissue}.trans_qtl_pairs.txt.gz`
- `exon.txt`
- `gene.txt`
- `medianGeneExpression.txt.gz`
- `ref/RatGTEx_rats.tsv`
- `ref/RatGTEx_samples.tsv`
- `ref/rats.html`: This HTML table must be copied into /about/samples/index.html
- `ref/samples.html`: This HTML table must be copied into /about/samples/index.html
- `singleTissueEqtl.zip`
- `tissueInfo.txt`
- `topExpressedGene.txt`

## Requirements

- pandas (Python)
- gtfparse (Python)
- tidyverse (R)

## Usage

Run `run.py`, which calls R and Python scripts that prepare the server data. The parameters are:

1. Path to the RatGTEx pipeline base directory
2. Output directory path
3. The remaining parameters are the list of tissues to include

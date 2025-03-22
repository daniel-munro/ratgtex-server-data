# ratgtex-server-data
Process RatGTEx pipeline results into data files for the RatGTEx site

Once the [RatGTEx pipeline](https://github.com/daniel-munro/ratgtex-pipeline) is run for all tissues, this code processes the results and assembles them into data files for [ratgtex.org](https://ratgtex.org).
Many of them are put into a convenient format to be used by the [API](github.com/daniel-munro/ratgtex), which is used mainly for the data visualizations.
Others are provided for download.
This code is designed to run only the latest RatGTEx data version, but files are still labeled with RatGTEx version and genome version for clarity.

## Inputs (in specified input directory):

- `tissue_info.txt`: Tab-delimited table of basic information about each tissue, e.g. full name and color to use in website visualizations.
- `geno/{geno_dataset}.vcf.gz{,.tbi}`: Genotype files, each including the individuals for one or more tissues.
- `geno/alleles.txt.gz`: Three-column tab-delimited file with variant ID, REF, and ALT alleles for the union of SNPs extracted from the genotype VCF files.
- `geno/genotyping_log.csv`: CSV metadata file containing sex and coat color info to add to the rat info table.
- `samples/{accession}_series_matrix.txt.gz`: GEO series metadata for each tissue whose RNA-Seq data has been added to GEO. When a tissue is added to GEO, download this file for the series (not the superseries in case of multitissue datasets) and add the accession to `sample_table.R`.
- `ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf`
- `v3/{tissue}/covar.txt`
- `v3/{tissue}/fastq_map.txt`
- `v3/{tissue}/rat_ids.txt`
- `v3/{tissue}/{tissue}.aFC.txt`
- `v3/{tissue}/{tissue}.cis_independent_qtl.txt.gz`
- `v3/{tissue}/{tissue}.cis_qtl_all_pvals.txt.gz`
- `v3/{tissue}/{tissue}.cis_qtl_signif.txt.gz`
- `v3/{tissue}/{tissue}.cis_qtl.txt.gz`
- `v3/{tissue}/{tissue}.expr.{iqn,log2,tpm}.bed.gz`
- `v3/{tissue}/{tissue}.trans_qtl_pairs.txt.gz`
- `v3/{tissue}/splice/{tissue}.leafcutter.bed.gz`
- `v3/{tissue}/splice/{tissue}.covar_splice.txt`
- `v3/{tissue}/splice/{tissue}_splice.cis_qtl.txt.gz`
- `v3/{tissue}/splice/{tissue}_splice.cis_independent_qtl.txt.gz`
- `v3/{tissue}/splice/{tissue}_splice.trans_qtl_pairs.txt.gz`

## Outputs (in specified output directory):

- `geno/{geno_dataset}.v3_rn7.vcf.gz{,.tbi}`
- `autocomplete.v3_rn7.json`
- `exon.v3_rn7.txt`
- `gene.v3_rn7.txt`
- `medianGeneExpression.v3_rn7.txt.gz`
- `singleTissueEqtl.v3_rn7.zip`
- `tissueInfo.v3_rn7.txt`
- `topExpressedGene.v3_rn7.txt`
- `cis_pvals/{tissue}.v3_rn7.zip`
- `covar/covar.{tissue}.v3_rn7.txt`
- `eqtl/eqtls_indep.v3_rn7.txt`
- `eqtl/top_assoc.v3_rn7.txt`
- `eqtl/cis_qtl_signif.{tissue}.v3_rn7.txt.gz`
- `eqtl/trans_qtl_pairs.{tissue}.v3_rn7.txt.gz`
- `expr/expr.{iqn,log2,tpm}.{tissue}.v3_rn7.bed.gz`
- `fastq_map/fastq_map.{tissue}.v3_rn7.txt`
- `rat_ids/rat_ids.{tissue}.v3_rn7.txt`
- `ref/RatGTEx_rats.v3_rn7.tsv`
- `ref/RatGTEx_samples.v3_rn7.tsv`
- `ref/rats.v3_rn7.html`: Jekyll will insert this HTML table into the About Samples page.
- `ref/samples.v3_rn7.html`: Jekyll will insert this HTML table into the About Samples page.
- `splice/top_assoc_splice.v3_rn7.txt`
- `splice/sqtls_indep.v3_rn7.txt`
- `splice/leafcutter.{tissue}.v3_rn7.bed.gz"`
- `splice/covar_splice.{tissue}.v3_rn7.txt`
- `splice/splice.trans_qtl_pairs.{tissue}.v3_rn7.txt.gz`

## Requirements

- pandas (Python)
- gtfparse (Python)
- tidyverse (R)

## Usage

Run `check_inputs.py` to check if all necessary input files exist. If so, run `run.py`, which calls R and Python scripts that prepare the server data. The parameters are:

1. Path to the RatGTEx pipeline base directory
2. Output directory path

Then, do these additional steps:

1. Run `check_outputs.py` to check if all necessary files have been created/copied.
2. `{outdir}/studies/{tissues}/` contain data directly used in the original publications, which may or may not differ from the data from the unified RatGTEx pipeline. These aren't used in the portal except to be available for download. Copy any files to this location and add links to them on the original study data page.

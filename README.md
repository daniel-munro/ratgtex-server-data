# ratgtex-server-data
Process RatGTEx pipeline results into data files for the RatGTEx site

Once the [RatGTEx pipeline](https://github.com/daniel-munro/ratgtex-pipeline) is run for all tissues, this code processes the results and assembles them into data files for [ratgtex.org](https://ratgtex.org). Most of them are put into a convenient format to be used by the [API](github.com/daniel-munro/ratgtex), which is used mainly for the data visualizations. This should be run separately for each genome version, and all output file names will contain the version.

## Inputs (in specified input directory):

- `geno/alleles.txt.gz`: Three-column tab-delimited file with variant ID, REF, and ALT alleles for the union of SNPs extracted from the genotype VCF files.
- `geno/genotyping_log.csv`: CSV metadata file containing sex and coat color info to add to the rat info table.
- `tissue_info.txt`: Tab-delimited table of basic information about each tissue, e.g. full name and color to use in website visualizations.
- `ref/GENES_RAT.txt`: Gene information downloaded from the Rat Genome Database. Used to show descriptive gene names.
- `samples/{accession}_series_matrix.txt.gz`: GEO series metadata for each tissue whose RNA-Seq data has been added to GEO. When a tissue is added to GEO, download this file for the series (not the superseries in case of multitissue datasets) and add the accession to `sample_table.R`.
- `ref/Rattus_norvegicus.Rnor_6.0.99.genes.bed`
- `ref/Rattus_norvegicus.Rnor_6.0.99.genes.gtf`
- `{tissue}/covar.txt`
- `{tissue}/fastq_map.txt`
- `{tissue}/rat_ids.txt`
- `{tissue}/{tissue}.aFC.txt`
- `{tissue}/{tissue}.cis_independent_qtl.txt.gz`
- `{tissue}/{tissue}.cis_qtl_all_pvals.txt.gz`
- `{tissue}/{tissue}.cis_qtl_signif.txt.gz`
- `{tissue}/{tissue}.cis_qtl.txt.gz`
- `{tissue}/{tissue}.expr.iqn.bed.gz`
- `{tissue}/{tissue}.expr.log2.bed.gz`
- `{tissue}/{tissue}.expr.tpm.bed.gz`
- `{tissue}/{tissue}.trans_qtl_pairs.txt.gz`
- `{tissue}/splice/{tissue}.leafcutter.bed.gz`
- `{tissue}/splice/{tissue}.covar_splice.txt`
- `{tissue}/splice/{tissue}_splice.cis_qtl.txt.gz`
- `{tissue}/splice/{tissue}_splice.cis_independent_qtl.txt.gz`
- `{tissue}/splice/{tissue}_splice.trans_qtl_pairs.txt.gz`

## Outputs (in specified output directory):

- `autocomplete.json`
- `cis_pvals/{tissue}.zip`
- `covar/{tissue}.covar.txt`
- `eqtl/eqtls_indep.txt`
- `eqtl/top_assoc.txt`
- `eqtl/{tissue}.cis_qtl_signif.txt.gz`
- `eqtl/{tissue}.trans_qtl_pairs.txt.gz`
- `exon.txt`
- `expr/{tissue}.expr.{iqn,log2,tpm}.bed.gz`
- `fastq_map/{tissue}.fastq_map.txt`
- `gene.txt`
- `medianGeneExpression.txt.gz`
- `rat_ids/{tissue}.rat_ids.txt`
- `ref/RatGTEx_rats.tsv`
- `ref/RatGTEx_samples.tsv`
- `ref/rats.html`: This HTML table must be copied into `/about/samples/index.html`
- `ref/samples.html`: This HTML table must be copied into `/about/samples/index.html`
- `singleTissueEqtl.zip`
- `tissueInfo.txt`
- `topExpressedGene.txt`
- `splice/top_assoc_splice.txt`
- `splice/sqtls_indep.txt`
- `splice/{tissue}.leafcutter.bed.gz"`
- `splice/{tissue}.covar_splice.txt`
- `splice/{tissue}.splice.trans_qtl_pairs.txt.gz`

## Requirements

- pandas (Python)
- gtfparse (Python)
- tidyverse (R)

## Usage

Run `run.py`, which calls R and Python scripts that prepare the server data. The parameters are:

1. Path to the RatGTEx pipeline base directory
2. Output directory path
3. The remaining parameters are the list of tissues to include

Then, do these additional steps:

1. Copy genotype files (`{tissues}.vcf.gz` and `{tissues}.vcf.gz.tbi`) for each study dataset from the input into `{outdir}/geno/`.
2. Run `check.py` to check if all necessary files have been created/copied.
3. Copy the contents of `ref/rats.html` and `ref/samples.html` into `/about/samples/index.html` in the site directory, replacing the old tables.
4. `{outdir}/studies/{tissues}/` contain data directly used in the original publications, which may or may not differ from the data from the unified RatGTEx pipeline. These aren't used in the portal except to be available for download. Copy any files to this location and add links to them on the original study data section of the downloads page.

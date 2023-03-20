# ratgtex-server-data
Process RatGTEx pipeline results into data files for the RatGTEx site

Once the [RatGTEx pipeline](https://github.com/daniel-munro/ratgtex-pipeline) is run for all tissues, this code processes the results and assembles them into data files for [ratgtex.org](https://ratgtex.org). Most of them are put into a convenient format to be used by the [API](github.com/daniel-munro/ratgtex), which is used mainly for the data visualizations. This should be run separately for each genome version, and all output file names will contain the version. Below, `{rn}` is `rn6` or `rn7`.

## Inputs (in specified input directory):

- `tissue_info.txt`: Tab-delimited table of basic information about each tissue, e.g. full name and color to use in website visualizations.
- `geno_{rn}/{tissues}.vcf.gz{,.tbi}`: Genotypes for each study, each consisting of one or more tissues.
- `geno_{rn}/alleles.txt.gz`: Three-column tab-delimited file with variant ID, REF, and ALT alleles for the union of SNPs extracted from the genotype VCF files.
- `geno_{rn}/genotyping_log.csv`: CSV metadata file containing sex and coat color info to add to the rat info table.
- `ref_{rn}/GENES_RAT.txt`: Gene information downloaded from the Rat Genome Database. Used to show descriptive gene names.
- `samples/{accession}_series_matrix.txt.gz`: GEO series metadata for each tissue whose RNA-Seq data has been added to GEO. When a tissue is added to GEO, download this file for the series (not the superseries in case of multitissue datasets) and add the accession to `sample_table.R`.
- `ref_rn6/Rattus_norvegicus.Rnor_6.0.99.genes.gtf` or `ref_rn7/Rattus_norvegicus.mRatBN7.2.108.genes.gtf`
- `{rn}/{tissue}/covar.txt`
- `{rn}/{tissue}/fastq_map.txt`
- `{rn}/{tissue}/rat_ids.txt`
- `{rn}/{tissue}/{tissue}.aFC.txt`
- `{rn}/{tissue}/{tissue}.cis_independent_qtl.txt.gz`
- `{rn}/{tissue}/{tissue}.cis_qtl_all_pvals.txt.gz`
- `{rn}/{tissue}/{tissue}.cis_qtl_signif.txt.gz`
- `{rn}/{tissue}/{tissue}.cis_qtl.txt.gz`
- `{rn}/{tissue}/{tissue}.expr.iqn.bed.gz`
- `{rn}/{tissue}/{tissue}.expr.log2.bed.gz`
- `{rn}/{tissue}/{tissue}.expr.tpm.bed.gz`
- `{rn}/{tissue}/{tissue}.trans_qtl_pairs.txt.gz`
- `{rn}/{tissue}/splice/{tissue}.leafcutter.bed.gz`
- `{rn}/{tissue}/splice/{tissue}.covar_splice.txt`
- `{rn}/{tissue}/splice/{tissue}_splice.cis_qtl.txt.gz`
- `{rn}/{tissue}/splice/{tissue}_splice.cis_independent_qtl.txt.gz`
- `{rn}/{tissue}/splice/{tissue}_splice.trans_qtl_pairs.txt.gz`

## Outputs (in specified output directory):

- `geno/{tissues}.{rn}.vcf.gz{,.tbi}`
- `{rn}.autocomplete.json`
- `{rn}.exon.txt`
- `{rn}.gene.txt`
- `{rn}.medianGeneExpression.txt.gz`
- `{rn}.singleTissueEqtl.zip`
- `{rn}.tissueInfo.txt`
- `{rn}.topExpressedGene.txt`
- `cis_pvals/{tissue}.{rn}.zip`
- `covar/{tissue}.{rn}.covar.txt`
- `eqtl/{rn}.eqtls_indep.txt`
- `eqtl/{rn}.top_assoc.txt`
- `eqtl/{tissue}.{rn}.cis_qtl_signif.txt.gz`
- `eqtl/{tissue}.{rn}.trans_qtl_pairs.txt.gz`
- `expr/{tissue}.{rn}.expr.{iqn,log2,tpm}.bed.gz`
- `fastq_map/{tissue}.{rn}.fastq_map.txt`
- `rat_ids/{tissue}.{rn}.rat_ids.txt`
- `ref/{rn}.RatGTEx_rats.tsv`
- `ref/{rn}.RatGTEx_samples.tsv`
- `ref/{rn}.rats.html`: This HTML table must be copied into `/about/samples/index.html`
- `ref/{rn}.samples.html`: This HTML table must be copied into `/about/samples/index.html`
- `splice/{rn}.top_assoc_splice.txt`
- `splice/{rn}.sqtls_indep.txt`
- `splice/{tissue}.{rn}.leafcutter.bed.gz"`
- `splice/{tissue}.{rn}.covar_splice.txt`
- `splice/{tissue}.{rn}.splice.trans_qtl_pairs.txt.gz`

## Requirements

- pandas (Python)
- gtfparse (Python)
- tidyverse (R)

## Usage

Run `check_inputs.py` to check if all necessary input files exist. If so, run `run.py`, which calls R and Python scripts that prepare the server data. The parameters are:

1. Path to the RatGTEx pipeline base directory
2. Genome version, rn6 or rn7
3. Output directory path

Then, do these additional steps:

1. Run `check_outputs.py` to check if all necessary files have been created/copied.
2. Copy the contents of `ref/rats.html` and `ref/samples.html` into `/about/samples/index.html` in the site directory, replacing the old tables.
3. `{outdir}/studies/{tissues}/` contain data directly used in the original publications, which may or may not differ from the data from the unified RatGTEx pipeline. These aren't used in the portal except to be available for download. Copy any files to this location and add links to them on the original study data section of the downloads page.

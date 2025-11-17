# ratgtex-server-data

Process RatGTEx pipeline results into data files for the RatGTEx site

Once the [RatGTEx pipeline](https://github.com/daniel-munro/ratgtex-pipeline) is run for all tissues, this code processes the results and assembles them into data files for [ratgtex.org](https://ratgtex.org).
Many of them are put into a convenient format to be used by the [API](github.com/daniel-munro/ratgtex), which is used mainly for the data visualizations.
Others are provided for download.
This code is designed to run only the latest RatGTEx data version, but files are still labeled with RatGTEx version, and sometimes genome version, for clarity.

## Inputs (in specified input directory):

- `tissue_info.tsv`: Tab-delimited table of basic information about each tissue, e.g. full name and color to use in website visualizations.
- `geno/{geno_dataset}.vcf.gz{,.tbi}`: Genotype files, each including the individuals for one or more tissues.
- `geno/alleles.txt.gz`: Three-column tab-delimited file with variant ID, REF, and ALT alleles for the union of SNPs extracted from the genotype VCF files.
- `geno/genotyping_log.csv`: CSV metadata file containing sex and coat color info to add to the rat info table.
- `samples/{accession}_series_matrix.txt.gz`: GEO series metadata for each tissue whose RNA-Seq data has been added to GEO. When a tissue is added to GEO, download this file for the series (not the superseries in case of multitissue datasets) and add the accession to `sample_table.R`.
- `ref/GCF_036323735.1_GRCr8_genomic.chr.gtf`
- `ref/GCF_036323735.1_GRCr8_genomic.chr.genes.gtf`
- `v4/{tissue}/fastq_map.txt`
- `v4/{tissue}/rat_ids.txt`
- `v4/{tissue}/{tissue}.aFC.txt`
- `v4/{tissue}/{tissue}.cis_qtl_all_pvals.txt.gz`
- `v4/{tissue}/{tissue}.{modality}.cis_qtl_signif.txt.gz`
- `v4/{tissue}/{tissue}.expr.log2.bed.gz`
- `v4/{tissue}/{tissue}.expression.trans_qtl_pairs.txt.gz`
- `v4/{tissue}/phenos/output/{modality}.bed.gz`
- `v4/{tissue}/phenos/output/unnorm/{modality}.bed`
- `v4/{tissue}/pheast/intermediate/covar/{modality}.covar.tsv`
- `v4/{tissue}/pheast/output/qtl/{modality}.cis_independent_qtl.txt.gz`
- `v4/{tissue}/pheast/output/qtl/{modality}.cis_qtl.txt.gz`

## Outputs (in specified output directory):

- `geno/{geno_dataset}.rn8.vcf.gz{,.tbi}`
- `autocomplete.v4.json`
- `exon.v4.tsv`
- `gene.v4.tsv`
- `medianGeneExpression.v4.tsv.gz`
- `singleTissueEqtl.v4.zip`
- `tissue_info.v4.tsv`
- `topExpressedGene.v4.tsv`
- `cis_pvals/{tissue}.v4.zip`
- `covar/covar.{tissue}.{modality}.v4.tsv`
- `eqtl/eqtls_indep.v4_rn8.tsv`
- `eqtl/top_assoc.v4_rn8.tsv`
- `expr/expr.{log2,tpm}.{tissue}.v4_rn8.bed.gz`
- `fastq_map/fastq_map.{tissue}.v4.txt`
- `phenos/phenos.{tissue}.{modality}.norm.v4_rn8.bed.gz`
- `phenos/phenos.{tissue}.{modality}.unnorm.v4_rn8.bed.gz`
- `rat_ids/rat_ids.{tissue}.v4.txt`
- `ref/RatGTEx_rats.v4.tsv`
- `ref/RatGTEx_samples.v4.tsv`
- `ref/rats.v4.html`: Jekyll will insert this HTML table into the About Samples page.
- `ref/samples.v4.html`: Jekyll will insert this HTML table into the About Samples page.
- `xqtl/cis_qtl_signif.{tissue}.{modality}.v4_rn8.txt.gz`
- `xqtl/top_assoc.{modality}.v4_rn8.tsv`
- `xqtl/trans_qtl_pairs.{tissue}.expression.v4_rn8.txt.gz`
- `xqtl/xqtls_indep.{modality}.v4_rn8.tsv`

## Requirements

- snakemake
- pandas (Python)
- gtfparse (Python)
- tidyverse (R)

## Usage

Run [Snakemake](https://snakemake.github.io) to generate the output files. Edit the Snakefile as needed to change which files are produced. It is recommended to first run as a dry run (with `-n`) to ensure it will run the steps you expect it to.

`{outdir}/studies/{tissues}/` contain data directly used in the original publications, which may or may not differ from the data from the unified RatGTEx pipeline. These aren't used in the portal except to be available for download. Copy any files to this location and add links to them on the original study data page.

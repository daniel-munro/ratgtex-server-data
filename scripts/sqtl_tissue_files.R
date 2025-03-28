# Process and aggregate sQTL results from tensorQTL
#
# Inputs:
#   {indir}/ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf
#   {indir}/v3/{tissue}/splice/{tissue}_splice.cis_qtl_signif.txt.gz
#   {indir}/v3/{tissue}/splice/{tissue}_splice.trans_qtl_pairs.txt.gz
#
# Outputs:
#   {outdir}/splice/splice.cis_qtl_signif.{tissue}.v3_rn7.txt.gz
#   {outdir}/splice/splice.trans_qtl_pairs.{tissue}.v3_rn7.txt.gz

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissue <- args[3]

version <- "v3"
v <- "v3_rn7"
anno <- "GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf"

genes <- read_tsv(str_glue("{indir}/ref/{anno}"),
                  col_types = "--cii-c-c",
                  col_names = c("type", "start", "end", "strand", "etc"),
                  comment = "#") |>
    filter(type == "gene") |>
    mutate(gene_id = str_match(etc, 'gene_id "([^"]+)"')[, 2],
           tss = if_else(strand == "+", start, end)) |>
    select(gene_id, strand, tss)
stopifnot(sum(duplicated(genes$gene_id)) == 0)

#################################
## Copy and modify signif file ##
#################################

infile <- str_glue("{indir}/{version}/{tissue}/splice/{tissue}_splice.cis_qtl_signif.txt.gz")
outfile <- str_glue("{outdir}/splice/splice.cis_qtl_signif.{tissue}.{v}.txt.gz")
read_tsv(infile, col_types = "cccccccccc") |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
                remove = FALSE) |>
    mutate(gene_id = str_match(phenotype_id, ":([^:]+)$")[, 2]) |>
    left_join(select(genes, gene_id, strand, tss), by = "gene_id", relationship = "many-to-one") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
    select(-gene_id, -chrom, -pos, -strand, -tss) |>
    write_tsv(outfile)

################################
## Copy and modify trans file ##
################################

infile <- str_glue("{indir}/{version}/{tissue}/splice/{tissue}_splice.trans_qtl_pairs.txt.gz")
outfile <- str_glue("{outdir}/splice/splice.trans_qtl_pairs.{tissue}.{v}.txt.gz")
read_tsv(infile, col_types = "cccccc") |>
    rename(slope = b,
            slope_se = b_se) |>
    write_tsv(outfile)

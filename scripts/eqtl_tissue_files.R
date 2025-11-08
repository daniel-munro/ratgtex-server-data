# Process eQTL results for one tissue

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissue <- args[3]

input_anno <- str_glue("{indir}/ref/GCF_036323735.1_GRCr8_genomic.chr.genes.gtf")
input_signif <- str_glue("{indir}/v4/{tissue}/{tissue}.cis_qtl_signif.txt.gz")
input_trans <- str_glue("{indir}/v4/{tissue}/{tissue}.trans_qtl_pairs.txt.gz")

output_signif <- str_glue("{outdir}/eqtl/cis_qtl_signif.{tissue}.v4_rn8.txt.gz")
output_trans <- str_glue("{outdir}/eqtl/trans_qtl_pairs.{tissue}.v4_rn8.txt.gz")

genes <- read_tsv(
  input_anno,
  col_types = "--cii-c-c",
  col_names = c("type", "start", "end", "strand", "etc"),
  comment = "#"
) |>
  filter(type == "gene") |>
  mutate(
    gene_id = str_match(etc, 'gene_id "([^"]+)"')[, 2],
    tss = if_else(strand == "+", start, end)
  ) |>
  select(gene_id, strand, tss)
stopifnot(sum(duplicated(genes$gene_id)) == 0)

#################################
## Copy and modify signif file ##
#################################

read_tsv(input_signif, col_types = "cccccccccc") |>
  rename(gene_id = phenotype_id) |>
  separate(
    variant_id, c("chrom", "pos"), sep = ":", convert = TRUE, remove = FALSE
  ) |>
  left_join(
    select(genes, gene_id, strand, tss),
    by = "gene_id",
    relationship = "many-to-one"
  ) |>
  mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
  select(-chrom, -pos, -strand, -tss) |>
  write_tsv(output_signif)

################################
## Copy and modify trans file ##
################################

read_tsv(input_trans, col_types = "cccccc") |>
  rename(
    gene_id = phenotype_id,
    slope = b,
    slope_se = b_se
  ) |>
  write_tsv(output_trans)

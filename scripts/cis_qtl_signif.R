# Process all significant cis-xQTL pair file for one tissue

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
input_signif <- args[1]
input_genes <- args[2]
output_signif <- args[3]

genes <- read_tsv(
  input_genes,
  col_types = cols(geneId = "c", strand = "c", tss = "i", .default = "-")
)

read_tsv(input_signif, col_types = "cccccccccc") |>
  separate(
    variant_id, c("chrom", "pos"), sep = ":", convert = TRUE, remove = FALSE
  ) |>
  mutate(geneId = str_split_i(phenotype_id, "__", 1)) |>
  left_join(genes, by = "geneId", relationship = "many-to-one") |>
  mutate(start_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
  rename(tss_distance = start_distance) |>
  select(-geneId, -chrom, -pos, -strand, -tss) |>
  write_tsv(output_signif)

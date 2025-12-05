# Generate tissue info table for web interface

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissues <- args[3:length(args)]

version <- "v4"

input_tissue <- str_glue("{indir}/tissue_info.tsv")
input_expr_template <- "{indir}/{version}/{tissue}/phenos/output/expression.bed.gz"
input_assoc <- str_glue("{outdir}/eqtl/top_assoc.{version}_rn8.tsv")

output <- str_glue("{outdir}/tissue_info.{version}.tsv")

samples <- tibble(tissue = tissues) |>
  reframe(
    {
      df <- read_tsv(
        str_glue(input_expr_template),
        col_types = cols(.default = "c"),
        n_max = 0
      )
      tibble(rnaSeqSampleCount = ncol(df) - 4)
    },
    .by = tissue
  )

expr <- tibble(tissue = tissues) |>
  reframe(
    read_tsv(
      str_glue(input_expr_template),
      col_types = cols(phenotype_id = "c", .default = "-")
    ) |>
      count(name = "expressedGeneCount"),
    .by = tissue
  )

egenes <- read_tsv(input_assoc, col_types = "ccicciccdiddddddd") |>
  filter(qval < 0.05) |>
  count(tissue, name = "eGeneCount")

d <- read_tsv(input_tissue, col_types = "cccccc")
stopifnot(all(tissues %in% d$tissue))
d <- d |>
  filter(tissue %in% tissues) |>
  mutate(tissueSiteDetailAbbr = tissue, .after = tissue) |>
  left_join(samples, by = "tissue", relationship = "one-to-one") |>
  left_join(expr, by = "tissue", relationship = "one-to-one") |>
  mutate(
    rnaSeqAndGenotypeSampleCount = rnaSeqSampleCount,
    .after = rnaSeqSampleCount
  ) |>
  left_join(egenes, by = "tissue", relationship = "one-to-one") |>
  mutate(
    datasetId = str_glue("ratgtex_{version}"),
    hasEGenes = "True"
  ) |>
  rename(tissueSiteDetailId = tissue)

write_tsv(d, output)

# Calculate median gene expression for web interface

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissues <- args[3:length(args)]

version <- "v4"

input_gene <- str_glue("{outdir}/gene.{version}.txt")
input_expr_template <- "{outdir}/expr/expr.tpm.{tissue}.{version}_rn8.bed.gz"

output_median <- str_glue("{outdir}/medianGeneExpression.{version}.txt.gz")
output_top <- str_glue("{outdir}/topExpressedGene.{version}.txt")
chroms <- read_tsv(
  input_gene,
  col_types = cols(geneId = "c", chromosome = "c", .default = "-")
)

expr <- tibble(tissue = tissues) |>
  reframe(
    read_tsv(str_glue(input_expr_template),
      col_types = c(
        `#chr` = "-", start = "-", end = "-",
        phenotype_id = "c", .default = "d"
      )
    ) |>
      pivot_longer(-phenotype_id, names_to = "rat_id", values_to = "tpm"),
    .by = tissue
  ) |>
  rename(
    tissueSiteDetailId = tissue,
    geneId = phenotype_id
  )

med <- expr |>
  group_by(tissueSiteDetailId, geneId) |>
  summarise(median = round(median(tpm), 3), .groups = "drop")

med |>
  pivot_wider(names_from = tissueSiteDetailId, values_from = median) |>
  write_tsv(output_median)

top <- med |>
  filter(geneId %in% chroms$geneId) |>
  left_join(chroms, by = "geneId", relationship = "many-to-one") |>
  mutate(mtGene = if_else(chromosome == "chrM", "True", "False")) |>
  arrange(tissueSiteDetailId, desc(median)) |>
  # slice(1:50) |>
  filter(
    cumsum(mtGene == "False") <= 50, # Need 50 non-MT genes for heatmap
    .by = tissueSiteDetailId
  ) |>
  mutate(
    datasetId = str_glue("ratgtex_{version}"),
    unit = "TPM"
  )

write_tsv(top, output_top)

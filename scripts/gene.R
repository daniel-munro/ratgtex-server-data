# Generate gene info table for web interface

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissues <- args[3:length(args)]

input_anno <- str_glue("{indir}/ref/GCF_036323735.1_GRCr8_genomic.chr.gtf")
input_signif_template <- "{indir}/v4/{tissue}/{tissue}.cis_qtl_signif.txt.gz"
input_expr_template <- "{indir}/v4/{tissue}/phenos/output/expression.bed.gz"
input_eqtl_indep <- str_glue("{outdir}/eqtl/eqtls_indep.v4_rn8.tsv")
input_top_assoc <- str_glue("{outdir}/eqtl/top_assoc.v4_rn8.tsv")

output_gene <- str_glue("{outdir}/gene.v4.tsv")
output_autocomplete <- str_glue("{outdir}/autocomplete.v4.json")

genes <- read_tsv(
  input_anno,
  col_types = "c-cii-c-c",
  col_names = c("chromosome", "type", "start", "end", "strand", "etc"),
  comment = "#"
) |>
  filter(
    type == "gene",
    chromosome %in% str_glue("chr{c(1:20, 'X', 'Y', 'M')}")
  ) |>
  mutate(
    geneId = str_match(etc, 'gene_id "([^"]+)"')[, 2],
    description = str_match(etc, 'description "([^"]+)"')[, 2],
    tss = if_else(strand == "+", start, end)
  ) |>
  replace_na(list(description = "")) |>
  select(geneId, chromosome, start, end, strand, tss, description)

# Not sure if this is exact same set as in top_assoc.txt, but this is to see if
# file of signif pairs exists for each gene:
signif <- tibble(tissue = tissues) |>
  reframe(
    read_tsv(str_glue(input_signif_template),
      col_types = "c---------"
    ),
    .by = tissue
  ) |>
  distinct(phenotype_id) |>
  pull()

## Add expression/eQTL status for gene page

expr <- tibble(tissue = tissues) |>
  reframe(
    read_tsv(str_glue(input_expr_template),
      col_types = c(phenotype_id = "c", .default = "-")
    ),
    .by = tissue
  ) |>
  rename(gene_id = phenotype_id)

expr_all <- expr |>
  distinct(gene_id) |>
  pull()

is_expr <- expr |>
  mutate(expressed = "True") |>
  complete(tissue, gene_id = genes$geneId, fill = list(expressed = "False")) |>
  pivot_wider(
    id_cols = gene_id, names_from = tissue, names_prefix = "expr_",
    values_from = expressed, values_fill = "False"
  )

was_tested_eqtl <- read_tsv(input_top_assoc, col_types = "cc----------------") |>
  mutate(testedEqtl = "True") |>
  complete(tissue, gene_id = genes$geneId, fill = list(testedEqtl = "False")) |>
  pivot_wider(
    id_cols = gene_id, names_from = tissue, names_prefix = "testedEqtl_",
    values_from = testedEqtl
  )

has_eqtl <- read_tsv(input_eqtl_indep, col_types = "cc---------------") |>
  filter(tissue %in% tissues) |>
  distinct(tissue, gene_id) |>
  mutate(eqtl = "True") |>
  complete(tissue, gene_id = genes$geneId, fill = list(eqtl = "False")) |>
  pivot_wider(
    id_cols = gene_id, names_from = tissue, names_prefix = "eqtl_",
    values_from = eqtl
  )

## Assemble

df <- genes |>
  mutate(hasEqtl = if_else(geneId %in% signif, "True", "False")) |>
  left_join(is_expr, by = c("geneId" = "gene_id"), relationship = "one-to-one") |>
  left_join(was_tested_eqtl, by = c("geneId" = "gene_id"), relationship = "one-to-one") |>
  left_join(has_eqtl, by = c("geneId" = "gene_id"), relationship = "one-to-one")

write_tsv(df, output_gene)

# Save list of all gene IDs for search autocomplete:
df |>
  pull(geneId) |>
  sort() |>
  jsonlite::toJSON() |>
  write_file(output_autocomplete)

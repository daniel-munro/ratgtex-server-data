# Calculate median gene expression for web interface
#
# Inputs:
#   {indir}/v3/{tissue}/{tissue}.expr.tpm.bed.gz
#   {outdir}/gene.v3.txt
#
# Outputs:
#   {outdir}/medianGeneExpression.v3.txt.gz
#   {outdir}/topExpressedGene.v3.txt

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissues <- args[3:length(args)]

version <- "v3"

chroms <- read_tsv(str_glue("{outdir}/gene.{version}.txt"),
    col_types = cols(geneId = "c", chromosome = "c", .default = "-")
)

expr <- tibble(tissue = tissues) |>
    reframe(
        read_tsv(str_glue("{indir}/{version}/{tissue}/{tissue}.expr.tpm.bed.gz"),
                 col_types = c(`#chr` = "-", start = "-", end = "-",
                               gene_id = "c", .default = "d")) |>
            pivot_longer(-gene_id, names_to = "rat_id", values_to = "tpm"),
        .by = tissue
    ) |>
    rename(tissueSiteDetailId = tissue,
           geneId = gene_id)

med <- expr |>
    group_by(tissueSiteDetailId, geneId) |>
    summarise(median = round(median(tpm), 3), .groups = "drop")

med |>
    pivot_wider(names_from = tissueSiteDetailId, values_from = median) |>
    write_tsv(str_glue("{outdir}/medianGeneExpression.{version}.txt.gz"))

top <- med |>
    filter(geneId %in% chroms$geneId) |>
    left_join(chroms, by = "geneId", relationship = "many-to-one") |>
    mutate(mtGene = if_else(chromosome == "chrM", "True", "False")) |>
    arrange(tissueSiteDetailId, desc(median)) |>
    group_by(tissueSiteDetailId) |>
    # slice(1:50) |>
    filter(cumsum(mtGene == "False") <= 50) |> # Need 50 non-MT genes for heatmap
    ungroup() |>
    mutate(datasetId = str_glue("ratgtex_{version}"),
           unit = "TPM")

write_tsv(top, str_glue("{outdir}/topExpressedGene.{version}.txt"))

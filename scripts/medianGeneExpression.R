# Calculate median gene expression for web interface
#
# Inputs:
#   {indir}/{version}/{tissue}/{tissue}.expr.tpm.bed.gz
#   {outdir}/gene.{v}.txt
#
# Outputs:
#   {outdir}/medianGeneExpression.{v}.txt.gz
#   {outdir}/topExpressedGene.{v}.txt

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissues <- args[3:length(args)]

version <- "v3"
v <- "v3_rn7"

genes <- read_tsv(str_glue("{outdir}/gene.{v}.txt"), col_types = "cc-------")

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
    summarise(median = median(tpm), .groups = "drop")

med |>
    pivot_wider(names_from = tissueSiteDetailId, values_from = median) |>
    write_tsv(str_glue("{outdir}/medianGeneExpression.{v}.txt.gz"))

top <- med |>
    filter(geneId %in% genes$geneId) |>
    left_join(genes, by = "geneId", relationship = "many-to-one") |>
    mutate(mtGene = if_else(str_sub(geneSymbol, 1, 3) == "Mt-", "True", "False")) |>
    arrange(tissueSiteDetailId, desc(median)) |>
    group_by(tissueSiteDetailId) |>
    # slice(1:50) |>
    filter(cumsum(mtGene == "False") <= 50) |> # Need 50 non-MT genes for heatmap
    ungroup() |>
    mutate(datasetId = "ratgtex_v1",
           unit = "TPM")

write_tsv(top, str_glue("{outdir}/topExpressedGene.{v}.txt"))

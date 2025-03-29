# Generate tissue info table for web interface
#
# Inputs:
#   {indir}/tissue_info.txt
#   {indir}/v3/{tissueSiteDetailId}/{tissueSiteDetailId}.expr.tpm.bed.gz
#   {outdir}/eqtl/top_assoc.v3_rn7.txt
#
# Outputs:
#   {outdir}/tissueInfo.v3_rn7.txt

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissues <- args[3:length(args)]

version <- "v3"
v <- "v3_rn7"

expr <- tibble(tissueSiteDetailId = tissues) |>
    reframe(
        read_tsv(
            str_glue("{indir}/{version}/{tissueSiteDetailId}/{tissueSiteDetailId}.expr.tpm.bed.gz"),
            col_types = c(
                `#chr` = "-", start = "-", end = "-",
                gene_id = "c", .default = "d"
            )
        ) |>
            pivot_longer(-gene_id, names_to = "rat_id", values_to = "tpm") |>
            group_by(gene_id) |>
            filter(sum(tpm > 0) > 1) |> # Genes tested for eQTLs
            ungroup() |>
            summarise(
                rnaSeqSampleCount = n_distinct(rat_id),
                expressedGeneCount = n_distinct(gene_id)
            ),
        .by = tissueSiteDetailId
    )

egenes <- read_tsv(str_glue("{outdir}/eqtl/top_assoc.{v}.txt"), col_types = "cccicciccdiddddddd") |>
    filter(qval < 0.05) |>
    count(tissue, name = "eGeneCount")

d <- read_tsv(str_glue("{indir}/tissue_info.txt"), col_types = "ccccc")
stopifnot(all(tissues %in% d$tissueSiteDetailId))
d <- d |>
    filter(tissueSiteDetailId %in% tissues) |>
    mutate(tissueSiteDetailAbbr = tissueSiteDetailId, .after = tissueSiteDetailId) |>
    left_join(expr, by = "tissueSiteDetailId", relationship = "one-to-one") |>
    mutate(rnaSeqAndGenotypeSampleCount = rnaSeqSampleCount, .after = rnaSeqSampleCount) |>
    left_join(egenes, by = c("tissueSiteDetailId" = "tissue"), relationship = "one-to-one") |>
    mutate(
        datasetId = str_glue("ratgtex_{version}"),
        hasEGenes = "True"
    )

write_tsv(d, str_glue("{outdir}/tissueInfo.{version}.txt"))

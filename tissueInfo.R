suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissues <- args[3:length(args)]

expr <- tibble(tissueSiteDetailId = tissues) |>
    group_by(tissueSiteDetailId) |>
    summarise(
        read_tsv(str_glue("{indir}/{tissueSiteDetailId}/{tissueSiteDetailId}.expr.tpm.bed.gz"),
                 col_types = c(`#chr` = "-", start = "-", end = "-",
                               gene_id = "c", .default = "d")) |>
            pivot_longer(-gene_id, names_to = "rat_id", values_to = "tpm") |>
            group_by(gene_id) |>
            filter(sum(tpm > 0) > 1) |> # Genes tested for eQTLs
            ungroup() |>
            summarise(rnaSeqSampleCount = n_distinct(rat_id),
                      expressedGeneCount = n_distinct(gene_id)),
        .groups = "drop"
    )

egenes <- read_tsv(str_glue("{outdir}/eqtl/top_assoc.txt"), col_types = "cccicciccdiddddddd") |>
    filter(qval < 0.05) |>
    count(tissue, name = "eGeneCount")

d <- read_tsv(str_glue("{indir}/tissue_info.txt"), col_types = "ccccc")
stopifnot(all(tissues %in% d$tissueSiteDetailId))
d <- d |>
    filter(tissueSiteDetailId %in% tissues) |>
    mutate(tissueSiteDetailAbbr = tissueSiteDetailId, .after = tissueSiteDetailId) |>
    left_join(expr, by = "tissueSiteDetailId") |>
    mutate(rnaSeqAndGenotypeSampleCount = rnaSeqSampleCount, .after = rnaSeqSampleCount) |>
    left_join(egenes, by = c("tissueSiteDetailId" = "tissue")) |>
    mutate(datasetId = "ratgtex_v1",
           hasEGenes = "True")

write_tsv(d, str_glue("{outdir}/tissueInfo.txt"))

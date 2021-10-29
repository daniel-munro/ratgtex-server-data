library(tidyverse)

genes <- read_tsv("server_data/gene.txt", col_types = "cc-------")

# # tmp: make TPM BED files for brain regions
# genes <- read_tsv("server_data/gene.txt", col_types = "c-c---i--")
# expr_br <- read.delim(
#     "~/Dropbox (Scripps Research)/HS-RNASeq/quantitation/EnsemblGene_v2/ensembl-gene_raw-tpm.txt",
#     check.names = FALSE
# ) |>
#     as_tibble(rownames = "geneId") |>
#     pivot_longer(-geneId, names_to = "library", values_to = "tpm") |>
#     # separate(sample, c("rat_id", "tissue"), sep = "_")
#     # mutate(tissue = str_sub(sample, 12, -1)) |>  # for speed
#     inner_join(samples, by = "library")
# for (tissue in unique(expr_br$tissueSiteDetailId)) {
#     x <- expr_br |>
#         filter(tissueSiteDetailId == tissue) |>
#         select(-tissueSiteDetailId) |>
#         mutate(library = str_split(library, "_", simplify = TRUE)[, 1]) |>
#         pivot_wider(id_cols = geneId, names_from = library, values_from = tpm)
#     y <- genes |>
#         mutate(start = tss - 1) |>
#         select(`#chr` = chromosome, start = tss, end = tss, geneId) |>
#         right_join(x, by = "geneId")
#     # print(y)
#     write_tsv(y, str_glue("data/expr/{tissue}.expr.tpm.bed"))
# }

tissues <- c("Eye", "IL", "LHb", "NAcc", "OFC", "PL")

expr <- tibble(tissue = tissues) |>
    group_by(tissue) |>
    summarise(
        read_tsv(str_glue("data/expr/{tissue}.expr.tpm.bed.gz"),
                 col_types = c(`#chr` = "-", start = "-", end = "-",
                               gene_id = "c", .default = "d")) |>
            pivot_longer(-gene_id, names_to = "rat_id", values_to = "tpm"),
        .groups = "drop"
    ) |>
    rename(tissueSiteDetailId = tissue,
           geneId = gene_id)

med <- expr |>
    group_by(tissueSiteDetailId, geneId) |>
    summarise(median = median(tpm), .groups = "drop")

med |>
    pivot_wider(names_from = tissueSiteDetailId, values_from = median) |>
    write_tsv("server_data/medianGeneExpression.txt.gz")

top <- med |>
    filter(geneId %in% genes$geneId) |>
    left_join(genes, by = "geneId") |>
    mutate(mtGene = if_else(str_sub(geneSymbol, 1, 3) == "Mt-", "True", "False")) |>
    arrange(tissueSiteDetailId, desc(median)) |>
    group_by(tissueSiteDetailId) |>
    # slice(1:50) |>
    filter(cumsum(mtGene == "False") <= 50) |>
    ungroup() |>
    mutate(datasetId = "ratgtex_v1",
           unit = "TPM")

write_tsv(top, "server_data/topExpressedGene.txt")

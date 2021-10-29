library(tidyverse)

descs <- read_tsv(
    "data/ref/GENES_RAT.txt",
    col_types = cols(ENSEMBL_ID = "c", NAME = "c", .default = "-"),
    skip = 83
) |>
    rename(geneId = ENSEMBL_ID,
           description = NAME) |>
    filter(!is.na(geneId)) |>
    separate_rows(geneId, sep=";") |>
    group_by(geneId) |>
    slice(1) |>
    ungroup()

tissues <- c("Eye", "IL", "LHb", "NAcc", "OFC", "PL")

# Not sure if this is exact same set as in top_assoc.txt, but this is to see if
# file of signif pairs exists for each gene:
signif <- tibble(tissue = tissues) |>
    group_by(tissue) |>
    summarise(
        read_tsv(str_glue("data/tensorqtl/{tissue}.cis_qtl_signif.txt.gz"),
                 col_types = "c---------"),
        .groups = "drop"
    ) |>
    distinct(gene_id) |>
    pull()

expr <- tibble(tissue = tissues) |>
    group_by(tissue) |>
    summarise(
        read_tsv(str_glue("data/expr/{tissue}.expr.tpm.bed.gz"),
                 col_types = c(`#chr` = "-", start = "-", end = "-",
                               gene_id = "c", .default = "d")) |>
            pivot_longer(-gene_id, names_to = "rat_id", values_to = "tpm") |>
            group_by(gene_id) |>
            summarise(expressed = sum(tpm > 0) > 1, # Genes tested for eQTLs
                      .groups = "drop"),
        .groups = "drop"
    ) |>
    filter(expressed) |>
    distinct(gene_id) |>
    pull()

genes <- read_tsv("data/ref/Rattus_norvegicus.Rnor_6.0.99.genes.bed",
                  col_types = "ciic-c---c",
                  col_names = c("chromosome", "start", "end", "geneId", "strand", "etc")) %>%
    mutate(geneSymbol = str_match(etc, 'gene_name "([^"]+)"')[, 2],
           tss = if_else(strand == "+", start, end)) %>%
    select(geneId, geneSymbol, chromosome, start, end, strand, tss) |>
    left_join(descs, by = "geneId") |>
    mutate(hasEqtl = if_else(geneId %in% signif, "True", "False")) |>
    replace_na(list(description = ""))

write_tsv(genes, "server_data/gene.txt")

# Save list of all (expressed) names and IDs for search autocomplete:
genes2 <- genes |>
    filter(geneId %in% expr)
c(genes2$geneId, genes2$geneSymbol) |>
    unique() |>
    sort() |>
    jsonlite::toJSON() |>
    write_file("server_data/autocomplete.json")

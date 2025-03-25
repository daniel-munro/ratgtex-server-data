suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissues <- args[3:length(args)]

descs <- read_tsv(
    str_glue("{indir}/ref_{rn}/GENES_RAT.txt"),
    col_types = cols(ENSEMBL_ID = "c", NAME = "c", .default = "-"),
    skip = c(rn6 = 83, rn7 = 85)[rn] # Some lines have '#' in the text, so I can't use '#' as comment char to skip pre-header
) |>
    rename(geneId = ENSEMBL_ID,
           description = NAME) |>
    filter(!is.na(geneId)) |>
    separate_rows(geneId, sep=";") |>
    group_by(geneId) |>
    slice(1) |>
    ungroup()

# Not sure if this is exact same set as in top_assoc.txt, but this is to see if
# file of signif pairs exists for each gene:
signif <- tibble(tissue = tissues) |>
    reframe(
        read_tsv(str_glue("{indir}/{rn}/{tissue}/{tissue}.cis_qtl_signif.txt.gz"),
                 col_types = "c---------"),
        .by = tissue
    ) |>
    distinct(phenotype_id) |>
    pull()

## Add expression/eQTL status for gene page

expr <- tibble(tissue = tissues) |>
    reframe(
        read_tsv(str_glue("{indir}/{rn}/{tissue}/{tissue}.expr.tpm.bed.gz"),
                 col_types = c(`#chr` = "-", start = "-", end = "-",
                               gene_id = "c", .default = "d")) |>
            pivot_longer(-gene_id, names_to = "rat_id", values_to = "tpm") |>
            group_by(gene_id) |>
            summarise(expressed = sum(tpm > 0) > 1,
                      .groups = "drop"),
        .by = tissue
    )
expr_all <- expr |>
    filter(expressed) |>
    distinct(gene_id) |>
    pull()
is_expr <- expr |>
    mutate(expressed = if_else(expressed, "True", "False")) |>
    pivot_wider(id_cols = gene_id, names_from = tissue, names_prefix = "expr_",
                values_from = expressed, values_fill = "False")

was_tested_eqtl <- read_tsv(str_glue("{outdir}/eqtl/{rn}.top_assoc.txt"),
                   col_types = "cc----------------") |>
    mutate(testedEqtl = "True") |>
    complete(tissue, gene_id = expr$gene_id, fill = list(testedEqtl = "False")) |>
    pivot_wider(id_cols = gene_id, names_from = tissue, names_prefix = "testedEqtl_",
                values_from = testedEqtl)

has_eqtl <- read_tsv(str_glue("{outdir}/eqtl/{rn}.eqtls_indep.txt"),
                     col_types = "cc---------------") |>
    filter(tissue %in% tissues) |>
    distinct(tissue, gene_id) |>
    mutate(eqtl = "True") |>
    complete(tissue, gene_id = expr$gene_id, fill = list(eqtl = "False")) |>
    pivot_wider(id_cols = gene_id, names_from = tissue, names_prefix = "eqtl_",
                values_from = eqtl)

## Add alt splicing/sQTL status for gene page

splice <- tibble(tissue = tissues) |>
    reframe(
        read_tsv(str_glue("{indir}/{rn}/{tissue}/splice/{tissue}.leafcutter.phenotype_groups.txt"),
                 skip = 1, col_names = c("phenotype_id", "gene_id"), col_types = "-c"),
        .by = tissue
    ) |>
    distinct(tissue, gene_id)
is_alt_spliced <- splice |>
    mutate(altSpliced = "True") |>
    complete(tissue, gene_id = expr$gene_id, fill = list(altSpliced = "False")) |>
    pivot_wider(id_cols = gene_id, names_from = tissue, names_prefix = "altSplice_",
                values_from = altSpliced)

was_tested_sqtl <- read_tsv(str_glue("{outdir}/splice/{rn}.top_assoc_splice.txt"),
                   col_types = "c-------------c----") |>
    distinct(tissue, gene_id) |>
    mutate(testedSqtl = "True") |>
    complete(tissue, gene_id = expr$gene_id, fill = list(testedSqtl = "False")) |>
    pivot_wider(id_cols = gene_id, names_from = tissue, names_prefix = "testedSqtl_",
                values_from = testedSqtl)

has_sqtl <- read_tsv(str_glue("{outdir}/splice/{rn}.sqtls_indep.txt"),
                     col_types = "c-------------c---") |>
    filter(tissue %in% tissues) |>
    distinct(tissue, gene_id) |>
    mutate(sqtl = "True") |>
    complete(tissue, gene_id = expr$gene_id, fill = list(sqtl = "False")) |>
    pivot_wider(id_cols = gene_id, names_from = tissue, names_prefix = "sqtl_",
                values_from = sqtl)

## Assemble

anno <- c(
    rn6 = "Rattus_norvegicus.Rnor_6.0.99.genes.gtf",
    rn7 = "Rattus_norvegicus.mRatBN7.2.108.genes.gtf"
)[rn]

genes <- read_tsv(
    str_glue("{indir}/ref_{rn}/{anno}"),
    col_types = "c-cii-c-c",
    col_names = c("chromosome", "type", "start", "end", "strand", "etc"),
    comment = "#"
) |>
    filter(type == "gene") |>
    mutate(geneId = str_match(etc, 'gene_id "([^"]+)"')[, 2],
           geneSymbol = str_match(etc, 'gene_name "([^"]+)"')[, 2],
           tss = if_else(strand == "+", start, end)) |>
    select(geneId, geneSymbol, chromosome, start, end, strand, tss) |>
    left_join(descs, by = "geneId", relationship = "one-to-one") |>
    replace_na(list(description = "")) |>
    mutate(hasEqtl = if_else(geneId %in% signif, "True", "False")) |>
    left_join(is_expr, by = c("geneId" = "gene_id"), relationship = "one-to-one") |>
    left_join(was_tested_eqtl, by = c("geneId" = "gene_id"), relationship = "one-to-one") |>
    left_join(has_eqtl, by = c("geneId" = "gene_id"), relationship = "one-to-one") |>
    left_join(is_alt_spliced, by = c("geneId" = "gene_id"), relationship = "one-to-one") |>
    left_join(was_tested_sqtl, by = c("geneId" = "gene_id"), relationship = "one-to-one") |>
    left_join(has_sqtl, by = c("geneId" = "gene_id"), relationship = "one-to-one")

write_tsv(genes, str_glue("{outdir}/{rn}.gene.txt"))

# Save list of all (expressed) names and IDs for search autocomplete:
genes2 <- genes |>
    filter(geneId %in% expr_all)
c(genes2$geneId, genes2$geneSymbol) |>
    unique() |>
    sort() |>
    jsonlite::toJSON() |>
    write_file(str_glue("{outdir}/{rn}.autocomplete.json"))

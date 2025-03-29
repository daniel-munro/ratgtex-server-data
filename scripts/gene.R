# Generate gene info table for web interface
#
# Inputs:
#   {indir}/ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf
#   {indir}/ref/GENES_RAT.txt
#   {indir}/{version}/{tissue}/{tissue}.cis_qtl_signif.txt.gz
#   {indir}/{version}/{tissue}/{tissue}.expr.tpm.bed.gz
#   {indir}/{version}/{tissue}/splice/{tissue}.leafcutter.phenotype_groups.txt
#   {outdir}/eqtl/eqtls_indep.v3_rn7.txt
#   {outdir}/eqtl/top_assoc.v3_rn7.txt
#   {outdir}/splice/sqtls_indep.v3_rn7.txt
#   {outdir}/splice/top_assoc_splice.v3_rn7.txt
#
# Outputs:
#   {outdir}/autocomplete.v3_rn7.json
#   {outdir}/gene.v3_rn7.txt

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissues <- args[3:length(args)]

version <- "v3"
v <- "v3_rn7"
anno <- "GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf"

ensembl_mappings <- read_tsv(
    str_glue("{indir}/ref/GENES_RAT.txt"),
    col_types = cols(SYMBOL = "c", ENSEMBL_ID = "c", .default = "-"),
    skip = 85 # Some lines have '#' in the text, so I can't use '#' as comment char to skip pre-header
) |>
    rename(geneId = SYMBOL,
           geneIdEnsembl = ENSEMBL_ID) |>
    filter(!is.na(geneId),
           !is.na(geneIdEnsembl)) |>
    separate_rows(geneIdEnsembl, sep=";") |>
    group_by(geneId) |>
    slice(1) |>
    ungroup()
# Check for duplicate gene IDs and warn if found
n_dups <- sum(duplicated(ensembl_mappings$geneId))
if (n_dups > 0) {
    warning(str_glue("Found {n_dups} duplicate gene IDs in ensembl_mappings"))
}

# Not sure if this is exact same set as in top_assoc.txt, but this is to see if
# file of signif pairs exists for each gene:
signif <- tibble(tissue = tissues) |>
    reframe(
        read_tsv(str_glue("{indir}/{version}/{tissue}/{tissue}.cis_qtl_signif.txt.gz"),
            col_types = "c---------"
        ),
        .by = tissue
    ) |>
    distinct(phenotype_id) |>
    pull()

## Add expression/eQTL status for gene page

expr <- tibble(tissue = tissues) |>
    reframe(
        read_tsv(str_glue("{indir}/{version}/{tissue}/{tissue}.expr.tpm.bed.gz"),
            col_types = c(
                `#chr` = "-", start = "-", end = "-",
                gene_id = "c", .default = "d"
            )
        ) |>
            pivot_longer(-gene_id, names_to = "rat_id", values_to = "tpm") |>
            group_by(gene_id) |>
            summarise(
                expressed = sum(tpm > 0) > 1,
                .groups = "drop"
            ),
        .by = tissue
    )

expr_all <- expr |>
    filter(expressed) |>
    distinct(gene_id) |>
    pull()

is_expr <- expr |>
    mutate(expressed = if_else(expressed, "True", "False")) |>
    pivot_wider(
        id_cols = gene_id, names_from = tissue, names_prefix = "expr_",
        values_from = expressed, values_fill = "False"
    )

was_tested_eqtl <- read_tsv(str_glue("{outdir}/eqtl/top_assoc.{v}.txt"),
    col_types = "cc----------------"
) |>
    mutate(testedEqtl = "True") |>
    complete(tissue, gene_id = expr$gene_id, fill = list(testedEqtl = "False")) |>
    pivot_wider(
        id_cols = gene_id, names_from = tissue, names_prefix = "testedEqtl_",
        values_from = testedEqtl
    )

has_eqtl <- read_tsv(str_glue("{outdir}/eqtl/eqtls_indep.{v}.txt"),
    col_types = "cc---------------"
) |>
    filter(tissue %in% tissues) |>
    distinct(tissue, gene_id) |>
    mutate(eqtl = "True") |>
    complete(tissue, gene_id = expr$gene_id, fill = list(eqtl = "False")) |>
    pivot_wider(
        id_cols = gene_id, names_from = tissue, names_prefix = "eqtl_",
        values_from = eqtl
    )

## Add alt splicing/sQTL status for gene page

splice <- tibble(tissue = tissues) |>
    reframe(
        read_tsv(str_glue("{indir}/{version}/{tissue}/splice/{tissue}.leafcutter.phenotype_groups.txt"),
            skip = 1, col_names = c("phenotype_id", "gene_id"), col_types = "-c"
        ),
        .by = tissue
    ) |>
    distinct(tissue, gene_id)

is_alt_spliced <- splice |>
    mutate(altSpliced = "True") |>
    complete(tissue, gene_id = expr$gene_id, fill = list(altSpliced = "False")) |>
    pivot_wider(
        id_cols = gene_id, names_from = tissue, names_prefix = "altSplice_",
        values_from = altSpliced
    )

was_tested_sqtl <- read_tsv(str_glue("{outdir}/splice/top_assoc_splice.{v}.txt"),
    col_types = "c-------------c----"
) |>
    distinct(tissue, gene_id) |>
    mutate(testedSqtl = "True") |>
    complete(tissue, gene_id = expr$gene_id, fill = list(testedSqtl = "False")) |>
    pivot_wider(
        id_cols = gene_id, names_from = tissue, names_prefix = "testedSqtl_",
        values_from = testedSqtl
    )

has_sqtl <- read_tsv(str_glue("{outdir}/splice/sqtls_indep.{v}.txt"),
                     col_types = "c-------------c---") |>
    filter(tissue %in% tissues) |>
    distinct(tissue, gene_id) |>
    mutate(sqtl = "True") |>
    complete(tissue, gene_id = expr$gene_id, fill = list(sqtl = "False")) |>
    pivot_wider(
        id_cols = gene_id, names_from = tissue, names_prefix = "sqtl_",
        values_from = sqtl
    )

## Assemble

genes <- read_tsv(
    str_glue("{indir}/ref/{anno}"),
    col_types = "c-cii-c-c",
    col_names = c("chromosome", "type", "start", "end", "strand", "etc"),
    comment = "#"
) |>
    filter(type == "gene") |>
    mutate(
        geneId = str_match(etc, 'gene_id "([^"]+)"')[, 2],
        description = str_match(etc, 'description "([^"]+)"')[, 2],
        tss = if_else(strand == "+", start, end)
    ) |>
    replace_na(list(description = "")) |>
    left_join(ensembl_mappings, by = "geneId", relationship = "one-to-one") |>
    select(geneId, geneIdEnsembl, chromosome, start, end, strand, tss, description) |>
    mutate(hasEqtl = if_else(geneId %in% signif, "True", "False")) |>
    left_join(is_expr, by = c("geneId" = "gene_id"), relationship = "one-to-one") |>
    left_join(was_tested_eqtl, by = c("geneId" = "gene_id"), relationship = "one-to-one") |>
    left_join(has_eqtl, by = c("geneId" = "gene_id"), relationship = "one-to-one") |>
    left_join(is_alt_spliced, by = c("geneId" = "gene_id"), relationship = "one-to-one") |>
    left_join(was_tested_sqtl, by = c("geneId" = "gene_id"), relationship = "one-to-one") |>
    left_join(has_sqtl, by = c("geneId" = "gene_id"), relationship = "one-to-one")

write_tsv(genes, str_glue("{outdir}/gene.{v}.txt"))

# Save list of all (expressed) names and IDs for search autocomplete:
genes2 <- genes |>
    filter(geneId %in% expr_all)
c(genes2$geneId, genes2$geneIdEnsembl) |>
    unique() |>
    sort() |>
    jsonlite::toJSON() |>
    write_file(str_glue("{outdir}/autocomplete.{v}.json"))

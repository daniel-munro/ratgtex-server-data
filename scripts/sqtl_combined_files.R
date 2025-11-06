# Process and aggregate sQTL results from tensorQTL
#
# Inputs:
#   {indir}/geno/alleles.txt.gz
#   {indir}/ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf
#   {indir}/v3/{tissue}/splice/{tissue}_splice.cis_qtl.txt.gz
#   {indir}/v3/{tissue}/splice/{tissue}_splice.cis_independent_qtl.txt.gz
#
# Outputs:
#   {outdir}/splice/top_assoc_splice.v3_rn7.txt
#   {outdir}/splice/sqtls_indep.v3_rn7.txt

suppressPackageStartupMessages(library(tidyverse))

load_tensorqtl <- function(tensorqtl_out) {
    read_tsv(tensorqtl_out, col_types = "ci----c----dddd-dcidd") |>
        rename(gene_id = group_id) |>
        separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
                 remove = FALSE)
}

load_tensorqtl_ind <- function(tensorqtl_out) {
    read_tsv(tensorqtl_out, col_types = "ci----c----dddd-dcii") |>
        rename(gene_id = group_id) |>
        separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
                 remove = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissues <- args[3:length(args)]

version <- "v3"
v <- "v3_rn7"
anno <- "GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf"

genes <- read_tsv(str_glue("{indir}/ref/{anno}"),
                  col_types = "--cii-c-c",
                  col_names = c("type", "start", "end", "strand", "etc"),
                  comment = "#") |>
    filter(type == "gene") |>
    mutate(gene_id = str_match(etc, 'gene_id "([^"]+)"')[, 2],
           tss = if_else(strand == "+", start, end)) |>
    select(gene_id, strand, tss)
stopifnot(sum(duplicated(genes$gene_id)) == 0)

alleles <- read_tsv(str_glue("{indir}/geno/alleles.txt.gz"), col_types = "ccc",
                    col_names = c("variant_id", "ref", "alt"))

top_assoc <- tibble(tissue = tissues) |>
    reframe(
        load_tensorqtl(str_glue("{indir}/{version}/{tissue}/splice/{tissue}_splice.cis_qtl.txt.gz")),
        .by = tissue
    ) |>
    left_join(genes, by = "gene_id", relationship = "many-to-one") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
    select(-strand, -tss) |>
    left_join(alleles, by = "variant_id", relationship = "many-to-one") |>
    relocate(tss_distance, .after = af) |>
    relocate(ref, alt, .after = pos)

write_tsv(top_assoc, str_glue("{outdir}/splice/top_assoc_splice.{v}.txt"))

sqtls_ind <- tibble(tissue = tissues) |>
    reframe(
        load_tensorqtl_ind(
            str_glue("{indir}/{version}/{tissue}/splice/{tissue}_splice.cis_independent_qtl.txt.gz")
        ),
        .by = tissue
    ) |>
    left_join(genes, by = "gene_id", relationship = "many-to-one") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
    select(-strand, -tss) |>
    left_join(alleles, by = "variant_id", relationship = "many-to-one") |>
    relocate(tss_distance, .after = af) |>
    relocate(ref, alt, .after = pos)

write_tsv(sqtls_ind, str_glue("{outdir}/splice/sqtls_indep.{v}.txt"))

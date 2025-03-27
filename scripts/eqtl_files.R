# Process eQTL results from tensorQTL for web interface
#
# Inputs:
#   {indir}/geno/alleles.txt.gz
#   {indir}/ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf
#   {indir}/v3/{tissue}/{tissue}.aFC.txt
#   {indir}/v3/{tissue}/{tissue}.cis_independent_qtl.txt.gz
#   {indir}/v3/{tissue}/{tissue}.cis_qtl.txt.gz
#   {indir}/v3/{tissue}/{tissue}.cis_qtl_signif.txt.gz
#   {indir}/v3/{tissue}/{tissue}.trans_qtl_pairs.txt.gz
#
# Outputs:
#   {outdir}/eqtl/cis_qtl_signif.{tissue}.v3_rn7.txt.gz
#   {outdir}/eqtl/eqtls_indep.{tissue}.v3_rn7.txt
#   {outdir}/eqtl/top_assoc.{tissue}.v3_rn7.txt
#   {outdir}/eqtl/trans_qtl_pairs.{tissue}.v3_rn7.txt.gz

suppressPackageStartupMessages(library(tidyverse))

load_tensorqtl <- function(tensorqtl_out) {
    read_tsv(tensorqtl_out, col_types = "ci----c----dddd-ddd") |>
        rename(gene_id = phenotype_id) |>
        separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
                 remove = FALSE)
}

load_tensorqtl_ind <- function(tensorqtl_out) {
    read_tsv(tensorqtl_out, col_types = "ci----c----dddd-cd") |>
        rename(gene_id = phenotype_id) |>
        separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
                 remove = FALSE)
}

load_afc <- function(afc_out) {
    read_tsv(afc_out, col_types = "cc--d--") |>
        rename(gene_id = pid,
               variant_id = sid) |>
        mutate(log2_aFC = signif(log2_aFC, 6))
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
           gene_name = str_match(etc, 'gene_name "([^"]+)"')[, 2],
           tss = if_else(strand == "+", start, end)) |>
    select(gene_id, gene_name, strand, tss)
stopifnot(sum(duplicated(genes$gene_id)) == 0)

alleles <- read_tsv(str_glue("{indir}/geno/alleles.txt.gz"), col_types = "ccc",
                    col_names = c("variant_id", "ref", "alt"))

afc <- tibble(tissue = tissues) |>
    reframe(
        load_afc(str_glue("{indir}/{version}/{tissue}/{tissue}.aFC.txt")),
        .by = tissue
    )

top_assoc <- tibble(tissue = tissues) |>
    reframe(
        load_tensorqtl(str_glue("{indir}/{version}/{tissue}/{tissue}.cis_qtl.txt.gz")),
        .by = tissue
    ) |>
    left_join(afc, by = c("tissue", "gene_id", "variant_id"), relationship = "one-to-one") |>
    left_join(genes, by = "gene_id", relationship = "many-to-one") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
    select(-strand, -tss) |>
    left_join(alleles, by = "variant_id", relationship = "many-to-one") |>
    relocate(gene_name, .after = gene_id) |>
    relocate(tss_distance, .after = af) |>
    relocate(ref, alt, .after = pos)

write_tsv(top_assoc, str_glue("{outdir}/eqtl/top_assoc.{v}.txt"))

eqtls_ind <- tibble(tissue = tissues) |>
    reframe(
        load_tensorqtl_ind(
            str_glue("{indir}/{version}/{tissue}/{tissue}.cis_independent_qtl.txt.gz")
        ),
        .by = tissue
    ) |>
    left_join(afc, by = c("tissue", "gene_id", "variant_id"), relationship = "one-to-one") |>
    left_join(genes, by = "gene_id", relationship = "many-to-one") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
    select(-strand, -tss) |>
    left_join(alleles, by = "variant_id", relationship = "many-to-one") |>
    relocate(gene_name, .after = gene_id) |>
    relocate(tss_distance, .after = af) |>
    relocate(ref, alt, .after = pos)

write_tsv(eqtls_ind, str_glue("{outdir}/eqtl/eqtls_indep.{v}.txt"))

##################################
## Copy and modify signif files ##
##################################

for (tissue in tissues) {
    infile <- str_glue("{indir}/{version}/{tissue}/{tissue}.cis_qtl_signif.txt.gz")
    outfile <- str_glue("{outdir}/eqtl/cis_qtl_signif.{tissue}.{v}.txt.gz")
    cat(str_glue("Making {outfile}"), "\n", sep = "")
    read_tsv(infile, col_types = "cccccccccc") |>
        rename(gene_id = phenotype_id) |>
        separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
                 remove = FALSE) |>
        left_join(select(genes, gene_id, strand, tss), by = "gene_id", relationship = "many-to-one") |>
        mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
        select(-chrom, -pos, -strand, -tss) |>
        write_tsv(outfile)
}

#################################
## Copy and modify trans files ##
#################################

for (tissue in tissues) {
    infile <- str_glue("{indir}/{version}/{tissue}/{tissue}.trans_qtl_pairs.txt.gz")
    outfile <- str_glue("{outdir}/eqtl/trans_qtl_pairs.{tissue}.{v}.txt.gz")
    cat(str_glue("Making {outfile}"), "\n", sep = "")
    read_tsv(infile, col_types = "cccccc") |>
        rename(gene_id = phenotype_id,
               slope = b,
               slope_se = b_se) |>
        write_tsv(outfile)
}

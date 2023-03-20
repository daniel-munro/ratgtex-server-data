suppressPackageStartupMessages(library(tidyverse))

load_tensorqtl <- function(tensorqtl_out) {
    read_tsv(tensorqtl_out, col_types = "ci----c---dddd-ddd") |>
        rename(gene_id = phenotype_id) |>
        separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
                 remove = FALSE) |>
        mutate(chrom = str_replace(chrom, "chr", ""))
}

load_tensorqtl_ind <- function(tensorqtl_out) {
    read_tsv(tensorqtl_out, col_types = "ci----c---dddd-cd") |>
        rename(gene_id = phenotype_id) |>
        separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
                 remove = FALSE) |>
        mutate(chrom = str_replace(chrom, "chr", ""))
}

load_afc <- function(afc_out) {
    read_tsv(afc_out, col_types = "cc--d--") |>
        rename(gene_id = pid,
               variant_id = sid) |>
        mutate(log2_aFC = signif(log2_aFC, 6))
}

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
rn <- args[2]
outdir <- args[3]
tissues <- args[4:length(args)]

anno <- c(
    rn6 = "Rattus_norvegicus.Rnor_6.0.99.genes.bed",
    rn7 = "Rattus_norvegicus.mRatBN7.2.108.genes.bed"
)[rn]

genes <- read_tsv(str_glue("{indir}/ref_{rn}/{anno}"),
                  col_types = "-iic-c---c",
                  col_names = c("start", "end", "gene_id", "strand", "etc")) |>
    mutate(gene_name = str_match(etc, 'gene_name "([^"]+)"')[, 2],
           tss = if_else(strand == "+", start, end)) |>
    select(gene_id, gene_name, strand, tss)

alleles <- read_tsv(str_glue("{indir}/geno_{rn}/alleles.txt.gz"), col_types = "ccc",
                    col_names = c("variant_id", "ref", "alt"))

afc <- tibble(tissue = tissues) |>
    group_by(tissue) |>
    summarise(
        load_afc(str_glue("{indir}/{rn}/{tissue}/{tissue}.aFC.txt")),
        .groups = "drop"
    )

top_assoc <- tibble(tissue = tissues) |>
    group_by(tissue) |>
    summarise(
        load_tensorqtl(str_glue("{indir}/{rn}/{tissue}/{tissue}.cis_qtl.txt.gz")),
        .groups = "drop"
    ) |>
    left_join(afc, by = c("tissue", "gene_id", "variant_id")) |>
    left_join(genes, by = "gene_id") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
    select(-strand, -tss) |>
    left_join(alleles, by = "variant_id") |>
    relocate(gene_name, .after = gene_id) |>
    relocate(tss_distance, .after = af) |>
    relocate(ref, alt, .after = pos)

write_tsv(top_assoc, str_glue("{outdir}/eqtl/{rn}.top_assoc.txt"))

eqtls_ind <- tibble(tissue = tissues) |>
    group_by(tissue) |>
    summarise(
        load_tensorqtl_ind(
            str_glue("{indir}/{rn}/{tissue}/{tissue}.cis_independent_qtl.txt.gz")
        ),
        .groups = "drop"
    ) |>
    left_join(afc, by = c("tissue", "gene_id", "variant_id")) |>
    left_join(genes, by = "gene_id") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
    select(-strand, -tss) |>
    left_join(alleles, by = "variant_id") |>
    relocate(gene_name, .after = gene_id) |>
    relocate(tss_distance, .after = af) |>
    relocate(ref, alt, .after = pos)

write_tsv(eqtls_ind, str_glue("{outdir}/eqtl/{rn}.eqtls_indep.txt"))

##################################
## Copy and modify signif files ##
##################################

for (tissue in tissues) {
    fname <- str_glue("{outdir}/eqtl/{tissue}.{rn}.cis_qtl_signif.txt.gz")
    cat(str_glue("Making {fname}"), "\n", sep = "")
    read_tsv(str_glue("{indir}/{rn}/{tissue}/{tissue}.cis_qtl_signif.txt.gz"),
             col_types = "cccccccccc") |>
        rename(gene_id = phenotype_id) |>
        separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
                 remove = FALSE) |>
        left_join(select(genes, gene_id, strand, tss), by = "gene_id") |>
        mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
        select(-chrom, -pos, -strand, -tss) |>
        write_tsv(fname)
}

#################################
## Copy and modify trans files ##
#################################

for (tissue in tissues) {
    fname <- str_glue("{outdir}/eqtl/{tissue}.{rn}.trans_qtl_pairs.txt.gz")
    cat(str_glue("Making {fname}"), "\n", sep = "")
    read_tsv(str_glue("{indir}/{rn}/{tissue}/{tissue}.trans_qtl_pairs.txt.gz"),
             col_types = "cccccc") |>
        rename(gene_id = phenotype_id,
               slope = b,
               slope_se = b_se) |>
        write_tsv(fname)
}

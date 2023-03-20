suppressPackageStartupMessages(library(tidyverse))

load_tensorqtl <- function(tensorqtl_out) {
    read_tsv(tensorqtl_out, col_types = "ci----c---dddd-dcidd") |>
        rename(gene_id = group_id) |>
        separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
                 remove = FALSE) |>
        mutate(chrom = str_replace(chrom, "chr", ""))
}

load_tensorqtl_ind <- function(tensorqtl_out) {
    read_tsv(tensorqtl_out, col_types = "ci----c---dddd-dcii") |>
        rename(gene_id = group_id) |>
        separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
                 remove = FALSE) |>
        mutate(chrom = str_replace(chrom, "chr", ""))
}

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
rn <- args[2]
outdir <- args[3]
tissues <- args[4:length(args)]

anno <- c(
    rn6 = "Rattus_norvegicus.Rnor_6.0.99.genes.gtf",
    rn7 = "Rattus_norvegicus.mRatBN7.2.108.genes.gtf"
)[rn]

genes <- read_tsv(str_glue("{indir}/ref_{rn}/{anno}"),
                  col_types = "--cii-c-c",
                  col_names = c("type", "start", "end", "strand", "etc"),
                  comment = "#") |>
    filter(type == "gene") |>
    mutate(gene_id = str_match(etc, 'gene_id "([^"]+)"')[, 2],
           gene_name = str_match(etc, 'gene_name "([^"]+)"')[, 2],
           tss = if_else(strand == "+", start, end)) |>
    select(gene_id, gene_name, strand, tss)
stopifnot(sum(duplicated(genes$gene_id)) == 0)

alleles <- read_tsv(str_glue("{indir}/geno_{rn}/alleles.txt.gz"), col_types = "ccc",
                    col_names = c("variant_id", "ref", "alt"))

top_assoc <- tibble(tissue = tissues) |>
    group_by(tissue) |>
    summarise(
        load_tensorqtl(str_glue("{indir}/{rn}/{tissue}/splice/{tissue}_splice.cis_qtl.txt.gz")),
        .groups = "drop"
    ) |>
    left_join(genes, by = "gene_id") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
    select(-strand, -tss) |>
    left_join(alleles, by = "variant_id") |>
    relocate(gene_name, .after = gene_id) |>
    relocate(tss_distance, .after = af) |>
    relocate(ref, alt, .after = pos)

write_tsv(top_assoc, str_glue("{outdir}/splice/{rn}.top_assoc_splice.txt"))

sqtls_ind <- tibble(tissue = tissues) |>
    group_by(tissue) |>
    summarise(
        load_tensorqtl_ind(
            str_glue("{indir}/{rn}/{tissue}/splice/{tissue}_splice.cis_independent_qtl.txt.gz")
        ),
        .groups = "drop"
    ) |>
    left_join(genes, by = "gene_id") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
    select(-strand, -tss) |>
    left_join(alleles, by = "variant_id") |>
    relocate(gene_name, .after = gene_id) |>
    relocate(tss_distance, .after = af) |>
    relocate(ref, alt, .after = pos)

write_tsv(sqtls_ind, str_glue("{outdir}/splice/{rn}.sqtls_indep.txt"))

#################################
## Copy and modify trans files ##
#################################

for (tissue in tissues) {
    fname <- str_glue("{outdir}/splice/{tissue}.{rn}.splice.trans_qtl_pairs.txt.gz")
    cat(str_glue("Making {fname}"), "\n", sep = "")
    read_tsv(str_glue("{indir}/{rn}/{tissue}/splice/{tissue}_splice.trans_qtl_pairs.txt.gz"),
             col_types = "cccccc") |>
        rename(slope = b,
               slope_se = b_se) |>
        write_tsv(fname)
}

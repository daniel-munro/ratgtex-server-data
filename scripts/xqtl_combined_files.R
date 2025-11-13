# Process and aggregate xQTL results from tensorQTL

suppressPackageStartupMessages(library(tidyverse))

load_assoc_ungrouped <- function(tensorqtl_out) {
  read_tsv(tensorqtl_out, col_types = "ci----c----cccc-ccc") |>
    mutate(gene_id = phenotype_id)
}

load_assoc_grouped <- function(tensorqtl_out) {
  read_tsv(tensorqtl_out, col_types = "ci----c----cccc-cc-cc") |>
    rename(gene_id = group_id)
}

load_indep_ungrouped <- function(tensorqtl_out) {
  read_tsv(tensorqtl_out, col_types = "ci----c----cccc-ci") |>
    mutate(gene_id = phenotype_id)
}

load_indep_grouped <- function(tensorqtl_out) {
  read_tsv(tensorqtl_out, col_types = "ci----c----cccc-cc-i") |>
    rename(gene_id = group_id)
}

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
modality <- args[3]
tissues <- args[4:length(args)]

input_genes <- str_glue("{outdir}/gene.v4.tsv")
input_alleles <- str_glue("{indir}/geno/alleles.tsv.gz")
input_xqtl_assoc_template <- "{indir}/v4/{tissue}/pheast/output/qtl/{modality}.cis_qtl.txt.gz"
input_xqtl_indep_template <- "{indir}/v4/{tissue}/pheast/output/qtl/{modality}.cis_independent_qtl.txt.gz"

output_assoc <- str_glue("{outdir}/xqtl/top_assoc.{modality}.v4_rn8.tsv")
output_indep <- str_glue("{outdir}/xqtl/xqtls_indep.{modality}.v4_rn8.tsv")

load_assoc <- if (modality %in% c("expression", "stability")) load_assoc_ungrouped else load_assoc_grouped
load_indep <- if (modality %in% c("expression", "stability")) load_indep_ungrouped else load_indep_grouped

genes <- read_tsv(
  input_genes,
  col_types = cols(geneId = "c", strand = "c", tss = "i", .default = "-")
) |>
  rename(gene_id = geneId)

alleles <- read_tsv(
  input_alleles,
  col_types = "ccc",
  col_names = c("variant_id", "ref", "alt")
)

top_assoc <- tibble(tissue = tissues) |>
  reframe(
    load_assoc(str_glue(input_xqtl_assoc_template)),
    .by = tissue
  ) |>
  separate(
    variant_id, c("chrom", "pos"), sep = ":", convert = TRUE, remove = FALSE
  ) |>
  left_join(genes, by = "gene_id", relationship = "many-to-one") |>
  mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
  left_join(alleles, by = "variant_id", relationship = "many-to-one") |>
  select(
    tissue, phenotype_id, gene_id, num_var, variant_id, chrom, pos, ref, alt,
    af, tss_distance, pval_nominal, slope, slope_se, pval_beta, qval,
    pval_nominal_threshold
  )

if (modality == "cross_modality") {
  top_assoc <- top_assoc |>
    separate_wider_delim(phenotype_id,  ":", names = c("modality", "phenotype_id"))
}

write_tsv(top_assoc, output_assoc)

xqtls_ind <- tibble(tissue = tissues) |>
  reframe(
    load_indep(str_glue(input_xqtl_indep_template)),
    .by = tissue
  ) |>
  separate(
    variant_id, c("chrom", "pos"), sep = ":", convert = TRUE, remove = FALSE
  ) |>
  left_join(genes, by = "gene_id", relationship = "many-to-one") |>
  mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
  left_join(alleles, by = "variant_id", relationship = "many-to-one") |>
  select(
    tissue, phenotype_id, gene_id, num_var, variant_id, chrom, pos, ref, alt,
    af, tss_distance, pval_nominal, slope, slope_se, pval_beta, rank
  )

if (modality == "cross_modality") {
  xqtls_ind <- xqtls_ind |>
    separate_wider_delim(phenotype_id,  ":", names = c("modality", "phenotype_id"))
}

write_tsv(xqtls_ind, output_indep)

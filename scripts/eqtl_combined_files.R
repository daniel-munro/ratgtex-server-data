# Process and aggregate eQTL results from tensorQTL

suppressPackageStartupMessages(library(tidyverse))

load_tensorqtl <- function(tensorqtl_out) {
  read_tsv(tensorqtl_out, col_types = "ci----c----dddd-ddd") |>
    rename(gene_id = phenotype_id) |>
    separate(variant_id, c("chrom", "pos"),
      sep = ":", convert = TRUE,
      remove = FALSE
    )
}

load_tensorqtl_ind <- function(tensorqtl_out) {
  read_tsv(tensorqtl_out, col_types = "ci----c----dddd-cd") |>
    rename(gene_id = phenotype_id) |>
    separate(variant_id, c("chrom", "pos"),
      sep = ":", convert = TRUE,
      remove = FALSE
    )
}

load_afc <- function(afc_out) {
  read_tsv(afc_out, col_types = "cc--d--") |>
    rename(
      gene_id = pid,
      variant_id = sid
    ) |>
    mutate(log2_aFC = signif(log2_aFC, 6))
}

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissues <- args[3:length(args)]

input_anno <- str_glue("{indir}/ref/GCF_036323735.1_GRCr8_genomic.chr.gtf")
input_alleles <- str_glue("{indir}/geno/alleles.txt.gz")
input_afc_template <- "{indir}/v4/{tissue}/{tissue}.aFC.txt"
input_eqtl_assoc_template <- "{indir}/v4/{tissue}/pheast/output/qtl/expression.cis_qtl.txt.gz"
input_eqtl_indep_template <- "{indir}/v4/{tissue}/pheast/output/qtl/expression.cis_independent_qtl.txt.gz"

output_assoc <- str_glue("{outdir}/eqtl/top_assoc.v4_rn8.tsv")
output_indep <- str_glue("{outdir}/eqtl/eqtls_indep.v4_rn8.tsv")

genes <- read_tsv(input_anno,
  col_types = "--cii-c-c",
  col_names = c("type", "start", "end", "strand", "etc"),
  comment = "#"
) |>
  filter(type == "gene") |>
  mutate(
    gene_id = str_match(etc, 'gene_id "([^"]+)"')[, 2],
    tss = if_else(strand == "+", start, end)
  ) |>
  select(gene_id, strand, tss)
stopifnot(sum(duplicated(genes$gene_id)) == 0)

alleles <- read_tsv(input_alleles,
  col_types = "ccc",
  col_names = c("variant_id", "ref", "alt")
)

afc <- tibble(tissue = tissues) |>
  reframe(
    load_afc(str_glue(input_afc_template)),
    .by = tissue
  )

top_assoc <- tibble(tissue = tissues) |>
  reframe(
    load_tensorqtl(str_glue(input_eqtl_assoc_template)),
    .by = tissue
  ) |>
  left_join(afc, by = c("tissue", "gene_id", "variant_id"), relationship = "one-to-one") |>
  left_join(genes, by = "gene_id", relationship = "many-to-one") |>
  mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
  select(-strand, -tss) |>
  left_join(alleles, by = "variant_id", relationship = "many-to-one") |>
  relocate(tss_distance, .after = af) |>
  relocate(ref, alt, .after = pos)

write_tsv(top_assoc, output_assoc)

eqtls_ind <- tibble(tissue = tissues) |>
  reframe(
    load_tensorqtl_ind(str_glue(input_eqtl_indep_template)),
    .by = tissue
  ) |>
  left_join(afc, by = c("tissue", "gene_id", "variant_id"), relationship = "one-to-one") |>
  left_join(genes, by = "gene_id", relationship = "many-to-one") |>
  mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
  select(-strand, -tss) |>
  left_join(alleles, by = "variant_id", relationship = "many-to-one") |>
  relocate(tss_distance, .after = af) |>
  relocate(ref, alt, .after = pos)

write_tsv(eqtls_ind, output_indep)

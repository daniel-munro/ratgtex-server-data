suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissue <- args[3]

d <- read_tsv(
    str_glue("{indir}/{tissue}/{tissue}.cis_qtl_all_pvals.txt.gz"),
    col_types = "ccd"
)
d <- split(d, ~ phenotype_id)
for (gene in names(d)) {
    d[[gene]] |>
        select(-phenotype_id) |>
        write_tsv(str_glue("{outdir}/cis_pvals/{tissue}/{gene}.txt"))
}

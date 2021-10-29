# library(RSQLite)
library(tidyverse)

for (tissue in c("Eye", "IL", "LHb", "NAcc", "OFC", "PL")) {
    cat(tissue, "\n", sep = "")
    d <- read_tsv(
        # str_glue("data/tensorqtl/all_cis_pvals/{str_sub(tissue, 1, 1)}QCT.all_cis_pvals.txt.gz"),
        str_glue("data/tensorqtl/{tissue}.cis_qtl_all_pvals.txt.gz"),
        col_types = "ccd"
    ) %>%
        rename(gene_id = phenotype_id)
    # con <- dbConnect(SQLite(), dbname = str_glue("~/code/portal_data/cis_pvals/{tissue}.db"))
    # dbWriteTable(con, "pvals", d)
    # dbDisconnect(con)
    d <- split(d, ~ gene_id)
    for (gene in names(d)) {
        d[[gene]] %>%
            select(-gene_id) %>%
            write_tsv(str_glue("~/code/portal_data/cis_pvals/{tissue}/{gene}.txt"))
    }
    rm(d); gc()
}

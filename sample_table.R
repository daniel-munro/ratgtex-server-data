suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissues <- args[3:length(args)]

# I downloaded GSE173141_series_matrix.txt.gz and edited it to make it readable.
geo <- read_tsv(str_glue("{indir}/samples/GSE173141_series_matrix.txt.gz"), comment = "!", 
                col_types = cols(.default = "c")) |>
    rename(field = `Sample_title`) |>
    pivot_longer(-field, names_to = "sample_id") |>
    pivot_wider(sample_id, names_from = field, values_from = value) |>
    separate(sample_id, c("rat_id", "tissue")) |>
    mutate(tissue = c(IL = "IL", LHB = "LHb", Acbc = "NAcc", VoLo = "OFC", PL = "PL")[tissue],
           BioSample = str_replace(
               BioSample,
               "BioSample: https://www.ncbi.nlm.nih.gov/biosample/",
               ""
           ),
           SRA = str_replace(
               SRA,
               fixed("SRA: https://www.ncbi.nlm.nih.gov/sra?term="),
               ""
           ),
           sex = str_replace(sex, "Sex: ", ""),
           rat_batch = str_replace(rat_batch, "rat batch: ", "")
    )

geo2 <- geo |>
    select(tissue,
           rat_id,
           GEO_accession = Sample_geo_accession,
           BioSample_accession = BioSample,
           SRA_accession = SRA)

geo_rats <- geo |>
    distinct(rat_id, sex, rat_batch)

samples <- tibble(tissue = tissues) |>
    group_by(tissue) |>
    summarise(
        tibble(rat_id = read_lines(str_glue("{indir}/{tissue}/rat_ids.txt"))),
        .groups = "drop"
    ) |>
    left_join(geo2, by = c("tissue", "rat_id"))

rats <- samples |>
    group_by(rat_id) |>
    summarise(tissues = str_c(tissue, collapse = ", ")) |>
    left_join(geo_rats, by = "rat_id")

write_tsv(samples, str_glue("{outdir}/ref/RatGTEx_samples.tsv"))

write_tsv(rats, str_glue("{outdir}/ref/RatGTEx_rats.tsv"))

###################
## Write to HTML ##
###################
# Copy and paste these files' contents into about/samples/index.html.

html_sam <- '<table id="sample-table"><thead><tr><th>Rat ID</th><th>Tissue</th><th>GEO Accession</th><th>BioSample Accession</th><th>SRA Accession</th></tr></thead><tbody>'
for (r in 1:nrow(samples)) {
    html_sam <- html_sam |> str_c(str_glue('<tr><td>{samples$rat_id[r]}</td><td>{samples$tissue[r]}</td>'))
    html_sam <- html_sam |> str_c(if (is.na(samples$GEO_accession[r])) {
        '<td>-</td>'
    } else {
        str_glue('<td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={samples$GEO_accession[r]}">{samples$GEO_accession[r]}</a></td>')
    })
    html_sam <- html_sam |> str_c(if (is.na(samples$BioSample_accession[r])) {
        '<td>-</td>'
    } else {
        str_glue('<td><a href="https://www.ncbi.nlm.nih.gov/biosample/{samples$BioSample_accession[r]}">{samples$BioSample_accession[r]}</a></td>')
    })
    html_sam <- html_sam |> str_c(if (is.na(samples$SRA_accession[r])) {
        '<td>-</td>'
    } else {
        str_glue('<td><a href="https://www.ncbi.nlm.nih.gov/sra?term={samples$SRA_accession[r]}">{samples$SRA_accession[r]}</a></td></tr>')
    })
}
html_sam <- html_sam |> str_c('</tbody></table>\n')

write_file(html_sam, str_glue("{outdir}/ref/samples.html"))

html_rat <- '<table id="rat-table"><thead><tr><th>Rat ID</th><th>Sex</th><th>Tissues</th></tr></thead><tbody>'
for (r in 1:nrow(rats)) {
    id <- rats$rat_id[r]
    html_rat <- html_rat |> str_c(str_glue('<tr><td>{id}</td>'))
    html_rat <- html_rat |> str_c(if (is.na(rats$sex[r])) {
        '<td>-</td>'
    } else {
        str_glue('<td>{rats$sex[r]}</td>')
    })
    tissue_text <- c()
    for (tis in str_split(rats$tissues[r], ", ")[[1]]) {
        sam <- filter(samples, tissue == tis, rat_id == id)
        tissue_text <- tissue_text |> c(if (is.na(sam$GEO_accession)) {
            str_glue('{tis}')
        } else {
            str_glue('<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={sam$GEO_accession}">{tis}</a>')
        })
    }
    html_rat <- html_rat |> str_c("<td>", str_c(tissue_text, collapse = ", "), "</td>")
}
html_rat <- html_rat |> str_c('</tbody></table>\n')

write_file(html_rat, str_glue("{outdir}/ref/rats.html"))

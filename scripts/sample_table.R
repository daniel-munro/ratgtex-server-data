# Generate sample and rat tables for web interface
#
# Inputs:
#   {indir}/geno/genotyping_log.csv
#   {indir}/samples/{geo_accession}_series_matrix.txt.gz
#   {indir}/{version}/{tissue}/rat_ids.txt
#
# Outputs:
#   {outdir}/ref/RatGTEx_rats.{v}.html
#   {outdir}/ref/RatGTEx_rats.{v}.tsv
#   {outdir}/ref/RatGTEx_samples.{v}.html
#   {outdir}/ref/RatGTEx_samples.{v}.tsv

suppressPackageStartupMessages(library(tidyverse))

field_values <- function(lines, field_name, keyword = NULL) {
    # Some field names are not unique, so a keyword can be supplied that must
    # be present anywhere in the line, e.g. "BioSample" or "SRA"
    if (!is.null(keyword)) {
        lines <- str_subset(lines, keyword)
    }
    lines |>
        str_subset(str_glue("^!{field_name}\\t")) |>
        str_replace(str_glue("^!{field_name}\\t"), "") |>
        str_replace_all('"', '') |>
        str_split(pattern = "\\t", simplify = TRUE) |>
        c()
}

load_geo <- function(filename, trim_ids = FALSE) {
    lines <- read_lines(filename)
    geo <- tibble(
        rat_id = field_values(lines, "Sample_title"),
        GEO_accession = field_values(lines, "Sample_geo_accession"),
        BioSample_accession = field_values(lines, "Sample_relation", "BioSample"),
        SRA_accession = field_values(lines, "Sample_relation", "SRA"),
    ) |>
        mutate(
            BioSample_accession = str_replace(
                BioSample_accession,
                "BioSample: https://www.ncbi.nlm.nih.gov/biosample/",
                ""
            ),
            SRA_accession = str_replace(
                SRA_accession,
                fixed("SRA: https://www.ncbi.nlm.nih.gov/sra?term="),
                ""
            ),
        )
    if (trim_ids) { # Remove tissue name from IDs in older submissions
        geo <- geo |>
            mutate(rat_id = str_replace(rat_id, "_.+$", ""))
    }
    geo
}

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
tissues <- args[3:length(args)]

version <- "v3"
v <- "v3_rn7"

accession <- c(
    Eye = "GSE201236",
    IL = "GSE173137",
    LHb = "GSE173138",
    NAcc = "GSE173136",
    OFC = "GSE173140",
    PL = "GSE173139"
)
to_trim <- c("IL", "LHb", "NAcc", "OFC", "PL") # IDs have tissue appended to rat ID

geo <- tibble(tissue = tissues) |>
    filter(tissue %in% names(accession)) |>
    reframe(
        load_geo(
            str_glue("{indir}/samples/{accession[tissue]}_series_matrix.txt.gz"),
            tissue %in% to_trim
        ),
        .by = tissue
    )

samples <- tibble(tissue = tissues) |>
    reframe(
        tibble(rat_id = read_lines(str_glue("{indir}/{version}/{tissue}/rat_ids.txt"))),
        .by = tissue
    ) |>
    left_join(geo, by = c("tissue", "rat_id"), relationship = "one-to-one")

rats <- samples |>
    group_by(rat_id) |>
    summarise(tissues = str_c(tissue, collapse = ", "))

meta <- read_csv(
    str_glue("{indir}/geno/genotyping_log.csv"),
    col_types = cols(rfid = "c", sex = "c", coatcolor = "c", .default = "-")
) |>
    rename(rat_id = rfid) |>
    mutate(
        coatcolor = coatcolor |> # Collapse equivalent labels
            str_to_upper() |>
            str_replace("BLK HOOD", "BLACKHOOD") |>
            str_replace("BRN HOOD", "BROWNHOOD")
    )

rats <- rats |>
    left_join(meta, by = "rat_id", relationship = "one-to-one") |>
    select(rat_id, sex, coatcolor, tissues)

# Adipose and Liver rats don't have sex/coat metadata, but are all males re correspondence with PI
rats_adi_liv <- samples |>
    filter(tissue %in% c("Adipose", "Liver")) |>
    distinct(rat_id) |>
    pull()
rats <- rats |>
    mutate(sex = if_else(is.na(sex) & rat_id %in% rats_adi_liv, "M", sex))

write_tsv(samples, str_glue("{outdir}/ref/RatGTEx_samples.{v}.tsv"))

write_tsv(rats, str_glue("{outdir}/ref/RatGTEx_rats.{v}.tsv"))

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

write_file(html_sam, str_glue("{outdir}/ref/RatGTEx_samples.{v}.html"))

html_rat <- '<table id="rat-table"><thead><tr><th>Rat ID</th><th>Sex</th><th>Coat Color</th><th>Tissues</th></tr></thead><tbody>'
for (r in 1:nrow(rats)) {
    id <- rats$rat_id[r]
    html_rat <- html_rat |> str_c(str_glue('<tr><td>{id}</td>'))
    html_rat <- html_rat |> str_c(if (is.na(rats$sex[r])) {
        '<td>-</td>'
    } else {
        str_glue('<td>{rats$sex[r]}</td>')
    })
    html_rat <- html_rat |> str_c(if (is.na(rats$coatcolor[r])) {
        '<td>-</td>'
    } else {
        str_glue('<td>{rats$coatcolor[r]}</td>')
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

write_file(html_rat, str_glue("{outdir}/ref/RatGTEx_rats.{v}.html"))

from pathlib import Path

configfile: "../ratgtex/config.yaml"

tissues_separate = sorted([tissue for tissue in config["tissues"] if Path(f"../ratgtex/v4/{tissue}").exists()])
tissues_all = sorted(tissues_separate + list(config["merged_tissues"].keys()))
tissues_merged = sorted([tissue for tissue in tissues_all if tissue not in sum(config["merged_tissues"].values(), [])])

modalities = ["alt_polyA", "alt_TSS", "expression", "isoforms", "splicing", "stability"]
modalities_cross = modalities + ["cross_modality"]

print(f"Separate tissues: {tissues_separate}")
print(f"Merged tissues: {tissues_merged}")
geno_datasets = list(set(config["tissues"][tissue]["geno_dataset"] for tissue in tissues_separate))

localrules:
    exon_table,
    sample_table,
    copy_meta_files,
    copy_covar_files,
    copy_geno_files,
    copy_expr_log2,
    copy_expr_tpm,
    copy_phenos_unnorm,
    copy_phenos_norm,
    median_gene_expression,

rule all:
    input:
        "data/tissue_info.v4.tsv",
        ##
        "data/gene.v4.tsv",
        "data/autocomplete.v4.json",
        ##
        "data/exon.v4.tsv",
        ##
        "data/ref/RatGTEx_rats.v4.tsv",
        "data/ref/RatGTEx_samples.v4.tsv",
        "data/ref/rats.v4.html",
        "data/ref/samples.v4.html",
        ##
        expand("data/covar/covar.{tissue}.{modality}.v4.tsv", tissue=tissues_merged, modality=modalities_cross),
        expand("data/fastq_map/fastq_map.{tissue}.v4.txt", tissue=tissues_merged),
        expand("data/rat_ids/rat_ids.{tissue}.v4.txt", tissue=tissues_merged),
        ##
        expand("data/geno/{geno_dataset}.rn8.{ext}", geno_dataset=geno_datasets, ext=["vcf.gz", "vcf.gz.tbi"]),
        ##
        expand("data/expr/expr.{units}.{tissue}.v4_rn8.bed.gz", units=["log2", "tpm"], tissue=tissues_merged),
        ##
        "data/expr/medianGeneExpression.v4.tsv.gz",
        "data/expr/topExpressedGene.v4.tsv",
        ##
        expand("data/phenos/phenos.{tissue}.{modality}.{norm}.v4_rn8.bed.gz", tissue=tissues_merged, modality=modalities, norm=["unnorm", "norm"]),
        ##
        "data/eqtl/top_assoc.v4_rn8.tsv",
        "data/eqtl/eqtls_indep.v4_rn8.tsv",
        ##
        "data/eqtl/singleTissueEqtl.v4.zip",
        ##
        expand("data/xqtl/top_assoc.{modality}.v4_rn8.tsv", modality=modalities + ["cross_modality"]),
        expand("data/xqtl/xqtls_indep.{modality}.v4_rn8.tsv", modality=modalities + ["cross_modality"]),
        ##
        expand("data/xqtl/cis_qtl_signif.{tissue}.{modality}.v4_rn8.txt.gz", tissue=tissues_merged, modality=modalities),
        expand("data/xqtl/trans_qtl_pairs.{tissue}.expression.v4_rn8.txt.gz", tissue=tissues_merged),
        ##
        expand("data/cis_pvals/{tissue}.v4.zip", tissue=tissues_merged),

rule create_directories:
    """Create directories"""
    output:
        directory("data/cis_pvals"),
        directory("data/covar"),
        directory("data/eqtl"),
        directory("data/expr"),
        directory("data/fastq_map"),
        directory("data/geno"),
        directory("data/rat_ids"),
        directory("data/ref"),
        directory("data/splice")
    shell:
        "mkdir -p {output}"

############################
## Metadata and reference ##
############################

rule tissue_info:
    """Process tissue info table"""
    input:
        info = "../ratgtex/tissue_info.tsv",
        expr = expand("../ratgtex/v4/{tissue}/phenos/output/expression.bed.gz", tissue=tissues_merged),
        eqtl = "data/eqtl/top_assoc.v4_rn8.txt",
    output:
        info = "data/tissue_info.v4.tsv"
    params:
        tissues = tissues_merged
    shell:
        "Rscript scripts/tissue_info.R ../ratgtex data {params.tissues}"

rule gene_info:
    """Process gene info table"""
    input:
        gtf = "../ratgtex/ref/GCF_036323735.1_GRCr8_genomic.chr.gtf",
        signif = expand("../ratgtex/v4/{tissue}/{tissue}.cis_qtl_signif.txt.gz", tissue=tissues_merged),
        expr = expand("../ratgtex/v4/{tissue}/phenos/output/expression.bed.gz", tissue=tissues_merged),
        eqtl = "data/eqtl/eqtls_indep.v4_rn8.txt",
        expr_assoc = "data/eqtl/top_assoc.v4_rn8.txt",
    output:
        gene = "data/gene.v4.tsv",
        autocomplete = "data/autocomplete.v4.json",
    params:
        tissues = tissues_merged
    shell:
        "Rscript scripts/gene.R ../ratgtex data {params.tissues}"

rule exon_table:
    """Prepare exon annotations for server API"""
    input:
        gtf = "../ratgtex/ref/GCF_036323735.1_GRCr8_genomic.chr.genes.gtf"
    output:
        exon = "data/exon.v4.tsv"
    shell:
        "python3 scripts/exon.py {input.gtf} {output.exon}"

rule sample_table:
    """Generate sample and rat tables for web interface"""
    input:
        geno_log = "../ratgtex/geno/genotyping_log.csv",
        geo = expand(
            "../ratgtex/samples/{accession}_series_matrix.txt.gz",
            accession=["GSE201236", "GSE173137", "GSE173138", "GSE173136", "GSE173140", "GSE173139"]
        ),
        rat_ids = expand("../ratgtex/v4/{tissue}/rat_ids.txt", tissue=tissues_separate),
    output:
        rats_tsv = "data/ref/RatGTEx_rats.v4.tsv",
        samples_tsv = "data/ref/RatGTEx_samples.v4.tsv",
        rats_html = "data/ref/rats.v4.html",
        samples_html = "data/ref/samples.v4.html",
    params:
        tissues = tissues_separate
    shell:
        "Rscript scripts/sample_table.R ../ratgtex data {params.tissues}"

rule copy_meta_files:
    """Copy fastq_map and rat_ids files"""
    input:
        "../ratgtex/v4/{tissue}/{ftype}.txt"
    output:
        "data/{ftype}/{ftype}.{tissue}.v4.txt"
    wildcard_constraints:
        ftype = "fastq_map|rat_ids"
    shell:
        "cp {input} {output}"

rule copy_covar_files:
    """Copy covariate files"""
    input:
        "../ratgtex/v4/{tissue}/pheast/intermediate/covar/{modality}.covar.tsv"
    output:
        "data/covar/covar.{tissue}.{modality}.v4.tsv"
    shell:
        "cp {input} {output}"

rule copy_geno_files:
    """Copy genotype files"""
    input:
        "../ratgtex/geno/{geno_dataset}.{ext}"
    output:
        "data/geno/{geno_dataset}.rn8.{ext}"
    wildcard_constraints:
        ext = "vcf.gz|vcf.gz.tbi"
    shell:
        "cp {input} {output}"

#######################################
## Expression and all RNA phenotypes ##
#######################################

rule copy_expr_log2:
    """Copy expression files"""
    input:
        "../ratgtex/v4/{tissue}/{tissue}.expr.log2.bed.gz"
    output:
        "data/expr/expr.log2.{tissue}.v4_rn8.bed.gz"
    shell:
        "cp {input} {output}"

rule copy_expr_tpm:
    """Copy expression files"""
    input:
        "data/phenos/phenos.{tissue}.expression.unnorm.v4_rn8.bed.gz"
    output:
        "data/expr/expr.tpm.{tissue}.v4_rn8.bed.gz"
    shell:
        "cp {input} {output}"

rule median_gene_expression:
    """Calculate median gene expression for web interface"""
    input:
        expr = expand("data/expr/expr.tpm.{tissue}.v4_rn8.bed.gz", tissue=tissues_merged),
        genes = "data/gene.v4.tsv",
    output:
        median = "data/expr/medianGeneExpression.v4.tsv.gz",
        top = "data/expr/topExpressedGene.v4.tsv",
    params:
        tissues = tissues_merged,
    shell:
        "Rscript scripts/medianGeneExpression.R ../ratgtex data {params.tissues}"

rule copy_phenos_unnorm:
    """Copy unnormalized phenotype files"""
    input:
        "../ratgtex/v4/{tissue}/phenos/output/unnorm/{modality}.bed"
    output:
        "data/phenos/phenos.{tissue}.{modality}.unnorm.v4_rn8.bed.gz"
    shell:
        "bgzip -c {input} > {output}"

rule copy_phenos_norm:
    """Copy normalized phenotype files"""
    input:
        "../ratgtex/v4/{tissue}/phenos/output/{modality}.bed.gz"
    output:
        "data/phenos/phenos.{tissue}.{modality}.norm.v4_rn8.bed.gz"
    shell:
        "cp {input} {output}"

#######################
## eQTL and all xQTL ##
#######################

rule eqtl_combined_files:
    """Process and aggregate eQTL results from tensorQTL"""
    input:
        alleles = "../ratgtex/geno/alleles.txt.gz",
        gtf = "../ratgtex/ref/GCF_036323735.1_GRCr8_genomic.chr.gtf",
        afc = expand("../ratgtex/v4/{tissue}/{tissue}.aFC.txt", tissue=tissues_merged),
        cis_indep = expand("../ratgtex/v4/{tissue}/pheast/output/qtl/expression.cis_independent_qtl.txt.gz", tissue=tissues_merged),
        cis_qtl = expand("../ratgtex/v4/{tissue}/pheast/output/qtl/expression.cis_qtl.txt.gz", tissue=tissues_merged),
    output:
        top_assoc = "data/eqtl/top_assoc.v4_rn8.txt",
        eqtls_indep = "data/eqtl/eqtls_indep.v4_rn8.txt",
    params:
        tissues = tissues_merged,
    resources:
        mem_mb = 32000,
    shell:
        "Rscript scripts/eqtl_combined_files.R ../ratgtex data {params.tissues}"

rule xqtl_combined_files:
    """Process and aggregate xQTL results from tensorQTL"""
    input:
        alleles = "../ratgtex/geno/alleles.txt.gz",
        genes = "data/gene.v4.tsv",
        cis_indep = expand("../ratgtex/v4/{tissue}/pheast/output/qtl/{{modality}}.cis_independent_qtl.txt.gz", tissue=tissues_merged),
        cis_qtl = expand("../ratgtex/v4/{tissue}/pheast/output/qtl/{{modality}}.cis_qtl.txt.gz", tissue=tissues_merged),
    output:
        top_assoc = "data/xqtl/top_assoc.{modality}.v4_rn8.txt",
        eqtls_indep = "data/xqtl/xqtls_indep.{modality}.v4_rn8.txt",
    params:
        tissues = tissues_merged,
    resources:
        mem_mb = 32000,
    shell:
        "Rscript scripts/xqtl_combined_files.R ../ratgtex data {wildcards.modality} {params.tissues}"

rule xqtl_signif_files:
    """Process all significant cis-xQTL pair file for one tissue"""
    input:
        signif = "../ratgtex/v4/{tissue}/{tissue}.{modality}.cis_qtl_signif.txt.gz",
        genes = "data/gene.v4.tsv",
    output:
        signif = "data/xqtl/cis_qtl_signif.{tissue}.{modality}.v4_rn8.txt.gz",
    resources:
        mem_mb = 32000,
    shell:
        "Rscript scripts/cis_qtl_signif.R {input.signif} {input.genes} {output.signif}"

rule xqtl_trans_files:
    """Process eQTL results for one tissue"""
    input:
        trans_qtl_pairs = "../ratgtex/v4/{tissue}/{tissue}.{modality}.trans_qtl_pairs.txt.gz",        
    output:
        trans_qtl_pairs = "data/xqtl/trans_qtl_pairs.{tissue}.{modality}.v4_rn8.txt.gz"
    resources:
        mem_mb = 32000,
    shell:
        """
        zcat {input.trans_qtl_pairs} | \
            sed '1s/\\bb\\b/slope/g; 1s/\\bb_se\\b/slope_se/g' | \
            gzip -c \
            > {output.trans_qtl_pairs}
        """

rule all_signif_eqtl:
    """Assemble all significant cis associations in zip archive of per-gene files"""
    input:
        signif = expand("../ratgtex/v4/{tissue}/{tissue}.expression.cis_qtl_signif.txt.gz", tissue=tissues_merged),
    output:
        zip = "data/eqtl/singleTissueEqtl.v4.zip"
    params:
        tissues = tissues_merged,
    resources:
        mem_mb = lambda w, attempt: 64000 * 2**(attempt-1),
    retries: 2
    shell:
        "python3 scripts/singleTissueEqtl.py ../ratgtex v4 data {params.tissues}"

rule cis_pvals:
    """Save all cis p-values in zip archive of per-gene files"""
    input:
        pvals = "../ratgtex/v4/{tissue}/{tissue}.cis_qtl_all_pvals.txt.gz",
    output:
        zip = "data/cis_pvals/{tissue}.v4.zip"
    shell:
        "python3 scripts/all_cis_pvals.py ../ratgtex v4 data {wildcards.tissue}"

# rule sqtl_combined_files:
#     """Process sQTL results from tensorQTL for web interface"""
#     input:
#         alleles = "../ratgtex/geno/alleles.txt.gz",
#         gtf = "../ratgtex/ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf",
#         cis_qtl = expand("../ratgtex/v3/{tissue}/splice/{tissue}_splice.cis_qtl.txt.gz", tissue=tissues_merged),
#         cis_indep = expand("../ratgtex/v3/{tissue}/splice/{tissue}_splice.cis_independent_qtl.txt.gz", tissue=tissues_merged),
#     output:
#         top_assoc = "data/splice/top_assoc_splice.v3_rn7.txt",
#         sqtls_indep = "data/splice/sqtls_indep.v3_rn7.txt",
#     params:
#         tissues = tissues_merged,
#     resources:
#         mem_mb = 32000,
#     shell:
#         "Rscript scripts/sqtl_combined_files.R ../ratgtex data {params.tissues}"

# rule sqtl_tissue_files:
#     """Process sQTL results for one tissue"""
#     input:
#         gtf = "../ratgtex/ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf",
#         cis_qtl_signif = "../ratgtex/v3/{tissue}/splice/{tissue}_splice.cis_qtl_signif.txt.gz",
#         trans_qtl_pairs = "../ratgtex/v3/{tissue}/splice/{tissue}_splice.trans_qtl_pairs.txt.gz",
#     output:
#         cis_qtl_signif = "data/splice/splice.cis_qtl_signif.{tissue}.v3_rn7.txt.gz",
#         trans_qtl_pairs = "data/splice/splice.trans_qtl_pairs.{tissue}.v3_rn7.txt.gz"
#     resources:
#         mem_mb = 32000,
#     shell:
#         "Rscript scripts/sqtl_tissue_files.R ../ratgtex data {wildcards.tissue}"


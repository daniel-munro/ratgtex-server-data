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
        "data/v4/tissue_info.v4.tsv",
        ##
        "data/v4/gene.v4.tsv",
        "data/v4/autocomplete.v4.json",
        ##
        "data/v4/exon.v4.tsv",
        ##
        "data/v4/ref/RatGTEx_rats.v4.tsv",
        "data/v4/ref/RatGTEx_samples.v4.tsv",
        "data/v4/ref/rats.v4.html",
        "data/v4/ref/samples.v4.html",
        ##
        expand("data/v4/covar/covar.{tissue}.{modality}.v4.tsv", tissue=tissues_merged, modality=modalities_cross),
        expand("data/v4/fastq_map/fastq_map.{tissue}.v4.txt", tissue=tissues_merged),
        expand("data/v4/rat_ids/rat_ids.{tissue}.v4.txt", tissue=tissues_merged),
        ##
        expand("data/v4/geno/{geno_dataset}.rn8.{ext}", geno_dataset=geno_datasets, ext=["vcf.gz", "vcf.gz.tbi"]),
        ##
        expand("data/v4/expr/expr.{units}.{tissue}.v4_rn8.bed.gz", units=["log2", "tpm"], tissue=tissues_merged),
        ##
        "data/v4/expr/medianGeneExpression.v4.tsv.gz",
        "data/v4/expr/topExpressedGene.v4.tsv",
        ##
        expand("data/v4/phenos/phenos.{tissue}.{modality}.{norm}.v4_rn8.bed.gz", tissue=tissues_merged, modality=modalities, norm=["unnorm", "norm"]),
        ##
        "data/v4/eqtl/top_assoc.v4_rn8.tsv",
        "data/v4/eqtl/eqtls_indep.v4_rn8.tsv",
        ##
        "data/v4/eqtl/all_gene_esnps.v4.zip",
        ##
        expand("data/v4/xqtl/top_assoc.{modality}.v4_rn8.tsv", modality=modalities + ["cross_modality"]),
        expand("data/v4/xqtl/xqtls_indep.{modality}.v4_rn8.tsv", modality=modalities + ["cross_modality"]),
        ##
        expand("data/v4/xqtl/cis_qtl_signif.{tissue}.{modality}.v4_rn8.txt.gz", tissue=tissues_merged, modality=modalities),
        expand("data/v4/xqtl/trans_qtl_pairs.{tissue}.expression.v4_rn8.txt.gz", tissue=tissues_merged),
        ##
        expand("data/v4/cis_pvals/{tissue}.v4.zip", tissue=tissues_merged),

############################
## Metadata and reference ##
############################

rule tissue_info:
    """Process tissue info table"""
    input:
        info = "../ratgtex/tissue_info.tsv",
        expr = expand("../ratgtex/v4/{tissue}/phenos/output/expression.bed.gz", tissue=tissues_merged),
        eqtl = "data/v4/eqtl/top_assoc.v4_rn8.tsv",
    output:
        info = "data/v4/tissue_info.v4.tsv"
    params:
        tissues = tissues_merged
    shell:
        "Rscript scripts/tissue_info.R ../ratgtex data/v4 {params.tissues}"

rule gene_info:
    """Process gene info table"""
    input:
        gtf = "../ratgtex/ref/GCF_036323735.1_GRCr8_genomic.chr.gtf",
        signif = expand("../ratgtex/v4/{tissue}/{tissue}.expression.cis_qtl_signif.txt.gz", tissue=tissues_merged),
        expr = expand("../ratgtex/v4/{tissue}/phenos/output/expression.bed.gz", tissue=tissues_merged),
        eqtl = "data/v4/eqtl/eqtls_indep.v4_rn8.tsv",
        expr_assoc = "data/v4/eqtl/top_assoc.v4_rn8.tsv",
    output:
        gene = "data/v4/gene.v4.tsv",
        autocomplete = "data/v4/autocomplete.v4.json",
    params:
        tissues = tissues_merged
    shell:
        "Rscript scripts/gene.R ../ratgtex data/v4 {params.tissues}"

rule exon_table:
    """Prepare exon annotations for server API"""
    input:
        gtf = "../ratgtex/ref/GCF_036323735.1_GRCr8_genomic.chr.genes.gtf"
    output:
        exon = "data/v4/exon.v4.tsv"
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
        rats_tsv = "data/v4/ref/RatGTEx_rats.v4.tsv",
        samples_tsv = "data/v4/ref/RatGTEx_samples.v4.tsv",
        rats_html = "data/v4/ref/rats.v4.html",
        samples_html = "data/v4/ref/samples.v4.html",
    params:
        tissues = tissues_separate
    shell:
        "Rscript scripts/sample_table.R ../ratgtex data/v4 {params.tissues}"

rule copy_meta_files:
    """Copy fastq_map and rat_ids files"""
    input:
        "../ratgtex/v4/{tissue}/{ftype}.txt"
    output:
        "data/v4/{ftype}/{ftype}.{tissue}.v4.txt"
    wildcard_constraints:
        ftype = "fastq_map|rat_ids"
    shell:
        "cp {input} {output}"

rule copy_covar_files:
    """Copy covariate files"""
    input:
        "../ratgtex/v4/{tissue}/pheast/intermediate/covar/{modality}.covar.tsv"
    output:
        "data/v4/covar/covar.{tissue}.{modality}.v4.tsv"
    shell:
        "cp {input} {output}"

rule copy_geno_files:
    """Copy genotype files"""
    input:
        "../ratgtex/geno/{geno_dataset}.{ext}"
    output:
        "data/v4/geno/{geno_dataset}.rn8.{ext}"
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
        "data/v4/expr/expr.log2.{tissue}.v4_rn8.bed.gz"
    shell:
        "cp {input} {output}"

rule copy_expr_tpm:
    """Copy expression files"""
    input:
        "data/v4/phenos/phenos.{tissue}.expression.unnorm.v4_rn8.bed.gz"
    output:
        "data/v4/expr/expr.tpm.{tissue}.v4_rn8.bed.gz"
    shell:
        "cp {input} {output}"

rule median_gene_expression:
    """Calculate median gene expression for web interface"""
    input:
        expr = expand("data/v4/expr/expr.tpm.{tissue}.v4_rn8.bed.gz", tissue=tissues_merged),
        genes = "data/v4/gene.v4.tsv",
    output:
        median = "data/v4/expr/medianGeneExpression.v4.tsv.gz",
        top = "data/v4/expr/topExpressedGene.v4.tsv",
    params:
        tissues = tissues_merged,
    shell:
        "Rscript scripts/medianGeneExpression.R ../ratgtex data/v4 {params.tissues}"

rule copy_phenos_unnorm:
    """Copy unnormalized phenotype files"""
    input:
        "../ratgtex/v4/{tissue}/phenos/output/unnorm/{modality}.bed"
    output:
        "data/v4/phenos/phenos.{tissue}.{modality}.unnorm.v4_rn8.bed.gz"
    shell:
        "bgzip -c {input} > {output}"

rule copy_phenos_norm:
    """Copy normalized phenotype files"""
    input:
        "../ratgtex/v4/{tissue}/phenos/output/{modality}.bed.gz"
    output:
        "data/v4/phenos/phenos.{tissue}.{modality}.norm.v4_rn8.bed.gz"
    shell:
        "cp {input} {output}"

#######################
## eQTL and all xQTL ##
#######################

rule eqtl_combined_files:
    """Process and aggregate eQTL results from tensorQTL"""
    input:
        alleles = "../ratgtex/geno/alleles.tsv.gz",
        gtf = "../ratgtex/ref/GCF_036323735.1_GRCr8_genomic.chr.gtf",
        afc = expand("../ratgtex/v4/{tissue}/{tissue}.aFC.tsv", tissue=tissues_merged),
        cis_indep = expand("../ratgtex/v4/{tissue}/pheast/output/qtl/expression.cis_independent_qtl.txt.gz", tissue=tissues_merged),
        cis_qtl = expand("../ratgtex/v4/{tissue}/pheast/output/qtl/expression.cis_qtl.txt.gz", tissue=tissues_merged),
    output:
        top_assoc = "data/v4/eqtl/top_assoc.v4_rn8.tsv",
        eqtls_indep = "data/v4/eqtl/eqtls_indep.v4_rn8.tsv",
    params:
        tissues = tissues_merged,
    resources:
        mem_mb = 32000,
    shell:
        "Rscript scripts/eqtl_combined_files.R ../ratgtex data/v4 {params.tissues}"

rule xqtl_combined_files:
    """Process and aggregate xQTL results from tensorQTL"""
    input:
        alleles = "../ratgtex/geno/alleles.tsv.gz",
        genes = "data/v4/gene.v4.tsv",
        cis_indep = expand("../ratgtex/v4/{tissue}/pheast/output/qtl/{{modality}}.cis_independent_qtl.txt.gz", tissue=tissues_merged),
        cis_qtl = expand("../ratgtex/v4/{tissue}/pheast/output/qtl/{{modality}}.cis_qtl.txt.gz", tissue=tissues_merged),
    output:
        top_assoc = "data/v4/xqtl/top_assoc.{modality}.v4_rn8.tsv",
        eqtls_indep = "data/v4/xqtl/xqtls_indep.{modality}.v4_rn8.tsv",
    params:
        tissues = tissues_merged,
    resources:
        mem_mb = 32000,
    shell:
        "Rscript scripts/xqtl_combined_files.R ../ratgtex data/v4 {wildcards.modality} {params.tissues}"

rule xqtl_signif_files:
    """Process all significant cis-xQTL pair file for one tissue"""
    input:
        signif = "../ratgtex/v4/{tissue}/{tissue}.{modality}.cis_qtl_signif.txt.gz",
        genes = "data/v4/gene.v4.tsv",
    output:
        signif = "data/v4/xqtl/cis_qtl_signif.{tissue}.{modality}.v4_rn8.txt.gz",
    resources:
        mem_mb = 32000,
    shell:
        "Rscript scripts/cis_qtl_signif.R {input.signif} {input.genes} {output.signif}"

rule xqtl_trans_files:
    """Process eQTL results for one tissue"""
    input:
        trans_qtl_pairs = "../ratgtex/v4/{tissue}/{tissue}.{modality}.trans_qtl_pairs.txt.gz",        
    output:
        trans_qtl_pairs = "data/v4/xqtl/trans_qtl_pairs.{tissue}.{modality}.v4_rn8.txt.gz"
    resources:
        mem_mb = 32000,
    shell:
        """
        zcat {input.trans_qtl_pairs} | \
            sed '1s/\\bb\\b/slope/g; 1s/\\bb_se\\b/slope_se/g' | \
            gzip -c \
            > {output.trans_qtl_pairs}
        """

rule all_gene_esnps:
    """Assemble all significant cis associations in zip archive of per-gene files"""
    input:
        signif = expand("../ratgtex/v4/{tissue}/{tissue}.expression.cis_qtl_signif.txt.gz", tissue=tissues_merged),
    output:
        zip = "data/v4/eqtl/all_gene_esnps.v4.zip"
    params:
        tissues = tissues_merged,
    resources:
        mem_mb = 128000,
    retries: 2
    shell:
        "python3 scripts/all_gene_esnps.py ../ratgtex v4 data/v4 {params.tissues}"

rule cis_pvals:
    """Save all expression cis p-values in zip archive of per-gene files"""
    input:
        pvals = "../ratgtex/v4/{tissue}/{tissue}.expression.cis_qtl_all_pvals.tsv.gz",
    output:
        zip = "data/v4/cis_pvals/{tissue}.v4.zip"
    shell:
        "python3 scripts/all_cis_pvals.py ../ratgtex v4 data/v4 {wildcards.tissue}"

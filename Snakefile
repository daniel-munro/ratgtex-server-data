from pathlib import Path

configfile: "../ratgtex/config.yaml"

tissues_separate = [tissue for tissue in config["tissues"] if (indir / version / tissue).exists()]
tissues_all = tissues_separate + config["merged_tissues"].keys()
tissues_merged = [tissue for tissue in tissues_all if tissue not in sum(config["merged_tissues"].values(), [])]
print(f"Separate tissues: {tissues_separate}")
print(f"Merged tissues: {tissues_merged}")
geno_datasets = list(set(config["tissues"][tissue]["geno_dataset"] for tissue in tissues_separate))

localrules:
    exon_table,
    sample_table,
    copy_meta_files,
    copy_geno_files,
    copy_expr_files,
    median_gene_expression,
    copy_splice_files,

rule all:
    input:
        "data/tissueInfo.v3.txt",
        ##
        "data/gene.v3_rn7.txt",
        "data/autocomplete.v3_rn7.json",
        ##
        "data/exon.v3_rn7.txt",
        ##
        "data/ref/RatGTEx_rats.v3.tsv",
        "data/ref/RatGTEx_samples.v3.tsv",
        "data/ref/rats.v3.html",
        "data/ref/samples.v3.html",
        ##
        expand("data/covar/covar.{tissue}.v3.txt", tissue=tissues_all),
        expand("data/fastq_map/fastq_map.{tissue}.v3.txt", tissue=tissues_all),
        expand("data/rat_ids/rat_ids.{tissue}.v3.txt", tissue=tissues_all),
        ##
        expand("data/geno/{geno_dataset}.rn7.{ext}", geno_dataset=geno_datasets, ext=["vcf.gz", "vcf.gz.tbi"]),
        ##
        expand("data/expr/expr.{units}.{tissue}.v3_rn7.bed.gz", units=["iqn", "log2", "tpm"], tissue=tissues_all),
        ##
        "data/medianGeneExpression.v3_rn7.txt.gz",
        "data/topExpressedGene.v3_rn7.txt",
        ##
        expand("data/splice/leafcutter.{tissue}.v3_rn7.bed.gz", tissue=tissues_all),
        expand("data/splice/covar_splice.{tissue}.v3.txt", tissue=tissues_all),
        ##
        "data/eqtl/top_assoc.v3_rn7.txt",
        "data/eqtl/eqtls_indep.v3_rn7.txt",
        ##
        expand("data/eqtl/cis_qtl_signif.{tissue}.v3_rn7.txt.gz", tissue=tissues_all),
        expand("data/eqtl/trans_qtl_pairs.{tissue}.v3_rn7.txt.gz", tissue=tissues_all),
        ##
        "data/singleTissueEqtl.v3.zip",
        ##
        expand("data/cis_pvals/{tissue}.v3.zip", tissue=tissues_merged),
        ##
        "data/splice/top_assoc_splice.v3_rn7.txt",
        "data/splice/sqtls_indep.v3_rn7.txt",
        ##
        expand("data/splice/splice.cis_qtl_signif.{tissue}.v3_rn7.txt.gz", tissue=tissues_all),
        expand("data/splice/splice.trans_qtl_pairs.{tissue}.v3_rn7.txt.gz", tissue=tissues_all),

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
        info = "../ratgtex/tissue_info.txt",
        expr = expand("../ratgtex/v3/{tissue}/{tissue}.expr.tpm.bed.gz", tissue=tissues_merged),
        eqtl = "data/eqtl/top_assoc.v3_rn7.txt",
    output:
        info = "data/tissueInfo.v3.txt"
    params:
        tissues = tissues_merged
    shell:
        "Rscript scripts/tissueInfo.R ../ratgtex data {params.tissues}"

rule gene_info:
    """Process gene info table"""
    input:
        gtf = "../ratgtex/ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf",
        genes = "../ratgtex/ref/GENES_RAT.txt",
        signif = expand("../ratgtex/v3/{tissue}/{tissue}.cis_qtl_signif.txt.gz", tissue=tissues_merged),
        expr = expand("../ratgtex/v3/{tissue}/{tissue}.expr.tpm.bed.gz", tissue=tissues_merged),
        leafcutter = expand("../ratgtex/v3/{tissue}/splice/{tissue}.leafcutter.phenotype_groups.txt", tissue=tissues_merged),
        eqtl = "data/eqtl/eqtls_indep.v3_rn7.txt",
        expr_assoc = "data/eqtl/top_assoc.v3_rn7.txt",
        sqtl = "data/splice/sqtls_indep.v3_rn7.txt",
        splice_assoc = "data/splice/top_assoc_splice.v3_rn7.txt",
    output:
        gene = "data/gene.v3_rn7.txt",
        autocomplete = "data/autocomplete.v3_rn7.json",
    params:
        tissues = tissues_merged
    shell:
        "Rscript scripts/gene.R ../ratgtex data {params.tissues}"

rule exon_table:
    """Prepare exon annotations for server API"""
    input:
        gtf = "../ratgtex/ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf"
    output:
        exon = "data/exon.v3_rn7.txt"
    shell:
        "python3 scripts/exon.py {input.gtf} {output.exon}"

rule sample_table:
    """Generate sample and rat tables for web interface"""
    input:
        geno_log = "../ratgtex/geno/genotyping_log.csv",
        geo = expand("../ratgtex/samples/{accession}_series_matrix.txt.gz", accession=["GSE201236", "GSE173137", "GSE173138", "GSE173136", "GSE173140", "GSE173139"]),
        rat_ids = expand("../ratgtex/v3/{tissue}/rat_ids.txt", tissue=tissues_separate),
    output:
        rats_tsv = "data/ref/RatGTEx_rats.v3.tsv",
        samples_tsv = "data/ref/RatGTEx_samples.v3.tsv",
        rats_html = "data/ref/rats.v3.html",
        samples_html = "data/ref/samples.v3.html",
    params:
        tissues = tissues_separate
    shell:
        "Rscript scripts/sample_table.R ../ratgtex data {params.tissues}"

rule copy_meta_files:
    """Copy covariate, fastq_map, and rat_ids files"""
    input:
        "../ratgtex/v3/{tissue}/{ftype}.txt"
    output:
        "data/{ftype}/{ftype}.{tissue}.v3.txt"
    wildcard_constraints:
        ftype = "covar|fastq_map|rat_ids"
    shell:
        "cp {input} {output}"

rule copy_geno_files:
    """Copy genotype files"""
    input:
        "../ratgtex/geno/{geno_dataset}.{ext}"
    output:
        "data/geno/{geno_dataset}.rn7.{ext}"
    wildcard_constraints:
        ext = "vcf.gz|vcf.gz.tbi"
    shell:
        "cp {input} {output}"

#############################
## Expression and splicing ##
#############################

rule copy_expr_files:
    """Copy expression files"""
    input:
        "../ratgtex/v3/{tissue}/{tissue}.expr.{units}.bed.gz"
    output:
        "data/expr/expr.{units}.{tissue}.v3_rn7.bed.gz"
    shell:
        "cp {input} {output}"

rule median_gene_expression:
    """Calculate median gene expression for web interface"""
    input:
        expr = expand("../ratgtex/v3/{tissue}/{tissue}.expr.tpm.bed.gz", tissue=tissues_merged),
        genes = "data/gene.v3_rn7.txt",
    output:
        median = "data/medianGeneExpression.v3_rn7.txt.gz",
        top = "data/topExpressedGene.v3_rn7.txt",
    params:
        tissues = tissues_merged,
    shell:
        "Rscript scripts/medianGeneExpression.R ../ratgtex data {params.tissues}"

rule copy_splice_files:
    """Copy leafCutter and splice covariate files"""
    input:
        leafcutter = "../ratgtex/v3/{tissue}/splice/{tissue}.leafcutter.bed.gz",
        covar = "../ratgtex/v3/{tissue}/splice/{tissue}.covar_splice.txt"
    output:
        leafcutter = "data/splice/leafcutter.{tissue}.v3_rn7.bed.gz",
        covar = "data/splice/covar_splice.{tissue}.v3.txt"
    shell:
        """
        cp {input.leafcutter} {output.leafcutter}
        cp {input.covar} {output.covar}
        """

###################
## eQTL and sQTL ##
###################

rule eqtl_combined_files:
    """Process and aggregate eQTL results from tensorQTL"""
    input:
        alleles = "../ratgtex/geno/alleles.txt.gz",
        gtf = "../ratgtex/ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf",
        afc = expand("../ratgtex/v3/{tissue}/{tissue}.aFC.txt", tissue=tissues_merged),
        cis_indep = expand("../ratgtex/v3/{tissue}/{tissue}.cis_independent_qtl.txt.gz", tissue=tissues_merged),
        cis_qtl = expand("../ratgtex/v3/{tissue}/{tissue}.cis_qtl.txt.gz", tissue=tissues_merged),
    output:
        top_assoc = "data/eqtl/top_assoc.v3_rn7.txt",
        eqtls_indep = "data/eqtl/eqtls_indep.v3_rn7.txt",
    params:
        tissues = tissues_merged,
    resources:
        mem_mb = 32000,
    shell:
        "Rscript scripts/eqtl_combined_files.R ../ratgtex data {params.tissues}"

rule eqtl_tissue_files:
    """Process eQTL results for one tissue"""
    input:
        gtf = "../ratgtex/ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf",
        cis_qtl_signif = "../ratgtex/v3/{tissue}/{tissue}.cis_qtl_signif.txt.gz",
        trans_qtl_pairs = "../ratgtex/v3/{tissue}/{tissue}.trans_qtl_pairs.txt.gz",        
    output:
        cis_qtl_signif = "data/eqtl/cis_qtl_signif.{tissue}.v3_rn7.txt.gz",
        trans_qtl_pairs = "data/eqtl/trans_qtl_pairs.{tissue}.v3_rn7.txt.gz"
    resources:
        mem_mb = 32000,
    shell:
        "Rscript scripts/eqtl_tissue_files.R ../ratgtex data {wildcards.tissue}"

rule all_signif_eqtl:
    """Assemble all significant cis associations in zip archive of per-gene files"""
    input:
        signif = expand("../ratgtex/v3/{tissue}/{tissue}.cis_qtl_signif.txt.gz", tissue=tissues_merged),
    output:
        zip = "data/singleTissueEqtl.v3.zip"
    params:
        tissues = tissues_merged,
    resources:
        mem_mb = 32000,
    shell:
        "python3 scripts/singleTissueEqtl.py ../ratgtex v3 data {params.tissues}"

rule cis_pvals:
    """Save all cis p-values in zip archive of per-gene files"""
    input:
        pvals = "../ratgtex/v3/{tissue}/{tissue}.cis_qtl_all_pvals.txt.gz",
    output:
        zip = "data/cis_pvals/{tissue}.v3.zip"
    shell:
        "python3 scripts/all_cis_pvals.py ../ratgtex v3 data {wildcards.tissue}"

rule sqtl_combined_files:
    """Process sQTL results from tensorQTL for web interface"""
    input:
        alleles = "../ratgtex/geno/alleles.txt.gz",
        gtf = "../ratgtex/ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf",
        cis_qtl = expand("../ratgtex/v3/{tissue}/splice/{tissue}_splice.cis_qtl.txt.gz", tissue=tissues_merged),
        cis_indep = expand("../ratgtex/v3/{tissue}/splice/{tissue}_splice.cis_independent_qtl.txt.gz", tissue=tissues_merged),
    output:
        top_assoc = "data/splice/top_assoc_splice.v3_rn7.txt",
        sqtls_indep = "data/splice/sqtls_indep.v3_rn7.txt",
    params:
        tissues = tissues_merged,
    resources:
        mem_mb = 32000,
    shell:
        "Rscript scripts/sqtl_combined_files.R ../ratgtex data {params.tissues}"

rule sqtl_tissue_files:
    """Process sQTL results for one tissue"""
    input:
        gtf = "../ratgtex/ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf",
        cis_qtl_signif = "../ratgtex/v3/{tissue}/splice/{tissue}_splice.cis_qtl_signif.txt.gz",
        trans_qtl_pairs = "../ratgtex/v3/{tissue}/splice/{tissue}_splice.trans_qtl_pairs.txt.gz",
    output:
        cis_qtl_signif = "data/splice/splice.cis_qtl_signif.{tissue}.v3_rn7.txt.gz",
        trans_qtl_pairs = "data/splice/splice.trans_qtl_pairs.{tissue}.v3_rn7.txt.gz"
    resources:
        mem_mb = 32000,
    shell:
        "Rscript scripts/sqtl_tissue_files.R ../ratgtex data {wildcards.tissue}"


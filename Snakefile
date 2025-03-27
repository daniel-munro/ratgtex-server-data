from pathlib import Path

configfile: "../ratgtex/config.yaml"

tissues = [tissue for tissue in config["tissues"] if (indir / version / tissue).exists()]
geno_datasets = list(set(config["tissues"][tissue]["geno_dataset"] for tissue in tissues))

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
        "data/tissueInfo.v3_rn7.txt",
        ##
        "data/gene.v3_rn7.txt",
        "data/autocomplete.v3_rn7.json",
        ##
        "data/exon.v3_rn7.txt",
        ##
        "data/ref/RatGTEx_rats.v3_rn7.tsv",
        "data/ref/RatGTEx_samples.v3_rn7.tsv",
        "data/ref/rats.v3_rn7.html",
        "data/ref/samples.v3_rn7.html",
        ##
        expand("data/covar/covar.{tissue}.v3_rn7.txt", tissue=tissues),
        expand("data/fastq_map/fastq_map.{tissue}.v3_rn7.txt", tissue=tissues),
        expand("data/rat_ids/rat_ids.{tissue}.v3_rn7.txt", tissue=tissues),
        ##
        expand("data/geno/{geno_dataset}.rn7.{ext}", geno_dataset=geno_datasets, ext=["vcf.gz", "vcf.gz.tbi"]),
        ##
        expand("data/expr/expr.{units}.{tissue}.v3_rn7.bed.gz", units=["iqn", "log2", "tpm"], tissue=tissues),
        ##
        "data/medianGeneExpression.v3_rn7.txt.gz",
        "data/topExpressedGene.v3_rn7.txt",
        ##
        expand("data/splice/leafcutter.{tissue}.v3_rn7.bed.gz", tissue=tissues),
        expand("data/splice/covar_splice.{tissue}.v3_rn7.txt", tissue=tissues),
        ##
        "data/eqtl/top_assoc.v3_rn7.txt",
        "data/eqtl/eqtls_indep.v3_rn7.txt",
        expand("data/eqtl/cis_qtl_signif.{tissue}.v3_rn7.txt.gz", tissue=tissues),
        expand("data/eqtl/trans_qtl_pairs.{tissue}.v3_rn7.txt.gz", tissue=tissues),
        # 
        "data/singleTissueEqtl.v3.zip",
        ##
        expand("data/cis_pvals/{tissue}.v3.zip", tissue=tissues),
        ##
        "data/splice/top_assoc_splice.v3_rn7.txt",
        "data/splice/sqtls_indep.v3_rn7.txt",
        expand("data/splice/splice.cis_qtl_signif.{tissue}.v3_rn7.txt.gz", tissue=tissues),
        expand("data/splice/splice.trans_qtl_pairs.{tissue}.v3_rn7.txt.gz", tissue=tissues),

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
        expr = expand("../ratgtex/v3/{tissue}/{tissue}.expr.tpm.bed.gz", tissue=tissues),
        eqtl = "data/eqtl/top_assoc.v3_rn7.txt",
    output:
        info = "data/tissueInfo.v3_rn7.txt"
    params:
        tissues = tissues
    shell:
        "Rscript scripts/tissueInfo.R ../ratgtex data {tissues}"

rule gene_info:
    """Process gene info table"""
    input:
        gtf = "../ratgtex/ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf",
        genes = "../ratgtex/ref/GENES_RAT.txt",
        signif = expand("../ratgtex/v3/{tissue}/{tissue}.cis_qtl_signif.txt.gz", tissue=tissues),
        expr = expand("../ratgtex/v3/{tissue}/{tissue}.expr.tpm.bed.gz", tissue=tissues),
        leafcutter = expand("../ratgtex/v3/{tissue}/splice/{tissue}.leafcutter.phenotype_groups.txt", tissue=tissues),
        eqtl = "data/eqtl/eqtls_indep.v3_rn7.txt",
        expr_assoc = "data/eqtl/top_assoc.v3_rn7.txt",
        sqtl = "data/splice/sqtls_indep.v3_rn7.txt",
        splice_assoc = "data/splice/top_assoc_splice.v3_rn7.txt",
    output:
        gene = "data/gene.v3_rn7.txt",
        autocomplete = "data/autocomplete.v3_rn7.json",
    params:
        tissues = tissues
    shell:
        "Rscript scripts/gene.R ../ratgtex data {tissues}"

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
        rat_ids = expand("../ratgtex/v3/{tissue}/rat_ids.txt", tissue=tissues),
    output:
        rats_tsv = "data/ref/RatGTEx_rats.v3_rn7.tsv",
        samples_tsv = "data/ref/RatGTEx_samples.v3_rn7.tsv",
        rats_html = "data/ref/rats.v3_rn7.html",
        samples_html = "data/ref/samples.v3_rn7.html",
    params:
        tissues = tissues
    shell:
        "Rscript scripts/sample_table.R ../ratgtex data {tissues}"

rule copy_meta_files:
    """Copy covariate, fastq_map, and rat_ids files"""
    input:
        "../ratgtex/v3/{tissue}/{ftype}.txt"
    output:
        "data/{ftype}/{ftype}.{tissue}.v3_rn7.txt"
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
        expr = expand("../ratgtex/v3/{tissue}/{tissue}.expr.tpm.bed.gz", tissue=tissues),
        genes = "data/gene.v3_rn7.txt",
    output:
        median = "data/medianGeneExpression.v3_rn7.txt.gz",
        top = "data/topExpressedGene.v3_rn7.txt",
    params:
        tissues = tissues,
    shell:
        "Rscript scripts/medianGeneExpression.R ../ratgtex data {tissues}"

rule copy_splice_files:
    """Copy leafCutter and splice covariate files"""
    input:
        leafcutter = "../ratgtex/v3/{tissue}/splice/{tissue}.leafcutter.bed.gz",
        covar = "../ratgtex/v3/{tissue}/splice/{tissue}.covar_splice.txt"
    output:
        leafcutter = "data/splice/leafcutter.{tissue}.v3_rn7.bed.gz",
        covar = "data/splice/covar_splice.{tissue}.v3_rn7.txt"
    shell:
        """
        cp {input.leafcutter} {output.leafcutter}
        cp {input.covar} {output.covar}
        """

###################
## eQTL and sQTL ##
###################

rule eqtl_files:
    """Process eQTL results from tensorQTL for web interface"""
    input:
        alleles = "../ratgtex/geno/alleles.txt.gz",
        gtf = "../ratgtex/ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf",
        afc = expand("../ratgtex/v3/{tissue}/{tissue}.aFC.txt", tissue=tissues),
        cis_indep = expand("../ratgtex/v3/{tissue}/{tissue}.cis_independent_qtl.txt.gz", tissue=tissues),
        cis_qtl = expand("../ratgtex/v3/{tissue}/{tissue}.cis_qtl.txt.gz", tissue=tissues),
        cis_qtl_signif = expand("../ratgtex/v3/{tissue}/{tissue}.cis_qtl_signif.txt.gz", tissue=tissues),
        trans_qtl_pairs = expand("../ratgtex/v3/{tissue}/{tissue}.trans_qtl_pairs.txt.gz", tissue=tissues),        
    output:
        top_assoc = "data/eqtl/top_assoc.v3_rn7.txt",
        eqtls_indep = "data/eqtl/eqtls_indep.v3_rn7.txt",
        cis_qtl_signif = expand("data/eqtl/cis_qtl_signif.{tissue}.v3_rn7.txt.gz", tissue=tissues),
        trans_qtl_pairs = expand("data/eqtl/trans_qtl_pairs.{tissue}.v3_rn7.txt.gz", tissue=tissues)
    params:
        tissues = tissues,
    shell:
        "Rscript scripts/eqtl_files.R ../ratgtex data {tissues}"

rule single_tissue_eqtl:
    """Assemble all significant cis associations in zip archive of per-gene files"""
    input:
        signif = expand("../ratgtex/v3/{tissue}/{tissue}.cis_qtl_signif.txt.gz", tissue=tissues),
    output:
        zip = "data/singleTissueEqtl.v3.zip"
    params:
        tissues = tissues,
    shell:
        "python3 scripts/singleTissueEqtl.py ../ratgtex v3 data {tissues}"

rule cis_pvals:
    """Save all cis p-values in zip archive of per-gene files"""
    input:
        pvals = expand("../ratgtex/v3/{tissue}/{tissue}.cis_qtl_all_pvals.txt.gz", tissue=tissues),
    output:
        zip = "data/cis_pvals/{tissue}.v3.zip"
    shell:
        "python3 scripts/all_cis_pvals.py ../ratgtex v3 data {wildcards.tissue}"

rule sqtl_files:
    """Process sQTL results from tensorQTL for web interface"""
    input:
        alleles = "../ratgtex/geno/alleles.txt.gz",
        gtf = "../ratgtex/ref/GCF_015227675.2_mRatBN7.2_genomic.chr.genes.gtf",
        cis_qtl = expand("../ratgtex/v3/{tissue}/splice/{tissue}_splice.cis_qtl.txt.gz", tissue=tissues),
        cis_indep = expand("../ratgtex/v3/{tissue}/splice/{tissue}_splice.cis_independent_qtl.txt.gz", tissue=tissues),
        cis_qtl_signif = expand("../ratgtex/v3/{tissue}/splice/{tissue}_splice.cis_qtl_signif.txt.gz", tissue=tissues),
        trans_qtl_pairs = expand("../ratgtex/v3/{tissue}/splice/{tissue}_splice.trans_qtl_pairs.txt.gz", tissue=tissues),
    output:
        top_assoc = "data/splice/top_assoc_splice.v3_rn7.txt",
        sqtls_indep = "data/splice/sqtls_indep.v3_rn7.txt",
        cis_qtl_signif = expand("data/splice/splice.cis_qtl_signif.{tissue}.v3_rn7.txt.gz", tissue=tissues),
        trans_qtl_pairs = expand("data/splice/splice.trans_qtl_pairs.{tissue}.v3_rn7.txt.gz", tissue=tissues)
    params:
        tissues = tissues,
    shell:
        "Rscript scripts/sqtl_files.R ../ratgtex data {tissues}"


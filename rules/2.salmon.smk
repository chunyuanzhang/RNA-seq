# ============================================================
# salmon quant (mapping mode, selective alignment)
#   --gcBias --seqBias --validateMappings: 现代标准用法
#   --libType A: 自动判断文库类型
#   索引由 build_salmon_index.sh 预先建好, 此处直接引用目录
# ============================================================

rule salmon_quant:
    input:
        cleandata_1 = "result/CleanData/{sample}.clean.1.fq.gz",
        cleandata_2 = "result/CleanData/{sample}.clean.2.fq.gz"
    output:
        quant =  "result/Salmon/{sample}/quant.sf"
    params:
        idxdir = lambda wildcards:  referenceDir + metainfo_dict_Genome[wildcards.sample] + "/SALMONIndex",
        outdir = "result/Salmon/{sample}",
    threads:
        32
    log:
        "logs/salmon_quant/{sample}.log"
    shell:
        """
        {salmon}
        salmon quant \
            -i {params.idxdir} \
            --libType A \
            -1 {input.cleandata_1} \
            -2 {input.cleandata_2} \
            --gcBias --seqBias \
            -p {threads} \
            -o {params.outdir} \
            2>{log}
        """
 

# salmon 使用的参考基因组和转录组序列一起创建的索引，其中转录组部分用的ID，
# 包含rna等字样，因此tx2gene文件需要单独创建
rule make_tx2tene:
    input:
        gff = lambda wildcards: get_genome_gff(wildcards.species) 
    output:
        tx2gene = "result/Salmon/{species}_trans2symbol.tsv"
    params:
        Species = "{species}"
    shell:
        """
        {R453}
        Rscript scripts/tx2gene.R --gfffile {input.gff}  --Species {params.Species}
        """


if design == "Species":
    # 直接挑选同源基因，任何参数都不修改
    rule DESeq_between_species:
        input:
            quant =  expand("result/Salmon/{sample}/quant.sf", sample = samples),
            tx2gene = expand("result/Salmon/{species}_trans2symbol.tsv", species = specieses),
            Orthologgenes = Orthologgenes
        output:
            all_genes = "result/Salmon/{pairname}.all.tsv",
            diffexp_genes = "result/Salmon/{pairname}.diffexp.tsv"
        params:
            infotable = infotable,
            design = design,
            lfc = lfc,
            padj = padj,
            untreated = untreated,
            confoundingvariable = confoundingvariable
        shell:
            """
            ~/tools/DEseq2/bin/Rscript scripts/DESeq_Salmon.R \
                --infotable {params.infotable} \
                --design {params.design} \
                --lfc {params.lfc} \
                --padj {params.padj} \
                --untreated {params.untreated} \
                --confoundingvariable {params.confoundingvariable} \
                --Orthologgenes {input.Orthologgenes} 
            """

else:
    rule DESeq:
        input:
            quant =  expand("result/Salmon/{sample}/quant.sf", sample = samples),
            tx2gene = expand("result/Salmon/{species}_trans2symbol.tsv", species = specieses)
        output:
            all_genes = expand("result/Salmon/{pairname}.all.tsv", pairname = pairnames),
            diffexp_genes = expand("result/Salmon/{pairname}.diffexp.tsv", pairname = pairnames)
        params:
            infotable = infotable,
            design = design,
            lfc = lfc,
            padj = padj,
            untreated = untreated,
            confoundingvariable = confoundingvariable
        shell:
            """
            ~/tools/DEseq2/bin/Rscript scripts/DESeq_Salmon.R \
                --infotable {params.infotable} \
                --design {params.design} \
                --lfc {params.lfc} \
                --padj {params.padj} \
                --untreated {params.untreated} \
                --confoundingvariable {params.confoundingvariable}
            """



# rule salmon_expression_matrix:
#     input:
#         quant =  expand("result/Salmon/{sample}/quant.sf", sample = samples)
#     output:
#         genecount = "result/Salmon/{species}_gene_counts.tsv",
#         genecount_scaled = "result/Salmon/{species}_gene_counts_scaled.tsv",
#         tpm = "result/Salmon/{species}_gene_lengthScaledTPM.tsv",
#         genelength = "result/Salmon/{species}_gene_length.tsv"
#     params:
#         infotable = infotable,
#         gff = lambda wildcards: get_genome_gff(wildcards.species),
#         Species = "{species}",
#         indir = "result/Salmon/",
#         outdir = "result/Salmon/"
#     log:
#         "logs/salmon_expression_matrix.{species}.log"
#     shell:
#         """
#         {R453}
#         Rscript scripts/RNAseq_Salmon_ExpressionMatrix.R \
#             --infotable {params.infotable} \
#             --Species {params.Species} \
#             --gff {params.gff} \
#             --indir {params.indir} \
#             --outdir {params.outdir} \
#             2>{log}
#         """




######## modules ######

include: "rules/0.Common.smk"

def outfiles():
    files = []

    if "qc" in step or "all" in step:
        ### 原始数据质量控制
        include: "rules/1.QC.smk"
        files.append(expand("result/CleanData/{sample}.clean.{num}.fq.gz", sample = samples,num = [1,2]))
    # if "rsem" in step or "all" in step:
    # ### rsem 一步完成比对和统计，注意比对同样适用STAR工具，两步法里面STAR的参数参考rsem的参数
    #     include: "rules/rsem.smk"
    #     files.append(expand("result/RSEM/{sample}.genes.results", sample = samples))
    #     files.append(expand("result/RSEM/{sample}.isoforms.results", sample = samples))
    if "salmon" in step or "all" in step:
        include: "rules/salmon.smk"
        files.append(expand("result/Salmon/{sample}/quant.sf", sample = samples)) 
        #files.append(expand("result/Salmon/{species}_gene_counts_scaled.tsv", species = list(metainfo_dict_Species_to_Genome.keys())))
        #files.append(expand("result/Salmon/merged_gene_counts_scaled.tsv"))
        if design == "Species":
            files.append(expand("result/Salmon/{species}_trans2symbol.tsv",  species = list(metainfo_dict_Species_to_Genome.keys())))
            files.append(expand("result/Salmon/{pairname}.diffexp.tsv", pairname = pairnames))

        # if "deseq" in step or "all" in step:
        #     include: "rules/4.DEseq.smk"
        #     #files.append("result/DEseq/PCA.pdf")
        #     #files.append("result/DEseq/dds.rds")
        #     files.append(expand("result/DEseq/{pairname}.all.csv", pairname = pairnames))
        #     #files.append(expand("result/DEseq/{pairname}.Volcano.pdf", pairname = pairnames))

    # if pipeline == "twostep":
    #     if "map" in step or "all" in step:
    #         ### STAR 将质控后的 reads 比对到参考基因组
    #         include: "rules/2.Mapping.smk"
    #         files.append(expand("result/STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam", sample = samples))
    #         files.append(expand("result/STAR/{sample}/{sample}_Aligned.toTranscriptome.out.bam", sample = samples))
    #         # files.append(expand("result/STAR/{sample}/{sample}-bamqc-qualimap-report/qualimapReport.html", sample = samples))
    #     if "count" in step or "all" in step:
    #         include: "rules/3.Count.smk"
    #         files.append(expand("result/Count/{sample}/{sample}.salmon_quant/quant.sf", sample = samples))
    #         #files.append(expand("test/{sample}", sample = samples))


    if "annotate" in step or "all" in step:
        include: "rules/5.GOandKEGG.smk"
        files.append(gokegg_outputs)
        # files.append("result/WGCNA/WGCNAfile.rds")
    #print(files)
    return files


rule all:
    input:
        outfiles()


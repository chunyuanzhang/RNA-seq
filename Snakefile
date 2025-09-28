
######## modules ######

include: "rules/0.Common.smk"

def outfiles():
    files = []

    if step == "qc" or step == "all":
        ### 原始数据质量控制
        include: "rules/1.QC.smk"
        #files.append("result/QC/Raw/multiqc_report.html")
        files.append(expand("result/CleanData/{sample}.clean.{num}.fq.gz", sample = samples,num = [1,2]))
        #files.append("result/QC/Clean/multiqc_report.html")
    if pipeline == "onestep":
        ### rsem 一步完成比对和统计，注意比对同样适用STAR工具，两步法里面STAR的参数参考rsem的参数
        include: "rules/rsem.smk"
        files.append(expand("result/RSEM/{sample}.genes.results", sample = samples))
    if pipeline == "twostep":
        if step == "map" or step == "all":
            ### STAR 将质控后的 reads 比对到参考基因组
            include: "rules/2.Mapping.smk"
            files.append(expand("result/STAR/{sample}/{sample}_Aligned.out.bam", sample = samples))
            #files.append(expand("result/STAR/{sample}/{sample}_Aligned.toTranscriptome.out.bam", sample = samples))
            # files.append(expand("result/STAR/{sample}/{sample}-bamqc-qualimap-report/qualimapReport.html", sample = samples))
        if step == "count" or step == "all":
            include: "rules/3.Count.smk"
            #files.append(expand("result/Count/{sample}/{sample}.salmon_quant/quant.sf", sample = samples))
            files.append(expand("test/{sample}", sample = samples))
    if step == "deseq" or step == "all":
        include: "rules/4.DEseq.smk"
        files.append("result/DEseq/PCA.pdf")
        files.append("result/DEseq/diffgenes.tsv")
    # if step == "annotate" or step == "all":
    #     include: "rules/5.GOandKEGG.smk"
    #     files.append()
        # files.append(expand("result/STAR/{sample}/{sample}-bamqc-qualimap-report/qualimapReport.html", sample = samples))
        
    # print(files)
    return files


rule all:
    input:
        outfiles()



#################
# Quantification
#################


# rule create_salmon_sample_file:
#     input:
#         count = expand("result/STAR/{sample}/{sample}.salmon_quant/quant.sf", sample = samples)
#     output:
#         file = "samples.txt"
#     run:
#         with open(output.file, "w") as f:
#             for s in samples:
#                 f.write(s + "\t" + "result/STAR/" + s + "/" + s + ".salmon_quant/quant.sf" + "\t" + metainfo_dict[s] + "\n")


# rule DEseq2:
#     input:
#         count = expand("result/STAR/{sample}/{sample}.salmon_quant/quant.sf", sample = samples)
#     output:
#         PCA = "PCA.pdf"
#     run:
#         mappedfiles = ','.join(input.count)
#         cmd = (f"~/tools/DEseq2/bin/Rscript RNA-seq.R --metafile {infotable} --lfc {lfc} --pval {pval} --gtf {ref_gtf} --untreated {untreated} --mappedfiles {mappedfiles}")
#         print(cmd)
#         os.system(cmd)




import pandas as pd
import numpy as np
import yaml
import json
import os

# === 元信息 =========================================
infotable = "metafile"

# infotable文件内容示意，为了代码重复使用方便，infotable不要列
# 第一列是样本名，第二列是组名
# White_1W_F1.sra,White
# White_1W_F2.sra,White
# White_1W_F3.sra,White
# Wild_1W_F1.sra,Wild
# Wild_1W_F2.sra,Wild
# Wild_1W_F3.sra,Wild


# === 指定参考基因组 ==================================
## 本流程使用STAR进行比对，注意在比对前需要先手动对参考基因组建立索引
genomeDir = "/home/zhangchunyuan/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b_Ensembl/STARindex"
ref_fa = "/home/zhangchunyuan/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b_Ensembl/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa"
ref_gtf = "/home/zhangchunyuan/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b_Ensembl/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gtf"
ref_transscripts = "/home/zhangchunyuan/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b_Ensembl/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.transcripts.fa"

## 在基因组所在文件夹执行下行命令 双端150bp测序时, --sjdbOverhang 149
## STAR --runThreadN 64 --runMode genomeGenerate --genomeDir ./STARindex --genomeFastaFiles ./STARindex/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.fna --sjdbGTFfile  ./GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.gtf --sjdbOverhang 149 --genomeSAindexNbases 13

# === 外部参数处理 ===================================

if os.path.exists(infotable):
    metainfo = pd.read_csv(infotable, sep = ",", dtype=str, index_col = False, header=None)
    metainfo.dropna(how='all', inplace=True) # 删除可能存在的空行
    metainfo.columns = [f'V{i+1}' for i in range(len(metainfo.columns))]
    metainfo = metainfo[ ~np.array([v.startswith("#") for v in metainfo.V1.to_list()])]  # 删除被#注释掉的行
    # metainfo.index = metainfo.V1

    metainfo_dict = dict(zip(metainfo.V1, metainfo.V2)) 
    
    samples = metainfo.V1.tolist()

    # 设置wildcard约束
    wildcard_constraints:
        samples = "|".join(samples)
    

# === rules =========================================
# Determine the {outname} for all output files
# (SAMPLES,) = glob_wildcards("result/fastq/{sample}_1.fastq")

files = []
files.append(expand("result/STAR/{sample}/{sample}_Aligned.out.bam", sample = samples))
files.append(expand("result/STAR/{sample}/{sample}_Aligned.toTranscriptome.out.bam", sample = samples))
files.append(expand("result/STAR/{sample}/{sample}_Aligned.sorted.bam", sample = samples))
files.append(expand("result/STAR/{sample}/{sample}_Aligned.sorted.bam.bai", sample = samples))
files.append(expand("result/STAR/{sample}/{sample}_Aligned.sorted.flagstat", sample = samples))
files.append(expand("result/STAR/{sample}/{sample}-rnaseq-qualimap-report/qualimapReport.html", sample = samples))
files.append(expand("result/STAR/{sample}/{sample}-bamqc-qualimap-report/qualimapReport.html", sample = samples))
files.append(expand("result/STAR/{sample}/{sample}.salmon_quant/quant.sf", sample = samples))
files.append("samples.txt")



rule all:
    input:
        files

# rule sra2fastq:
#     input:
#         sra = "result/sra/{sample}.sra"
#     output:
#         fastq_1 = "result/fastq/{sample}_1.fastq",
#         fastq_2 = "result/fastq/{sample}_2.fastq"
#     shell:
#         """
#         fastq-dump --split-files --outdir result/fastq/ {input.sra}
#         """

# rule quality_evaluation:
#     input:
#         fastq_1 = "result/fastq/{sample}_1.fastq",
#         fastq_2 = "result/fastq/{sample}_2.fastq"
#     output:
#         fastq_1 = "result/fastq/{sample}_1_fastqc.html",
#         fastq_2 = "result/fastq/{sample}_2_fastqc.html"
#     threads:
#         16
#     shell:
#         """
#         fastqc {input.fastq_1} --threads {threads} --extract --delete
#         fastqc {input.fastq_2} --threads {threads} --extract --delete
#         """

# rule quality_summary:
#     input:
#         fastq_1 = expand("result/fastq/{sample}_1_fastqc.html", sample = samples ),
#         fastq_2 = expand("result/fastq/{sample}_2_fastqc.html", sample = samples )
#     output:    
#         summary = "result/fastq/multiqc_report.html"
#     shell:
#         """
#          multiqc --outdir result/fastq/ result/fastq/
#         """

# rule quality_control:
#     input:
#         fastq_1 = "result/fastq/{sample}_1.fastq",
#         fastq_2 = "result/fastq/{sample}_2.fastq",
#         summary = "result/fastq/multiqc_report.html"
#     output:
#         fastp_1 = "result/CleanData/{sample}_1.clean.fq.gz",
#         fastp_2 = "result/CleanData/{sample}_2.clean.fq.gz",
#         html = "result/CleanData/{sample}.html",
#         json = "result/CleanData/{sample}.json"
#     threads:
#         16
#     shell:
#         """
#         fastp -i {input.fastq_1} -o {output.fastp_1} \
#                 -I {input.fastq_2} -O {output.fastp_2} \
#                 -w {threads} \
#                 -h {output.html} \
#                 -j {output.json} \
#                 --detect_adapter_for_pe -l 25
#         """

# rule quality_reevaluation:
#     input:
#         fastp_1 = "result/CleanData/{sample}_1.clean.fq.gz",
#         fastp_2 = "result/CleanData/{sample}_2.clean.fq.gz"
#     output:
#         fastp_1 = "result/CleanData/{sample}_1.clean.html",
#         fastp_2 = "result/CleanData/{sample}_2.clean.html"
#     threads:
#         16
#     shell:
#         """
#         fastqc {input.fastp_1}  --threads {threads} --extract --delete
#         fastqc {input.fastp_2}  --threads {threads} --extract --delete
#         """

# rule quality_resummary:
#     input:
#         fastp_1 = expand("result/CleanData/{sample}_1.clean.html", sample = samples ),
#         fastp_2 = expand("result/CleanData/{sample}_2.clean.html", sample = samples )
#     output:
#         summary = "result/CleanData/multiqc_report.html"
#     shell:
#         """
#          multiqc --outdir result/CleanData/ result/CleanData/
#         """




rule map_reads_to_reference:
    input:
        fastp_1 = "result/CleanData/{sample}_1.clean.fq.gz",
        fastp_2 = "result/CleanData/{sample}_2.clean.fq.gz"
    output:
        bam = "result/STAR/{sample}/{sample}_Aligned.out.bam",
        bam2 = "result/STAR/{sample}/{sample}_Aligned.toTranscriptome.out.bam"
    params:
        "result/STAR/{sample}/{sample}_"
    threads:
        16
    shell:
        """
        STAR --runThreadN {threads} \
                --genomeDir {genomeDir} \
                --readFilesIn {input.fastp_1} {input.fastp_2} \
                --readFilesCommand zcat \
                --outFileNamePrefix {params} \
                --outSAMtype BAM Unsorted \
                --quantTranscriptomeBan Singleend \
                --outFilterType BySJout \
                --alignSJoverhangMin 8 --outFilterMultimapNmax 20 \
                --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
                --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 \
                --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
                --quantMode TranscriptomeSAM --outSAMattributes NH HI AS NM MD 
        """



#################
# Quality summary
#################



rule sort_bam:
    input:
        bam = "result/STAR/{sample}/{sample}_Aligned.out.bam"
    output:
        sortedbam = "result/STAR/{sample}/{sample}_Aligned.sorted.bam"
    params:
        prefix = "result/STAR/{sample}"
    threads:
        16
    shell:
        """
        # 该步骤同时跑太多会报错，将threads设置高一点，保证并行程序数量
        samtools sort -@ {threads} -T {params.prefix} {input.bam} > {output.sortedbam}
        
        """

rule index_bam:
    input:
        sortedbam = "result/STAR/{sample}/{sample}_Aligned.sorted.bam"
    output:
        bai = "result/STAR/{sample}/{sample}_Aligned.sorted.bam.bai",
        flagstat = "result/STAR/{sample}/{sample}_Aligned.sorted.flagstat"
    threads:
        16
    shell: 
        """
        samtools index -@ {threads} {input.sortedbam}
        samtools flagstat -@ {threads} {input.sortedbam} > {output.flagstat}
        """

rule qualimap_bamQC_report:
    input:
        sortedbam = "result/STAR/{sample}/{sample}_Aligned.sorted.bam"
    output:
        html = "result/STAR/{sample}/{sample}-bamqc-qualimap-report/qualimapReport.html"
    params:
        outdir = "result/STAR/{sample}/{sample}-bamqc-qualimap-report"
    threads:
        32
    shell:
        """
        qualimap bamqc -bam {input.sortedbam} -gff {ref_gtf} -outdir {params.outdir} --java-mem-size=16G -nt {threads} 
        """

rule create_sampletxt_file:
    input: 
        html = expand("result/STAR/{sample}/{sample}-bamqc-qualimap-report/qualimapReport.html", sample = samples)
    output:
        sampletxt = "sample.txt"
    run:
        with open(output.sampletxt, "w") as f:
            for s in samples:
                f.write(s + "\t" + "result/STAR/" + s + "/" + s + "-bamqc-qualimap-report" + "\t" + metainfo_dict[s] + "\n")
        os.system('qualimap multi-bamqc -d sample.txt')






rule qualimap_rnaseqQC_report:
    input:
        sortedbam = "result/STAR/{sample}/{sample}_Aligned.sorted.bam"
    output:
        html = "result/STAR/{sample}/{sample}-rnaseq-qualimap-report/qualimapReport.html"
    params:
        outdir = "result/STAR/{sample}/{sample}-rnaseq-qualimap-report"
    threads:
        32
    shell:
        """
        # 该步骤没有修改线程的参数，但是看起来比较吃资源,因此用较高的线程控制不要占太多
        qualimap rnaseq -bam {input.sortedbam} -gtf {ref_gtf} -outdir {params.outdir} --java-mem-size=16G
        """




#################
# Quantification
#################

rule salmon_quant:
    input:
        bam = "result/STAR/{sample}/{sample}_Aligned.toTranscriptome.out.bam"
    output:
        quant = "result/STAR/{sample}/{sample}.salmon_quant/quant.sf"
    threads:
        16
    params:
        ourdir = "result/STAR/{sample}/{sample}.salmon_quant"
    shell:
        """
        salmon quant -t {ref_transscripts} --libType A -a {input.bam} -o {params.ourdir} --threads {threads}
        """

rule create_salmon_sample_file:
    input:
        count = expand("result/STAR/{sample}/{sample}.salmon_quant/quant.sf", sample = samples)
    output:
        file = "samples.txt"
    run:
        with open(output.file, "w") as f:
            for s in samples:
                f.write(s + "\t" + "result/STAR/" + s + "/" + s + ".salmon_quant/quant.sf" + "\t" + metainfo_dict[s] + "\n")





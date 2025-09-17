
######## modules ######

include: "rules/0.Common.smk"


files = []
### 原始数据质量控制
files.append(expand("result/QC/{sample}_1_fastqc.html", sample = samples))

files.append(expand("result/STAR/{sample}/{sample}_Aligned.out.bam", sample = samples))
files.append(expand("result/STAR/{sample}/{sample}_Aligned.toTranscriptome.out.bam", sample = samples))
files.append(expand("result/STAR/{sample}/{sample}_Aligned.sorted.bam", sample = samples))
files.append(expand("result/STAR/{sample}/{sample}_Aligned.sorted.bam.bai", sample = samples))
files.append(expand("result/STAR/{sample}/{sample}_Aligned.sorted.flagstat", sample = samples))
files.append(expand("result/STAR/{sample}/{sample}-rnaseq-qualimap-report/qualimapReport.html", sample = samples))
files.append(expand("result/STAR/{sample}/{sample}-bamqc-qualimap-report/qualimapReport.html", sample = samples))
files.append(expand("result/STAR/{sample}/{sample}.salmon_quant/quant.sf", sample = samples))
files.append("PCA.pdf")



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
#         fastq_1 = "result/RowData/{sample}_1.fastq",
#         fastq_2 = "result/RowData/{sample}_2.fastq"
#     output:
#         fastq_1 = "result/QC/{sample}_1_fastqc.html",
#         fastq_2 = "result/QC/{sample}_2_fastqc.html"
#     threads:
#         16
#     shell:
#         """
#         fastqc {input.fastq_1} --threads {threads} --extract --delete --outdir result/QC
#         fastqc {input.fastq_2} --threads {threads} --extract --delete --outdir result/QC
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

# rule create_salmon_sample_file:
#     input:
#         count = expand("result/STAR/{sample}/{sample}.salmon_quant/quant.sf", sample = samples)
#     output:
#         file = "samples.txt"
#     run:
#         with open(output.file, "w") as f:
#             for s in samples:
#                 f.write(s + "\t" + "result/STAR/" + s + "/" + s + ".salmon_quant/quant.sf" + "\t" + metainfo_dict[s] + "\n")


rule DEseq2:
    input:
        count = expand("result/STAR/{sample}/{sample}.salmon_quant/quant.sf", sample = samples)
    output:
        PCA = "PCA.pdf"
    run:
        mappedfiles = ','.join(input.count)
        cmd = (f"~/tools/DEseq2/bin/Rscript RNA-seq.R --metafile {infotable} --lfc {lfc} --pval {pval} --gtf {ref_gtf} --untreated {untreated} --mappedfiles {mappedfiles}")
        print(cmd)
        os.system(cmd)



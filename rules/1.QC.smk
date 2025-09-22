

# rule sra2fastq:
#     input:
#         sra = "result/sra/{sample}.sra"
#     output:
#         fastq_1 = "result/RawData/{sample}.raw.1.fq.gz",
#         fastq_2 = "result/RawData/{sample}.raw.2.fq.gz"
#     shell:
#         """
#         fastq-dump --split-files --outdir result/fastq/ {input.sra} 
#         """


rule quality_evaluation:
    input:
        rawdata_1 = "result/RawData/{sample}.raw.1.fq.gz",
        rawdata_2 = "result/RawData/{sample}.raw.2.fq.gz"
    output:
        QC_1 = temp("result/QC/Raw/{sample}.raw.1_fastqc.html"),
        QC_2 = temp("result/QC/Raw/{sample}.raw.2_fastqc.html")
    threads:
        16
    params:
        outdir = "result/QC/Raw"
    shell:
        """
        fastqc {input.rawdata_1} --threads {threads} --extract --delete --outdir {params.outdir}
        fastqc {input.rawdata_2} --threads {threads} --extract --delete --outdir {params.outdir}
        """

rule quality_summary:
    input:
        QC_1 = expand("result/QC/Raw/{sample}.raw.1_fastqc.html", sample = samples ),
        QC_2 = expand("result/QC/Raw/{sample}.raw.2_fastqc.html", sample = samples )
    output:    
        summary = "result/QC/Raw/multiqc_report.html"
    shell:
        """
         multiqc --outdir result/QC/Raw result/QC/Raw
        """

rule quality_control:
    input:
        rawdata_1 = "result/RawData/{sample}.raw.1.fq.gz",
        rawdata_2 = "result/RawData/{sample}.raw.2.fq.gz",
        summary = "result/QC/Raw/multiqc_report.html"
    output:
        cleandata_1 = "result/CleanData/{sample}.clean.1.fq.gz",
        cleandata_2 = "result/CleanData/{sample}.clean.2.fq.gz",
        html = "result/QC/fastp/{sample}.clean.html",
        json = "result/QC/fastp/{sample}.clean.json"
    threads:
        16
    shell:
        """
        fastp -i {input.rawdata_1} -o {output.cleandata_1} \
                -I {input.rawdata_2} -O {output.cleandata_2} \
                -w {threads} \
                -h {output.html} \
                -j {output.json} \
                --detect_adapter_for_pe -l 25
        """


rule quality_reevaluation:
    input:
        cleandata_1 = "result/CleanData/{sample}.clean.1.fq.gz",
        cleandata_2 = "result/CleanData/{sample}.clean.2.fq.gz",
    output:
        QC_1 = temp("result/QC/Clean/{sample}.clean.1_fastqc.html"),
        QC_2 = temp("result/QC/Clean/{sample}.clean.2_fastqc.html")
    threads:
        16
    params:
        outdir = "result/QC/Clean"
    shell:
        """
        fastqc {input.cleandata_1}  --threads {threads} --extract --delete --outdir {params.outdir}
        fastqc {input.cleandata_2}  --threads {threads} --extract --delete --outdir {params.outdir}
        """

rule quality_resummary:
    input:
        QC_1 = expand("result/QC/Clean/{sample}.clean.1_fastqc.html", sample = samples ),
        QC_2 = expand("result/QC/Clean/{sample}.clean.2_fastqc.html", sample = samples )
    output:
        summary = "result/QC/Clean/multiqc_report.html"
    shell:
        """
         multiqc --outdir result/QC/Clean result/QC/Clean
        """

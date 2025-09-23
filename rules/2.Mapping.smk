

rule map_reads_to_reference:
    input:
        cleandata_1 = "result/CleanData/{sample}.clean.1.fq.gz",
        cleandata_2 = "result/CleanData/{sample}.clean.2.fq.gz",
    output:
        bam1 = "result/STAR/{sample}/{sample}_Aligned.out.bam"
        # bam2 = "result/STAR/{sample}/{sample}_Aligned.toTranscriptome.out.bam"
    params:
        prefix = "result/STAR/{sample}/{sample}_",
        genomeDir =  lambda wildcards:  referenceDir + metainfo_dict_Genome[wildcards.sample] + "/STARindex"
    threads:
        100
    shell:
        """
        STAR --genomeDir {params.genomeDir} \
            --outSAMunmapped Within \
            --outFilterType BySJout \
            --outSAMattributes NH HI AS NM MD \
            --outFilterMultimapNmax 20 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --sjdbScore 1 \
            --runThreadN 96 \
            --genomeLoad NoSharedMemory \
            --outSAMtype BAM Unsorted \
            --quantMode TranscriptomeSAM \
            --outSAMheaderHD @HD VN:1.4 SO:unsorted \
            --outFileNamePrefix {params.prefix} \
            --readFilesCommand zcat \
            --readFilesIn {input.cleandata_1} {input.cleandata_2}
        """

# 对比对质量进行评估
# rule qualimap_bamQC_report1:
#     input:
#         bam1 = "result/STAR/{sample}/{sample}_Aligned.out.bam"
#     output:
#         html = "result/STAR/{sample}/{sample}-bamqc-qualimap-report/qualimapReport.html"
#     params:
#         outdir = "result/STAR/{sample}/{sample}-bamqc-qualimap-report",
#         referenceGtf = lambda wildcards: glob.glob(referenceDir + metainfo_dict_Genome[wildcards.sample] + "/*.gtf" )
#     threads:
#         32
#     shell:
#         """
#         qualimap bamqc -bam {input.bam1} -gff {params.referenceGtf} -outdir {params.outdir} --java-mem-size=16G -nt {threads} 
#         """


# rule qualimap_bamQC_report2:
#     input:
#         bam1 = "result/STAR/{sample}/{sample}_Aligned.toTranscriptome.out.bam"
#     output:
#         html = "result/STAR/{sample}/{sample}-bamqc-qualimap-report/qualimapReport.html"
#     params:
#         outdir = "result/STAR/{sample}/{sample}-bamqc-qualimap-report",
#         referenceGtf = lambda wildcards: glob.glob(referenceDir + metainfo_dict_Genome[wildcards.sample] + "/*.gtf" )
#     threads:
#         32
#     shell:
#         """
#         qualimap bamqc -bam {input.bam1} -gff {params.referenceGtf} -outdir {params.outdir} --java-mem-size=16G -nt {threads} 
#         """


# rule create_sampletxt_file:
#     input: 
#         html = expand("result/STAR/{sample}/{sample}-bamqc-qualimap-report/qualimapReport.html", sample = samples)
#     output:
#         sampletxt = "sample.txt"
#     run:
#         with open(output.sampletxt, "w") as f:
#             for s in samples:
#                 f.write(s + "\t" + "result/STAR/" + s + "/" + s + "-bamqc-qualimap-report" + "\t" + metainfo_dict[s] + "\n")
#         os.system('qualimap multi-bamqc -d sample.txt')



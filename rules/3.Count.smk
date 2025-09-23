## 使用奉化水鸭统计表达量时，STAR检测到400多个参考基因组没有的转录本，约2.2%的reads组成这些新转录本
# rule salmon_quant:
#     input:
#         bam = "result/STAR/{sample}/{sample}_Aligned.toTranscriptome.out.bam"
#     output:
#         quant = "result/Count/{sample}/{sample}.salmon_quant/quant.sf"
#     threads:
#         64
#     params:
#         ourdir = "result/Count/{sample}/{sample}.salmon_quant",
#         ref_transscripts = lambda wildcards: glob.glob(referenceDir + metainfo_dict_Genome[wildcards.sample] + "/*transcripts.fna" )
#     shell:
#         """
#         salmon quant -t {params.ref_transscripts} --libType A -a {input.bam} -o {params.ourdir} --threads {threads}
#         """


rule rsem_count:
    input:
        bam = "result/STAR/{sample}/{sample}_Aligned.toTranscriptome.out.bam"
    output:
        test = "test/{sample}"
    threads:
        64
    params:
        prefix = "result/RSEM/{sample}"
    shell:
        """
        ~/tools/rsem/bin/rsem-calculate-expression --alignments  --paired-end \
            --no-bam-output \
            -p {threads} \
            {input.bam} \
            ~/zhangchunyuan/reference/IASCAAS_PekinDuck_T2T/RSEMindex/IASCAAS_PekinDuck_T2T \
            {params.prefix}
        """



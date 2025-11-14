rule rsem_mapping_and_counting:
    input:
        cleandata_1 = "result/CleanData/{sample}.clean.1.fq.gz",
        cleandata_2 = "result/CleanData/{sample}.clean.2.fq.gz"
    output:
        genes = "result/RSEM/{sample}.genes.results",
        isoforms = "result/RSEM/{sample}.isoforms.results",
        transcript = "result/RSEM/{sample}.transcript.bam"
    threads:
        96
    params:
        rsem_reference = lambda wildcards:  referenceDir + metainfo_dict_Genome[wildcards.sample] + "/RSEMindex/" + metainfo_dict_Genome[wildcards.sample] 
    log:
        "logs/rsem_mapping_and_counting/{sample}.log"
    shell:
        """
        {STAR}
         ~/tools/rsem/bin/rsem-calculate-expression \
            -p {threads} \
            --paired-end \
            --star \
            --star-gzipped-read-file \
            --estimate-rspd \
            --append-names \
            {input.cleandata_1} {input.cleandata_2} \
            {params.rsem_reference} \
            result/RSEM/{wildcards.sample} 2>{log}
        """

# 对比对质量进行评估
rule qualimap_bamQC_report1:
    input:
        bam1 = "result/STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        html = "result/STAR/{sample}/{sample}-bamqc-qualimap-report/qualimapReport.html"
    params:
        outdir = "result/STAR/{sample}/{sample}-bamqc-qualimap-report",
        referenceGtf = lambda wildcards: glob.glob(referenceDir + metainfo_dict_Genome[wildcards.sample] + "/*.gtf" )
    threads:
        32
    log:
        "logs/qualimap_bamQC_report1/{sample}.log"
    shell:
        """
        {qualimap}
        qualimap bamqc -bam {input.bam1} -gff {params.referenceGtf} -outdir {params.outdir} --java-mem-size=16G -nt {threads} 2>{log}
        """

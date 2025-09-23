
if pipeline == "onestep":
    rule DEseq_analysis:
        input:  
            isoforms = expand("result/RSEM/{sample}.isoforms.results", sample = samples),
            genes = expand("result/RSEM/{sample}.genes.results", sample = samples)
        output:
            PCA = "result/DEseq/PCA.pdf"
        params:
            isoforms = lambda wildcards, input: ",".join(input.isoforms),
            genes = lambda wildcards, input: ",".join(input.genes)
        shell:
            """
            ~/tools/DEseq2/bin/Rscript scripts/RNA-seq.R \
                --infotable {infotable} \
                --lfc {lfc}  --pval {pval}  --untreated {untreated} \
                --CountingMethod rsem \
                --CountingFiles.isoforms {params.isoforms} \
                --CountingFiles.genes {params.genes}
            """
        

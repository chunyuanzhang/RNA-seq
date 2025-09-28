
if interspecies is not True:
    if pipeline == "onestep":
        rule DEseq_analysis:
            input:  
                isoforms = expand("result/RSEM/{sample}.isoforms.results", sample = samples)
            output:
                PCA = "result/DEseq/PCA.pdf",
                diffgenes = "result/DEseq/diffgenes.tsv"
            params:
                isoforms = lambda wildcards, input: ",".join(input.isoforms)
            shell:
                """
                ~/tools/DEseq2/bin/Rscript scripts/DEseq.R  \
                    --infotable {infotable} \
                    --lfc {lfc}  --pval {pval}  --untreated {untreated} \
                    --CountingMethod rsem \
                    --CountingFiles.isoforms {params.isoforms} \
                    --CountingFiles.genes {params.genes}
                """

else:
    if pipeline == "onestep":
        rule DEseq_analysis_interspecies:
            input:
                isoforms = expand("result/RSEM/{sample}.isoforms.results", sample = samples)
            output:
                PCA = "result/DEseq/PCA.pdf",
                diffgenes = "result/DEseq/diffgenes.tsv"
            params:
                isoforms = lambda wildcards, input: ",".join(input.isoforms)
            shell:
                """
                ~/tools/DEseq2/bin/Rscript scripts/DEseq-interspecies.R \
                    --infotable {infotable} \
                    --lfc {lfc} \
                    --pval {pval} \
                    --untreated {untreated} \
                    --CountingFiles.isoforms {params.isoforms} \
                    --orthologgenes result/Ortho_chicken_vs_duck/one_to_one_orthologgenes.txt
                """

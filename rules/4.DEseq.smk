
if pipeline == "onestep":
    rule DEseq_analysis_interspecies:
        input:
            isoforms = expand("result/RSEM/{sample}.isoforms.results", sample = samples)
        output:
            PCA = "result/DEseq/PCA.pdf",
            Volcano = expand("result/DEseq/{pairname}.Volcano.pdf", pairname = pairnames),
            diffgenes = expand("result/DEseq/{pairname}.all.csv", pairname = pairnames)
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
                --orthologgenes {Orthologgenes}
            """





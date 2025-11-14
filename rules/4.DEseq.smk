
if pipeline == "onestep":
    rule DEseq_analysis:
        input:
            isoforms = expand("result/RSEM/{sample}.isoforms.results", sample = samples)
        output:
            rds = "result/DEseq/dds.rds",
            PCA = "result/DEseq/PCA.pdf",
            Volcano = expand("result/DEseq/{pairname}.Volcano.pdf", pairname = pairnames),
            diffgenes = expand("result/DEseq/{pairname}.all.csv", pairname = pairnames)
        params:
            isoforms = lambda wildcards, input: ",".join(input.isoforms),
            Orthologgenes = lambda wildcards: globals().get('Orthologgenes', 'None')
        log:
            "logs/DEseq_analysis.log"
        shell:
            """
            ~/tools/DEseq2/bin/Rscript scripts/DEseq.R \
                --infotable {infotable} \
                --lfc {lfc} \
                --pval {pval} \
                --design {design} \
                --untreated {untreated} \
                --CountingFiles.isoforms {params.isoforms} \
                --orthologgenes {params.Orthologgenes} 2>{log}
            """





if pipeline == "onestep":
    rule GOandKEGG:
        input:
            diffgenes = "result/DEseq/diffgenes.tsv"
        output:
            "test.done
        shell:
            """
            ~/tools/DEseq2/bin/Rscript scripts/GOandKEGG.R \
                
            """
        
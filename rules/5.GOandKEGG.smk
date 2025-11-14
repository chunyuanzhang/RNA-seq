
rule GOandKEGG:
    input: 
        diffgenes = "result/DEseq/{pairname}.all.csv"
    output:
        KEGG = "result/DEseq/{pairname}.{loop}.KEGG.csv",
        GO = "result/DEseq/{pairname}.{loop}.GO.csv"
    params:
        referenceGtf = get_gtf_by_design_value,
        emapperannotations = get_emapper_by_design_value,
        pval = pval
    log:
        "logs/GOandKEGG/{pairname}.{loop}.log"
    shell:
        """
        ~/tools/rstudio/bin/Rscript scripts/anno.R \
            --target {wildcards.loop} \
            --gtf {params.referenceGtf} \
            --diffgenes {input.diffgenes} \
            --backgroundgenes {input.diffgenes} \
            --pval {params.pval} \
            --emapperannotations {params.emapperannotations} 2>{log}
        """

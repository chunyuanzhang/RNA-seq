if interspecies is not True:
    if pipeline == "onestep":
        rule GOandKEGG:
            input:
                diffgenes = "result/DEseq/diffgenes.csv"
            output:
                "test.done"
            shell:
                """
                ~/tools/DEseq2/bin/Rscript scripts/anno.R \
                    
                """
else:
    if pipeline == "onestep":
        rule GOandKEGG_interspecies:
            input: 
                diffgenes = "result/DEseq/diffgenes.csv"
            output:
                KEGG = "result/DEseq/{specie}.KEGG.csv"
            params:
                referenceGtf = lambda wildcards: glob.glob(referenceDir + metainfo_dict_Species_to_Genome.get(wildcards.specie) + "/*.gtf"),
                emapperannotations = lambda wildcards: glob.glob(referenceDir + metainfo_dict_Species_to_Genome.get(wildcards.specie) + "/*.emapper.annotations")
            shell:
                """
                ~/tools/DEseq2/bin/Rscript scripts/anno.R \
                    --infotable {infotable} \
                    --lfc {lfc} \
                    --pval {pval} \
                    --untreated {untreated} \
                    --target {wildcards.specie} \
                    --gtf {params.referenceGtf} \
                    --diffgenes {input.diffgenes} \
                    --emapperannotations {params.emapperannotations}
                """

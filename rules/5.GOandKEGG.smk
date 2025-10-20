if design == "Species":
    if pipeline == "onestep":
        rule GOandKEGG_interspecies:
            input: 
                diffgenes = "result/DEseq/{pairname}.all.csv"
            output:
                KEGG = "result/DEseq/{pairname}.{loop}.KEGG.csv",
                GO = "result/DEseq/{pairname}.{loop}.GO.csv"
            params:
                referenceGtf = lambda wildcards: glob.glob(referenceDir + metainfo_dict_Species_to_Genome.get(wildcards.loop) + "/*.gtf"),
                emapperannotations = lambda wildcards: glob.glob(referenceDir + metainfo_dict_Species_to_Genome.get(wildcards.loop) + "/*.emapper.annotations")
            shell:
                """
                ~/tools/rstudio/bin/Rscript scripts/anno.R \
                    --target {wildcards.loop} \
                    --gtf {params.referenceGtf} \
                    --diffgenes {input.diffgenes} \
                    --backgroundgenes {input.diffgenes} \
                    --pval {pval} \
                    --emapperannotations {params.emapperannotations}
                """
else:
    if pipeline == "onestep":
        rule GOandKEGG:
            input:
                diffgenes = "result/DEseq/{pairname}.all.csv"
            output:
                KEGG = "result/DEseq/{pairname}.{loop}.KEGG.csv",
                GO = "result/DEseq/{pairname}.{loop}.GO.csv"
            shell:
                """
                ~/tools/DEseq2/bin/Rscript scripts/anno.R \
                    
                """

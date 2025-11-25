
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

    rule WGCNA:
        input:
            rds = "result/DEseq/dds.rds",
            gtfchicken = "/home/zhangchunyuan/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf",
            emapperannotationschicken = "/home/zhangchunyuan/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b/chicken.emapper.annotations",
            gtfduck = "/home/zhangchunyuan/zhangchunyuan/reference/IASCAAS_PekinDuck_T2T/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.gtf",
            emapperannotationsduck = "/home/zhangchunyuan/zhangchunyuan/reference/IASCAAS_PekinDuck_T2T/duck.emapper.annotations"
        output:
            WGCNA = "result/WGCNA/WGCNAfile.rds"
        shell:
            """
            ~/tools/rstudio/bin/Rscript scripts/WGCNA.R \
                --rds {input.rds} \
                --gtf.chicken {input.gtfchicken} \
                --gtf.duck {input.gtfduck} \
                --emapperannotations.chicken {input.emapperannotationschicken} \
                --emapperannotations.duck {input.emapperannotationsduck}
            """



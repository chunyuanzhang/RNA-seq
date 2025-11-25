suppressMessages(library(optparse))
suppressMessages(library(dynamicTreeCut))
suppressMessages(library(fastcluster))
suppressMessages(library(WGCNA))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
# 启用多线程
enableWGCNAThreads()

source("scripts/functions.R")

# ==============================================================
# 读取外部参数
# ==============================================================

option_list <- list(
  make_option("--rds", type="character", default="result/DEseq/dds.rds", help="差异分析的DEseq对象"),
  make_option("--gtf.chicken",type="character", default = NULL, help = "gtf文件"),
  make_option("--emapperannotations.chicken", type="character", default=NULL, help="emapper比对结果文件，在没有现成的注释数据库时使用"),
  make_option("--gtf.duck",type="character", default = NULL, help = "gtf文件"),
  make_option("--emapperannotations.duck", type="character", default=NULL, help="emapper比对结果文件，在没有现成的注释数据库时使用"),
  make_option("--outdir", type = "character", default = "result/WGCNA/", help = "差异分析结果存储路径")
)

args <- parse_args(OptionParser(option_list=option_list))

rds <- args$rds
gtf.chicken <- args$gtf.chicken
emapperannotations.chicken <- args$emapperannotations.chicken
gtf.duck <- args$gtf.duck
emapperannotations.duck <- args$emapperannotations.duck
outdir <- args$outdir



# ==============================================================
# 读取dds数据
# ==============================================================

dds <- readRDS(rds)

# ==============================================================
# 在 WGCNA 网络分析之前，需要先过滤掉低方差的基因，这些基因往往没有特定的生物学功能，还会增加分析噪音
# ==============================================================

datExpr.chicken <- prepare_data_for_wgcna(dds = dds[, dds$Species == "Chicken"], quantile = 0.75, min_count = 5)
datExpr.duck <- prepare_data_for_wgcna(dds = dds[, dds$Species == "Duck"], quantile = 0.75, min_count = 5)


# ==============================================================
# 通过绘制样本树筛选样本
# ==============================================================

# sampleTree <- hclust(dist(t(datExpr)), method = "average")
# plot(sampleTree)

removesample <- c("MA-C-4", "BR-C-1")
datExpr.chicken <- datExpr.chicken[, !colnames(datExpr.chicken) %in% "BR-C-1"]
datExpr.duck <- datExpr.duck[, !colnames(datExpr.duck) %in% "MA-C-4"]

# ==============================================================
# 通过绘制样本树筛选样本
# ==============================================================

softPower.chicken <- auto_select_soft_Power(datExpr = t(datExpr.chicken))
softPower.duck <- auto_select_soft_Power(datExpr = t(datExpr.duck))

# ==============================================================
# 构建网络，并绘制模块图
# ==============================================================

net.chicken <- net_and_plot(datExpr = datExpr.chicken, softPower = softPower.chicken, labeltext = "Chicken",
                            pdffile = paste0(outdir, "Chicken.module.pdf"))
net.duck <- net_and_plot(datExpr = datExpr.duck, softPower = softPower.duck, labeltext = "Duck",
                         pdffile = paste0(outdir, "Duck.module.pdf" ))


# ==============================================================
# 模块的相关性
# ==============================================================

MEs.chicken <- moduleEigengenes(t(datExpr.chicken), colors = labels2colors(net.chicken$colors))$eigengenes %>% orderMEs()
colnames(MEs.chicken) <- colnames(MEs.chicken) %>% gsub("ME","", .)
pdf(paste0(outdir, "Chicken.cormodule.pdf"), width = 4, height = 4)
pheatmap::pheatmap(mat = cor(MEs.chicken))
dev.off()

MEs.duck <- moduleEigengenes(t(datExpr.duck), colors = labels2colors(net.duck$colors))$eigengenes %>% orderMEs()
colnames(MEs.duck) <- colnames(MEs.duck) %>% gsub("ME","", .)
pdf(paste0(outdir,"Duck.cormodule.pdf"), width = 4, height = 4)
pheatmap::pheatmap(mat = cor(MEs.duck))
dev.off()

# ==============================================================
#  Module preservation analysis
# ==============================================================

multiData <- list(
  Chicken = list(data = t(datExpr.chicken)), # Transpose to make samples columns
  Duck = list(data = t(datExpr.duck)) # Transpose to make samples columns
)

multiColor <- list(
  Chicken = labels2colors(net.chicken$colors) %>% setNames(net.chicken$colors),
  Duck = labels2colors(net.duck$colors) %>% setNames(net.duck$colors)
)

cat("\n\n\nPreservation analysis Chicken #################\n\n\n")

preservation_chicken_vs_duck <- WGCNA::modulePreservation(
  multiData = multiData,
  multiColor = multiColor,
  referenceNetworks = 1,  # 将第一个数据集作为 reference
  nPermutations = 200, 
  randomSeed = 10086,
  parallelCalculation = TRUE,
  verbose = 3
)

cat("\n\n\nPreservation analysis Duck #################\n\n\n")

preservation_duck_vs_chicken <- WGCNA::modulePreservation(
  multiData = multiData,
  multiColor = multiColor,
  referenceNetworks = 2,  # 将第一个数据集作为 reference
  nPermutations = 200, 
  randomSeed = 10086,
  parallelCalculation = TRUE,
  verbose = 3
)

preservation_stats <- preservation_chicken_vs_duck$preservation$Z$ref.Chicken$inColumnsAlsoPresentIn.Duck %>%
  mutate(Module = rownames(.)) %>%
  dplyr::filter(Module != "gold") 
p.pre.chicken <- ggplot(data = preservation_stats, aes(x = Module, y = Zsummary.pres, fill = Module))  + geom_bar(stat = "identity") +
  scale_fill_manual(values = preservation_stats$Module) + ggtitle("Chicken as ref")
ggsave(filename = paste0(outdir, "Chicken.preservation.pdf"), plot = p.pre.chicken, width = 6, height = 4)


preservation_stats.duck <- preservation_duck_vs_chicken$preservation$Z$ref.Duck$inColumnsAlsoPresentIn.Chicken %>%
  mutate(Module = rownames(.)) %>%
  dplyr::filter(Module != "gold") 
p.pre.duck <- ggplot(data = preservation_stats.duck, aes(x = Module, y = Zsummary.pres, fill = Module))  + geom_bar(stat = "identity") +
  scale_fill_manual(values = preservation_stats.duck$Module) + ggtitle("Duck as ref")
ggsave(filename = paste0(outdir, "Duck.preservation.pdf"), plot = p.pre.duck, width = 6, height = 4)


# ==============================================================
#  用列表储存模块基因
# ==============================================================

genes_in_module.chicken <- split(net.chicken$colors %>% names, f = labels2colors(net.chicken$colors))
genes_in_module.duck <- split(net.duck$colors %>% names, f = labels2colors(net.duck$colors))

# ==============================================================
#  GO 和 KEGG 注释
# ==============================================================

# 各模块 GO 注释
enrichment_GO.chicken <- list()
enrichment_KEGG.chicken <- list()

for (module in genes_in_module.chicken %>% names()) {
  
  cat("\t模块 ",module," 功能注释")
  enrichment_GO.chicken[[module]] <- GO_enrichment(
    diffgenes = extract_gene_from_gene_vs_gene_pair(genes_in_module.chicken[[module]], side = 1), 
    backgroundgenes = extract_gene_from_gene_vs_gene_pair(SummarizedExperiment::assay(dds) %>% rownames(), side = 1),
    GTF = gtf.chicken, emapperannotations = emapperannotations.chicken
  )
  
  enrichment_KEGG.chicken[[module]] <- KEGG_enrichment(
    diffgenes = extract_gene_from_gene_vs_gene_pair(genes_in_module.chicken[[module]], side = 1), 
    backgroundgenes = extract_gene_from_gene_vs_gene_pair(SummarizedExperiment::assay(dds) %>% rownames(), side = 1),
    GTF = gtf.chicken, emapperannotations = emapperannotations.chicken
  )
  
}

# 各模块 KEGG 注释
enrichment_GO.duck <- list()
enrichment_KEGG.duck <- list()
for (module in genes_in_module.duck %>% names()) {

  enrichment_GO.duck[[module]] <- GO_enrichment(
    diffgenes = extract_gene_from_gene_vs_gene_pair(genes_in_module.duck[[module]], side = 2), 
    backgroundgenes = extract_gene_from_gene_vs_gene_pair(SummarizedExperiment::assay(dds) %>% rownames(), side = 2),
    GTF = gtf.duck, emapperannotations = emapperannotations.duck
  )
  
  enrichment_KEGG.duck[[module]] <- KEGG_enrichment(
    diffgenes = extract_gene_from_gene_vs_gene_pair(genes_in_module.duck[[module]], side = 2), 
    backgroundgenes = extract_gene_from_gene_vs_gene_pair(SummarizedExperiment::assay(dds) %>% rownames(), side = 2),
    GTF = gtf.duck, emapperannotations = emapperannotations.duck
  ) 
  
}



### 绘制点图

for (module in genes_in_module.chicken %>% names()) {

  p.GO <- anno_plot(object = enrichment_GO.chicken[[module]], title = paste0("GO: Chicken module ", module))
  p.KEGG <- anno_plot(object = enrichment_KEGG.chicken[[module]], title = paste0("KEGG: Chicken module ", module))
  p <- ggpubr::ggarrange(p.GO, p.KEGG, labels = c("A", "B"), ncol = 2)
  ggsave(filename = paste0(outdir, "Chicken.",module,".anno.pdf" ), plot = p, width = 14, height = 6)
}

for (module in genes_in_module.duck %>% names()) {
  cat("模块 ",module, "\n")
  p.GO <- anno_plot(object = enrichment_GO.duck[[module]], title = paste0("GO: Duck module ", module))
  p.KEGG <- anno_plot(object = enrichment_KEGG.duck[[module]], title = paste0("KEGG: Duck module ", module))
  p <- ggpubr::ggarrange(p.GO, p.KEGG, labels = c("A", "B"), ncol = 2)
  ggsave(filename = paste0(outdir, "Duck.",module,".anno.pdf" ), plot = p, width = 14, height = 6)
}

# ==============================================================
#  将数据作成list对象输出，方便后续调用
# ==============================================================


WGCNAfile <- list(
  "Chicken" = list(
    "genes" = genes_in_module.chicken,
    "enrichment_GO" = enrichment_GO.chicken,
    "enrichment_KEGG" = enrichment_KEGG.chicken,
    "preservation" = preservation_chicken_vs_duck
  ),
  "Duck" = list(
    "genes" = genes_in_module.duck,
    "enrichment_GO" = enrichment_GO.duck,
    "enrichment_KEGG" = enrichment_KEGG.duck,
    "preservation" = preservation_duck_vs_chicken
  )
)

saveRDS(WGCNAfile, file = paste0(outdir,"WGCNAfile.rds"))



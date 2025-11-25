
library(dynamicTreeCut)
library(fastcluster)
library(WGCNA)
library(dplyr)
library(ggplot2)
# library(DESeq2)
# library(pheatmap)
# library(RColorBrewer)
# library(corrplot)
# library(igraph)
# library(networkD3)
# library(tibble)
# library(reshape2)
# library(clusterProfiler)
# 启用多线程（如果支持）
enableWGCNAThreads()

source("~/Desktop/02.Coding/RNA-seq/scripts/functions.R")

dds <- readRDS("~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/dds.rds")

# datExpr <- counts(dds, normalized = TRUE)
# chickens <- sample_info %>% as.data.frame() %>% filter(Species == "Chicken") %>% pull(SampleID) 
# ducks <- sample_info %>% as.data.frame() %>% filter(Species == "Duck") %>% pull(SampleID) 


# ===============================
# 在 WGCNA 网络分析之前，需要先过滤掉低方差的基因，这些基因往往没有特定的生物学功能，还会增加分析噪音
# ===============================

datExpr <- prepare_data_for_wgcna(dds = dds, quantile = 0.75, min_count = 5)
datExpr.chicken <- prepare_data_for_wgcna(dds = dds[, dds$Species == "Chicken"], quantile = 0.75, min_count = 5)
datExpr.duck <- prepare_data_for_wgcna(dds = dds[, dds$Species == "Duck"], quantile = 0.75, min_count = 5)


# ===============================
# 通过绘制样本树筛选样本
# ===============================

sampleTree <- hclust(dist(t(datExpr)), method = "average")
plot(sampleTree)

removesample <- c("MA-C-4", "BR-C-1")
datExpr <- datExpr[,!colnames(datExpr) %in% removesample] 
datExpr.chicken <- datExpr.chicken[, !colnames(datExpr.chicken) %in% "BR-C-1"]
datExpr.duck <- datExpr.duck[, !colnames(datExpr.duck) %in% "MA-C-4"]

### 查看三个数据集共享的基因
VennDiagram::venn.diagram(
  list(
    Combine = rownames(datExpr),
    Chicken = rownames(datExpr.chicken),
    Duck = rownames(datExpr.duck)
  ), filename = "~/Downloads/test.tiff",
  fill = c( "green", "gold", "darkorchid1")
)

VennDiagram::venn.diagram(
  list(
    Chicken = rownames(datExpr.chicken),
    Duck = rownames(datExpr.duck)
  ), filename = "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/质控后基因数量.tiff",
  fill = c( "green", "gold")
)


# ===============================
# 选择合适的软阈值
# ===============================

# 有不合适的样本时，没有自动选到软阈值，删除不合适的阈值后就得到了合适的软阈值

softPower <- auto_select_soft_Power(datExpr = t(datExpr))
softPower.chicken <- auto_select_soft_Power(datExpr = t(datExpr.chicken))
softPower.duck <- auto_select_soft_Power(datExpr = t(datExpr.duck))

softPower.chicken <- 10  #自动检测时选择14，有点过高了，因此稍微降低一点

# 软阈值差异 ≤ 1   →  统一使用同一个软阈值 ✅
# 软阈值差异 2-3   →  使用平均值作为折中 ⚖️  
# 软阈值差异 ≥ 4   →  各自使用最优软阈值 🔬

# ===============================
# 构建网络，并绘制模块图
# 【注意】WGCNA网络构建过程和DESeq2包存在冲突，加载后会导致程序运行错误
# ===============================

net.combine <- net_and_plot(datExpr = datExpr, softPower = softPower, labeltext = "Combine")
net.chicken <- net_and_plot(datExpr = datExpr.chicken, softPower = softPower.chicken, labeltext = "Chicken",
                            pdffile = "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Chicken.模块.pdf")
net.duck <- net_and_plot(datExpr = datExpr.duck, softPower = softPower.duck, labeltext = "Duck",
                         pdffile = "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Duck.模块.pdf" )

# ===============================
# 模块的相关性
# ===============================

MEs.chicken <- moduleEigengenes(t(datExpr.chicken), colors = labels2colors(net.chicken$colors))$eigengenes %>% orderMEs()
colnames(MEs.chicken) <- colnames(MEs.chicken) %>% gsub("ME","", .)
pdf("~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Chicken.模块相关性.pdf", width = 4, height = 4)
pheatmap::pheatmap(mat = cor(MEs.chicken))
dev.off()

MEs.duck <- moduleEigengenes(t(datExpr.duck), colors = labels2colors(net.duck$colors))$eigengenes %>% orderMEs()
colnames(MEs.duck) <- colnames(MEs.duck) %>% gsub("ME","", .)
pdf("~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Duck.模块相关性.pdf", width = 4, height = 4)
pheatmap::pheatmap(mat = cor(MEs.duck))
dev.off()

# ===============================
#  Module preservation analysis
# ===============================


multiData <- list(
  Chicken = list(data = t(datExpr.chicken)), # Transpose to make samples columns
  Duck = list(data = t(datExpr.duck)) # Transpose to make samples columns
)

multiColor <- list(
  Chicken = labels2colors(net.chicken$colors) %>% setNames(net.chicken$colors),
  Duck = labels2colors(net.duck$colors) %>% setNames(net.duck$colors)
)

preservation_chicken_vs_duck <- WGCNA::modulePreservation(
  multiData = multiData,
  multiColor = multiColor,
  referenceNetworks = 1,  # 将第一个数据集作为 reference
  nPermutations = 100, 
  randomSeed = 10086,
  parallelCalculation = TRUE,
  verbose = 3
)


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
ggsave(filename = "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Chicken.保守性分析.pdf", plot = p.pre.chicken, width = 6, height = 4)

# turquoise 和 brown 模块的复现性很高，查看一下是否是功能保守的模块


preservation_stats.duck <- preservation_duck_vs_chicken$preservation$Z$ref.Duck$inColumnsAlsoPresentIn.Chicken %>%
  mutate(Module = rownames(.)) %>%
  dplyr::filter(Module != "gold") 

p.pre.duck <- ggplot(data = preservation_stats.duck, aes(x = Module, y = Zsummary.pres, fill = Module))  + geom_bar(stat = "identity") +
  scale_fill_manual(values = preservation_stats.duck$Module) + ggtitle("Duck as ref")
ggsave(filename = "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Duck.保守性分析.pdf", plot = p.pre.duck, width = 6, height = 4)


### 用列表储存模块基因

genes_in_module.chicken <- split(net.chicken$colors %>% names, f = labels2colors(net.chicken$colors))
genes_in_module.duck <- split(net.duck$colors %>% names, f = labels2colors(net.duck$colors))


# ===============================
#  GO 和 KEGG 注释
# ===============================

# 各模块 GO 注释
enrichment_GO.chicken <- list()
enrichment_KEGG.chicken <- list()
gtf.chicken <- "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/data/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf"
emapperannotations.chicken <- "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/data/chicken.emapper.annotations"

for (module in genes_in_module.chicken %>% names()) {
  
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
gtf.duck <- "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/data/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.gtf"
emapperannotations.duck <- "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/data/duck.emapper.annotations"
for (module in genes_in_module.duck %>% names()) {
  
  cat("模块 ", module, "\n")
  
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
  cat("模块 ", module, "\n")
  p.GO <- anno_plot(object = enrichment_GO.chicken[[module]], title = paste0("GO: Chicken module ", module))
  p.KEGG <- anno_plot(object = enrichment_KEGG.chicken[[module]], title = paste0("KEGG: Chicken module ", module))
  p <- ggpubr::ggarrange(p.GO, p.KEGG, labels = c("A", "B"), ncol = 2)
  ggsave(filename = paste0("~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Chicken.",module,".anno.pdf" ), plot = p, width = 14, height = 6)
}

for (module in genes_in_module.duck %>% names()) {
  cat("模块 ",module, "\n")
  p.GO <- anno_plot(object = enrichment_GO.duck[[module]], title = paste0("GO: Duck module ", module))
  p.KEGG <- anno_plot(object = enrichment_KEGG.duck[[module]], title = paste0("KEGG: Duck module ", module))
  p <- ggpubr::ggarrange(p.GO, p.KEGG, labels = c("A", "B"), ncol = 2)
  ggsave(filename = paste0("~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Duck.",module,".anno.pdf" ), plot = p, width = 14, height = 6)
}


# ==============================================================================
#  Expression Profiles
# ==============================================================================

for (module in genes_in_module.chicken %>% names()) {
  
  submod_df <- datExpr.chicken[genes_in_module.chicken[[module]],] %>% 
    as.data.frame() %>% 
    mutate(gene_id = rownames(.)) %>%
    reshape2::melt() 
  
  p <- ggplot(submod_df, aes(x=variable, y=value, group = gene_id)) +
    geom_line(alpha = 0.3) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90)
    ) + 
    labs(x = "treatment",
         y = "normalized expression") +
    ggtitle(paste0(module, " module expression profiles in Chicken"))
  
  filename <- paste0("~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Chicken.",module, ".expressionprofiles.pdf" )
  ggsave(filename = filename,  width = 14, height = 4)
  
}


for (module in genes_in_module.duck %>% names()) {
  
  submod_df <- datExpr.duck[genes_in_module.duck[[module]],] %>% 
    as.data.frame() %>% 
    mutate(gene_id = rownames(.)) %>%
    reshape2::melt() 
  
  p <- ggplot(submod_df, aes(x=variable, y=value, group = gene_id)) +
    geom_line(alpha = 0.3) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90)
    ) + 
    labs(x = "treatment",
         y = "normalized expression") +
    ggtitle(paste0(module, " module expression profiles in Duck"))
  
  filename <- paste0("~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Duck.",module, ".expressionprofiles.pdf" )
  ggsave(filename = filename, plot = p,  width = 14, height = 4)
  
}


VennDiagram::venn.diagram(
  list(
    Chicken.Brown = genes_in_module.chicken$brown,
    Duck.Red = genes_in_module.duck$red,
    Duck.Pink = genes_in_module.duck$pink
  ), filename = "~/Downloads/test.tiff",
  fill = c( "brown", "red", "pink")
)



submod_df <- DESeq2::getVarianceStabilizedData(dds)[genes_in_module.chicken$green, colnames(datExpr.duck)]  %>% 
  as.data.frame() %>% 
  mutate(gene_id = rownames(.)) %>%
  reshape2::melt() 

p <- ggplot(submod_df, aes(x=variable, y=value, group = gene_id)) +
  geom_line(alpha = 0.3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) + 
  labs(x = "treatment",
       y = "normalized expression") +
  ggtitle("鸡细胞周期基因在鸭中的表达情况")

ggsave(filename = "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/鸡细胞周期基因在鸭中的表达情况.pdf", plot = p,  width = 14, height = 4)




# ==============================================================================
#  查看一下鸡Brown 模块特有的基因在鸭中是什么情况
# ==============================================================================


ChickenUniueImmuneGenes <- setdiff( setdiff(genes_in_module.chicken$brown, genes_in_module.duck$red), genes_in_module.duck$pink)

ChickenUniueImmuneGenes

# 把鸡免疫模块相对鸭特有的基因找出来，并进行功能注释，查看这些基因是什么功能
ChickenUniueImmuneGenes_InDuck <- list()
for (module in genes_in_module.duck %>% names()) {
  ChickenUniueImmuneGenes_InDuck[[module]]  <- genes_in_module.duck[[module]][genes_in_module.duck[[module]] %in% ChickenUniueImmuneGenes]
}

enrichment_GO.test <- GO_enrichment(
  diffgenes = ChickenUniueImmuneGenes_InDuck$brown %>% extract_gene_from_gene_vs_gene_pair(side = 2), 
  backgroundgenes = extract_gene_from_gene_vs_gene_pair(SummarizedExperiment::assay(dds) %>% rownames(), side = 2),
  GTF = gtf.duck, emapperannotations = emapperannotations.duck
)

anno_plot(object = enrichment_GO.test, title = "unique gene anno by duck")

enrichment_GO.test@result -> test


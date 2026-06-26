
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

# datExpr <- prepare_data_for_wgcna(dds = dds, quantile = 0.6, min_count = 5)
datExpr.chicken <- prepare_data_for_wgcna(dds = dds[, dds$Species == "Chicken"], quantile = 0.6, min_count = 5)
datExpr.duck <- prepare_data_for_wgcna(dds = dds[, dds$Species == "Duck"], quantile = 0.6, min_count = 5)


# ===============================
# 通过绘制样本树筛选样本
# ===============================

# sampleTree <- hclust(dist(t(datExpr)), method = "average")
# plot(sampleTree)
# 
# removesample <- c("MA-C-4","MA-C-1","MA-C-2", "BR-C-1",  "BR-C-2")
# datExpr <- datExpr[,!colnames(datExpr) %in% removesample] 
datExpr.chicken <- datExpr.chicken[, !colnames(datExpr.chicken) %in% c("BR-C-1", "BR-C-2")]
datExpr.duck <- datExpr.duck[, !colnames(datExpr.duck) %in% c("MA-C-4", "MA-C-1","MA-C-2")] #, "MA-C-32", "MA-C-33", "MA-C-36"

datExpr.duck <- datExpr.duck[, !colnames(datExpr.duck) %in% c("MA-C-4", "MA-C-1","MA-C-2","MA-C-32", "MA-C-33", "MA-C-34", "MA-C-35", "MA-C-36","PK-C-1", "PK-C-2", "PK-C-3")] 


### 查看三个数据集共享的基因
# VennDiagram::venn.diagram(
#   list(
#     Combine = rownames(datExpr),
#     Chicken = rownames(datExpr.chicken),
#     Duck = rownames(datExpr.duck)
#   ), filename = "~/Downloads/test.tiff",
#   fill = c( "green", "gold", "darkorchid1")
# )
# 
# VennDiagram::venn.diagram(
#   list(
#     Chicken = rownames(datExpr.chicken),
#     Duck = rownames(datExpr.duck)
#   ), filename = "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/质控后基因数量.tiff",
#   fill = c( "green", "gold")
# )


# ===============================
# 选择合适的软阈值
# ===============================

# 有不合适的样本时，没有自动选到软阈值，删除不合适的阈值后就得到了合适的软阈值

# softPower <- auto_select_soft_Power(datExpr = t(datExpr))
softPower.chicken <- auto_select_soft_Power(datExpr = t(datExpr.chicken))
softPower.duck <- auto_select_soft_Power(datExpr = t(datExpr.duck))

# softPower.chicken <- 10  #自动检测时选择14，有点过高了，因此稍微降低一点
# softPower.duck <- 10
# 软阈值差异 ≤ 1   →  统一使用同一个软阈值 ✅
# 软阈值差异 2-3   →  使用平均值作为折中 ⚖️  
# 软阈值差异 ≥ 4   →  各自使用最优软阈值 🔬

# ==============================================================================
# 构建网络，并绘制模块图
# 【注意】WGCNA网络构建过程和DESeq2包存在冲突，加载后会导致程序运行错误
# ==============================================================================

# net.combine <- net_and_plot(datExpr = datExpr, softPower = softPower, labeltext = "Combine")
net.chicken <- net_and_plot(datExpr = datExpr.chicken, softPower = softPower.chicken, labeltext = "Chicken",
                            pdffile = "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Chicken.模块.pdf")
net.duck <- net_and_plot(datExpr = datExpr.duck, softPower = softPower.duck, labeltext = "Duck",
                         pdffile = "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Duck.模块.pdf" )

net.chicken$colors %>% table()
net.duck$colors %>% table()


# mergeCutHeight = 0.25 意味着相关性 > 0.75 的模块会被合并
# mergeCutHeight = 0.15 意味着相关性 > 0.85 的模块才会被合并


# ==============================================================================
### 用列表储存模块基因
# ==============================================================================

genes_in_module.chicken <- split(net.chicken$colors %>% names, f = net.chicken$colors)
genes_in_module.duck <- split(net.duck$colors %>% names, f = net.duck$colors)


# ==============================================================================
# 模块的相关性
# ==============================================================================

MEs.chicken <- moduleEigengenes(t(datExpr.chicken), colors = net.chicken$colors)$eigengenes %>% orderMEs()
# colnames(MEs.chicken) <- colnames(MEs.chicken) %>% gsub("ME","", .)
pdf("~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Chicken.模块相关性.pdf", width = 8, height = 8)
pheatmap::pheatmap(mat = cor(MEs.chicken), display_numbers = TRUE, main = "Chicken")
dev.off()

MEs.duck <- moduleEigengenes(t(datExpr.duck), colors = net.duck$colors)$eigengenes %>% orderMEs()
# colnames(MEs.duck) <- colnames(MEs.duck) %>% gsub("ME","", .)
pdf("~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Duck.模块相关性.pdf", width = 8, height = 8)
pheatmap::pheatmap(mat = cor(MEs.duck), display_numbers = TRUE, main = "Duck")
dev.off()


# ==============================================================================
# 评估模块和网络的质量
# ==============================================================================

### 统计特征向量对模块的解释度

datKME.chicken <- signedKME(datExpr = t(datExpr.chicken), datME = MEs.chicken)
datKME.duck <- signedKME(datExpr = t(datExpr.duck), datME = MEs.duck)
calculate_variance_explained(datKME = datKME.chicken, genes_in_module = genes_in_module.chicken)
calculate_variance_explained(datKME = datKME.duck, genes_in_module = genes_in_module.duck)


# 好的模块：VE > 50%
# 中等模块：VE 30-50%
# 差的模块：VE < 30%

### 模块成员度 - 每个基因与模块特征值的相关性

# 质量标准
# 优秀模块：>80%基因 kME > 0.7
# 良好模块：>60%基因 kME > 0.7
# 较差模块：<50%基因 kME > 0.7

for (module in genes_in_module.chicken %>% names()) {
  hist(datKME.chicken[genes_in_module.chicken[[module]], paste0("kME", module)], 
       main=paste("kME distribution in", module, "module"),
       xlab="Module Membership (kME)",
       col="lightblue")
  abline(v=0.7, col="red", lwd=2, lty=2)
  percent <- sum((datKME.chicken[genes_in_module.chicken[[module]], paste0("kME", module)] > 0.7 )) / length(genes_in_module.chicken[[module]])
  cat("模块 ",module," 有百分之 ", percent, " 的基因 kME > 0.7\n")
}


for (module in genes_in_module.duck %>% names()) {
  hist(datKME.duck[genes_in_module.duck[[module]], paste0("kME", module)], 
       main=paste("kME distribution in", module, "module"),
       xlab="Module Membership (kME)",
       col="lightblue")
  abline(v=0.7, col="red", lwd=2, lty=2)
  percent <- sum((datKME.duck[genes_in_module.duck[[module]], paste0("kME", module)] > 0.7 )) / length(genes_in_module.duck[[module]])
  cat("模块 ",module," 有百分之 ", percent, " 的基因 kME > 0.7\n")
}


### 模块内连接度

# 好的模块特征：
# 1. 有明显的hub基因（高连接度）
# 2. 连接度分布呈幂律分布

# 可视化
# 计算模块内连接度

for (module in genes_in_module.chicken %>% names()) {
  IMC <- intramodularConnectivity.fromExpr(datExpr = t(datExpr.chicken), net.chicken$colors  )
  hist( IMC$kWithin[which(labels2colors(net.chicken$colors) == module)],
        main=paste("Connectivity in", module, "module"),
        xlab="Intramodular Connectivity")
}


### Bootstrap重抽样检验
# 
# nBootstrap <- 100
# moduleStability <- matrix(0, nrow=nBootstrap, ncol=length(genes_in_module.chicken %>% names()))
# 
# for(i in 1:nBootstrap) {
#   # 重采样
#   bootSamples <- sample(1:ncol(t(datExpr.chicken)), replace=TRUE)
#   bootExpr <- t(datExpr.chicken)[, bootSamples]
#   
#   # 重新构建网络
#   bootNet <- blockwiseModules(bootExpr, 
#                               power=softPower.chicken,
#                               networkType="unsigned",
#                               minModuleSize=30)
#   
#   # 比较模块保守性
#   # ...
# }
# 






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
  parallelCalculation = FALSE,
  verbose = 3
)


preservation_duck_vs_chicken <- WGCNA::modulePreservation(
  multiData = multiData,
  multiColor = multiColor,
  referenceNetworks = 2,  # 将第一个数据集作为 reference
  nPermutations = 200, 
  randomSeed = 10086,
  parallelCalculation = FALSE,
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


# ===============================
#  GO 和 KEGG 注释
# ===============================

# 各模块 GO 注释
enrichment_GO.chicken <- list()
enrichment_KEGG.chicken <- list()
gtf.chicken <- "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/data/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf"
emapperannotations.chicken <- "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/data/chicken.emapper.annotations"

for (module in genes_in_module.chicken %>% names()) {
  
  cat("######### Chicken 模块 ", module, "\n")
  
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
  
  cat("######### Duck 模块 ", module, "\n")
  
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
  cat("######### Chicken 模块 ", module, "\n")
  tryCatch({
    p.GO <- anno_plot(object = enrichment_GO.chicken[[module]], title = paste0("GO: Chicken module ", module))
    p.KEGG <- anno_plot(object = enrichment_KEGG.chicken[[module]], title = paste0("KEGG: Chicken module ", module))
    p <- ggpubr::ggarrange(p.GO, p.KEGG, labels = c("A", "B"), ncol = 2)
    ggsave(filename = paste0("~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Chicken.",module,".anno.pdf" ), plot = p, width = 14, height = 6)}, error = function(e) {
      cat("模块", module, "处理时出错:", conditionMessage(e), "\n")
    })
}


for (module in genes_in_module.duck %>% names()) {
  cat("######### Duck 模块 ",module, "\n")
  tryCatch({
    p.GO <- anno_plot(object = enrichment_GO.duck[[module]], title = paste0("GO: Duck module ", module))
    p.KEGG <- anno_plot(object = enrichment_KEGG.duck[[module]], title = paste0("KEGG: Duck module ", module))
    p <- ggpubr::ggarrange(p.GO, p.KEGG, labels = c("A", "B"), ncol = 2)
    ggsave(filename = paste0("~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/Duck.",module,".anno.pdf" ), plot = p, width = 14, height = 6)
  }, error = function(e) {
    cat("模块", module, "处理时出错:", conditionMessage(e), "\n")
  })
}




# ==============================================================================
#  保存分析结果
# ==============================================================================


WGCNAfile <- list(
  "Chicken" = list(
    "softPower" =  softPower.chicken,
    # "preservation" = preservation_chicken_vs_duck,
    "genes" = genes_in_module.chicken,
    "enrichment_GO" = enrichment_GO.chicken,
    "enrichment_KEGG" = enrichment_KEGG.chicken
  ),
  "Duck" = list(
    "softPower" = softPower.duck,
    # "preservation" = preservation_duck_vs_chicken,
    "genes" = genes_in_module.duck,
    "enrichment_GO" = enrichment_GO.duck,
    "enrichment_KEGG" = enrichment_KEGG.duck
  )
)

saveRDS(WGCNAfile, file = "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/WGCNAfile.rds")


# ==============================================================================
#  查看指定模块基因在另一个物种中的哪些模块中
# ==============================================================================


find_gene_in_which_module <- function(targetgenes, net){
  targetgenes <- targetgenes
  net <- net
  
  targetgenes_in_anotherspecies <- net$colors[targetgenes] 
  targetgenes_in_anotherspecies <- targetgenes_in_anotherspecies[!targetgenes_in_anotherspecies %>% is.na()] 
  targetgenes_in_anotherspecies <- split(targetgenes_in_anotherspecies %>% names(), f = targetgenes_in_anotherspecies %>% unname())
  
  return(targetgenes_in_anotherspecies)
}


find_gene_in_which_module(targetgenes = genes_in_module.chicken$black, net = net.duck)

find_gene_in_which_module(targetgenes = genes_in_module.chicken$tan, net = net.duck)




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
    Duck.Yellow = genes_in_module.duck$yellow
  ), filename = "~/Desktop/04.湘湖实验室/家禽病毒课题/RNA-seq/WGCAN/鸡和鸭免疫模块基因共享韦恩图.tiff",
  fill = c( "brown", "yellow")
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



submod_df.red <- datExpr.duck[c(genes_in_module.duck$red),] %>% 
  as.data.frame() %>% 
  mutate(gene_id = rownames(.)) %>%
  reshape2::melt() %>% 
  mutate("module" = "red")

submod_df.pink <- datExpr.duck[c(genes_in_module.duck$pink),] %>% 
  as.data.frame() %>% 
  mutate(gene_id = rownames(.)) %>%
  reshape2::melt() %>% 
  mutate("module" = "pink")

submod_df <- rbind(submod_df.red, submod_df.pink)

ggplot(submod_df, aes(x=variable, y=value, group = gene_id, color = module)) +
  geom_smooth(se = FALSE) +
  # geom_line(alpha = 0.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) + 
  labs(x = "treatment",
       y = "normalized expression") +
  ggtitle("Red and Pink modules")


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
  diffgenes = ChickenUniueImmuneGenes_InDuck$turquoise %>% extract_gene_from_gene_vs_gene_pair(side = 2), 
  backgroundgenes = extract_gene_from_gene_vs_gene_pair(SummarizedExperiment::assay(dds) %>% rownames(), side = 2),
  GTF = gtf.duck, emapperannotations = emapperannotations.duck
)

anno_plot(object = enrichment_GO.test, title = "Turquoise")

enrichment_GO.test@result -> test


for (module in ChickenUniueImmuneGenes_InDuck %>% names) {
  genes <- paste0(ChickenUniueImmuneGenes_InDuck[[module]], collapse = ", ")
  cat(module, "\t", genes, "\n")
}







#setwd("/home/zhangchunyuan/zhangchunyuan/RNA-seq")
#.libPaths("/home/zhangchunyuan/tools/rstudio/lib/R/library")

#source("scripts/functions.R")


suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))

option_list <- list(
  make_option("--infotable", type="character", default=NULL, help="差异倍数"),
  make_option("--lfc", type="double", default=1, help="差异倍数"),
  make_option("--pval", type="double", default=0.05, help="显著性阈值，转录组一般使用padj"),
  make_option("--untreated", type="character", default=NULL, help="指定对照组"),
  make_option("--CountingMethod", type="character", default="rsem", help="表达量统计工具，影响数据导入的格式"),
  make_option("--CountingFiles.isoforms", type="character", default=NULL, help="当使用一步法时使用该参数，基于转录本的表达量统计文件"),
  make_option("--orthologgenes", type="character", default=NULL, help="1:1的同源基因"),
  make_option("--outdir", type = "character", default = "result/DEseq/", help = "差异分析结果存储路径")
)

args <- parse_args(OptionParser(option_list=option_list))

infotable <- args$infotable
lfc <- args$lfc
pval <- args$pval
out <- args$out
untreated <- args$untreated
CountingFiles.isoforms <- args$CountingFiles.isoforms
CountingFiles.isoforms <- CountingFiles.isoforms %>% strsplit(",") %>% unlist()
CountingMethod <- args$CountingMethod
orthologgenes <- args$orthologgenes
outdir <- args$outdir


# lfc = 1
# pval = 0.05
# untreated <- "Chicken"
# infotable <- "/home/zhangchunyuan/zhangchunyuan/RNA-seq/infotable.csv"
# CountingMethod <- "rsem"
# outdir <- "result/DEseq"
# CountingFiles.isoforms <- list.files(path = "result/RSEM", pattern = ".isoforms.results", full.names = T)
# orthologgenes <- "/home/zhangchunyuan/zhangchunyuan/RNA-seq/result/Ortho_chicken_vs_duck/one_to_one_orthologgenes.txt"


#====================================================================
# 读取元信息
#====================================================================

meta_table <- read.table(infotable, sep = ",", header = T, row.names = 1)
meta_table_list <- split(x = meta_table, f = ~Species)
species <- names(meta_table_list)

treated <- setdiff(meta_table$Species %>% as.vector() %>% unique(), untreated)


#====================================================================
# 读取1:1的同源基因，并给名称不同的同源基因重命名，用于后续提取基因
#====================================================================

orthologgenes <- read.delim(file = orthologgenes)
orthologgenes <- orthologgenes %>% mutate(
  genename =  ifelse(!!as.name(species[1]) == !!as.name(species[2]),!!as.name(species[1])  ,  paste0(!!as.name(species[1]), "_vs_", !!as.name(species[2])))
) 


#====================================================================
# 按照meta_table的顺序读取文件，注意处理同源基因
#====================================================================

counts.genes <- list()
for (s in names(meta_table_list)) {
  
  sub_table <- meta_table_list[[s]]
  pattern <- rownames(sub_table) %>% paste0(collapse = "|")
  sub_isoforms <- CountingFiles.isoforms[grep(pattern = pattern, x = CountingFiles.isoforms)] # 按照meta_table的顺序进行调整
  
  tx2gene <- read.delim(file = sub_isoforms[1], header = T) %>% select(transcript_id, gene_id)
  temp <- tximport::tximport(files = as.character(sub_isoforms), type = CountingMethod, tx2gene = tx2gene )
  
  for (arr in names(temp)[1:3]) {
    temp[[arr]] <- temp[[arr]][orthologgenes[,s],]
    rownames(temp[[arr]]) <- orthologgenes[,"genename"]
    counts.genes[[arr]] <- cbind(counts.genes[[arr]], temp[[arr]])
  }
  
  counts.genes[[names(temp)[4]]] <- temp[[names(temp)[4]]]
}


#====================================================================
# 提取各物种同源基因的子集，给不同名字的基因名用新命名的基因重命名，好用于后续差异分析
#====================================================================

meta_table$Species <- as.factor(meta_table$Species)

###创建DESeq2对象，指定谁是control，运行DEseq
dds <- DESeq2::DESeqDataSetFromTximport(counts.genes, colData = meta_table, design = ~Species)
dds$Species <- relevel(dds$Species, ref = untreated)
dds <- DESeq2::DESeq(dds)  # 这一步自动进行了标准化


#====================================================================
# 质量检查，如果样本量大于3，并且没有明显的技术性差异，则可以直接使用DESeq(dds)输出的结果
#====================================================================
### 检查测序深度是否均匀
pdf(file = file.path(outdir, "barplot.pdf"), width = 6, height = 5)
barplot(colSums(counts(dds)), names.arg = colnames(dds), las = 2)
dev.off()

### PCA图 - 样本按预期分组
vsd <- vst(dds, blind = FALSE)
p.PCA <- plotPCA(vsd, intgroup = "Species")
ggsave(filename = file.path(outdir, "PCA.pdf"), plot = p.PCA, width = 6, height = 5)

### 样本热图，看组内相关性
library(pheatmap)
cor_matrix <- cor(assay(vsd))
pdf(file = file.path(outdir, "cor_samples.pdf"), width = 6, height = 5)
pheatmap(cor_matrix, annotation_col = as.data.frame(colData(dds)["Species"]))
dev.off()

### 输出Disp图
pdf(file = file.path(outdir, "Disp.pdf"),  width = 6, height = 5)
DESeq2::plotDispEsts(dds)
dev.off()


#====================================================================
# 绘制火山图，输出结果
#====================================================================


diffgenes <- c()

for (t in treated) {
  
  ### 从dds文件中提取结果 ===============================
  pairname <- paste0("Species","_", t, "_vs_", untreated)
  result <- DESeq2::results(dds, name = pairname)
  
  ###绘制火山图 ========================================
  p.vol <- EnhancedVolcano::EnhancedVolcano(result,
                                            lab = rownames(result),
                                            x = 'log2FoldChange',
                                            y = 'padj',
                                            title = pairname,
                                            subtitle = NULL,
                                            pCutoff = pval,
                                            FCcutoff = lfc,
                                            labSize = 3,
                                            gridlines.major = FALSE,
                                            gridlines.minor = FALSE,
                                            legendPosition = "right"
  )
  
  
  filename <- paste0(pairname, ".Volcano.pdf")
  
  ggsave(filename = file.path(outdir,filename), plot = p.vol, width = 10, height = 7, create.dir = TRUE)
  
  
  ### 输出差异分析表格 ==================================
  result <- result[order(result$padj),]
  result_dataframe <- as.data.frame(result) %>% 
    tibble::rownames_to_column(var = "genename") %>%
    mutate(
      Up = case_when(
        padj < pval & log2FoldChange > lfc ~ t,
        padj < pval & log2FoldChange < -lfc ~ untreated,
        TRUE ~ "NotSig"
      )
    ) %>%
  left_join(., orthologgenes, by = "genename")
  
  write.csv(result_dataframe, file.path(outdir,paste0(pairname, "_DEseq2_results.csv")),
    quote = F,  row.names = F, )
  

  ### 将差异基因都整理到一个表格中，便于控制输入和输出文件========
  tmpdata <- result_dataframe %>% filter(Up != "NotSig") %>% mutate(compare = pairname)
  diffgenes <- rbind(diffgenes,tmpdata )
  remove(tmpdata)
  
  
  ### 绘制热图 ==========================================
  sig_genes <- subset(result, padj < pval & abs(log2FoldChange) > lfc) %>% rownames()
  groupA <- meta_table %>% filter(Species == t) %>% rownames()
  groupB <- meta_table %>% filter(Species == untreated) %>% rownames() 
  #pheatmap( assay(dds)[sig_genes,c(groupA, groupB)], scale = "row" )
  p.heatmap <- pheatmap::pheatmap(assay(vst(dds))[sig_genes,c(groupA, groupB)], scale = "row", show_rownames = F)
  filename <- paste0(pairname, ".Heatmap.pdf")
  ggsave(plot = p.heatmap, filename = file.path(outdir,filename), width = 6, height = 6, create.dir = TRUE )
  

}


write.csv(x = diffgenes, file = file.path(outdir, "diffgenes.csv"), 
             quote = F, row.names = F, )



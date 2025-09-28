#setwd("/home/zhangchunyuan/zhangchunyuan/lipan/RNA-seq")
#.libPaths("/home/zhangchunyuan/tools/rstudio/lib/R/library")

# ====================================================================
# suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
suppressMessages(library(tximport))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
# suppressMessages(library(pheatmap))
# suppressMessages(library(apeglm))
# suppressMessages(library(ggrepel))
# suppressMessages(library(EnhancedVolcano))
suppressMessages(library(optparse))


# ====================================================================
#  命令行参数传递
# ====================================================================

option_list <- list(
  make_option("--infotable", type="character", default=NULL, help="差异倍数"),
  make_option("--lfc", type="double", default=1, help="差异倍数"),
  make_option("--pval", type="double", default=0.05, help="显著性阈值，转录组一般使用padj"),
  make_option("--untreated", type="character", default=NULL, help="指定对照组"),
  make_option("--CountingMethod", type="character", default="rsem", help="表达量统计工具，影响数据导入的格式"),
  make_option("--CountingFiles.isoforms", type="character", default=NULL, help="基于转录本的表达量统计文件"),
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
outdir <- args$outdir


#Set your log-fold-change and p-value thresholds
# lfc = 1
# pval = 0.05
# untreated <- "mKO"
# infotable <- "/home/zhangchunyuan/zhangchunyuan/lipan/RNA-seq/infotable.csv"
# CountingMethod <- "rsem"
# CountingFiles.isoforms <- c("result/RSEM/576-F4.isoforms.results", "result/RSEM/576-F9.isoforms.results",  "result/RSEM/576-W5.isoforms.results",  "result/RSEM/589-1-GF.isoforms.results",  "result/RSEM/589-4-GF.isoforms.results",  "result/RSEM/W2-576.isoforms.results")
outdir <- "result/DEseq"


# ====================================================================
#   读取元信息
# ====================================================================

meta_table <- read.table(infotable, sep = ",", header = T, row.names = 1)
meta_table$GroupID <- as.factor(meta_table$GroupID)
treated <- setdiff(meta_table$GroupID %>% as.vector() %>% unique(), untreated)


### 将比对到转录组的表达数据汇总到基因水平
pattern <- rownames(meta_table) %>% paste0(collapse = "|") 
CountingFiles.isoforms <- CountingFiles.isoforms[grep(pattern = pattern, x = CountingFiles.isoforms)]

tx2gene <- read.delim(file = CountingFiles.isoforms[1], header = T) %>% select(transcript_id, gene_id)
counts.genes <- tximport::tximport(files = as.character(CountingFiles.isoforms), type = CountingMethod, tx2gene = tx2gene )

###创建DESeq2对象 
dds <- DESeq2::DESeqDataSetFromTximport(counts.genes, colData = meta_table, design = ~GroupID )

# ====================================================================
#  DESeq2分析
# ====================================================================

# 制定比较组，将mKO指定为case
# 不指定的话会根据字符顺序比较
dds$GroupID <- relevel(dds$GroupID, ref = untreated)

# 运行DESeq2
dds <- DESeq(dds)  # 这一步自动进行了标准化


### 检查测序深度是否均匀
pdf(file = file.path(outdir, "barplot.pdf"), width = 6, height = 5)
barplot(colSums(counts(dds)), names.arg = colnames(dds), las = 2)
dev.off()

### PCA图 - 样本按预期分组
vsd <- vst(dds, blind = FALSE)
p.PCA <- plotPCA(vsd, intgroup = "GroupID")
ggsave(filename = file.path(outdir, "PCA.pdf"), plot = p.PCA, width = 6, height = 5)


### 样本热图，看组内相关性
cor_matrix <- cor(assay(vsd))
pdf(file = file.path(outdir, "cor_samples.pdf"), width = 6, height = 5)
pheatmap::pheatmap(cor_matrix, annotation_col = as.data.frame(colData(dds)["GroupID"]))
dev.off()


# 输出Disp图
pdf(file = file.path(outdir, "Disp.pdf"),  width = 6, height = 5)
DESeq2::plotDispEsts(dds)
dev.off()
# resultsNames(dds)

# ====================================================================
# PCA
# ====================================================================
### 一般分析只使用top的差异基因，但是我们这里使用所有的基因数据进行分析

### 样本量比较少的时候和比较多的时候所使用的方法不同
### 样本量比较少时使用以下方法
# rld <- DESeq2::rlog(dds, blind = TRUE) #This function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size.
# rld_mat <- assay(rld)
# pca <- prcomp(t(rld_mat))
# pca <- pca$x %>% as.data.frame() 
# pca$GroupID <- meta_table[rownames(pca),"GroupID"]
# p.pca <- ggplot(data = pca, aes(x = PC1, y = PC2, colour = GroupID)) + geom_point()
# ggsave(filename = file.path(outdir, "PCA.pdf"), plot = p.pca, width = 6, height = 5, create.dir = TRUE)


# ====================================================================
#  每个比较对信息提取
# ====================================================================

diffgenes <- c()

for (t in treated) {
  
  ### 从dds文件中提取结果 ===============================
  pairname <- paste0("GroupID","_", t, "_vs_", untreated)
  
  ### 对不可靠的倍数变化估计进行贝叶斯收缩 
  ### 当样本量较少时这种收缩是很有必要的
  min_sample_size <- meta_table$GroupID %>% as.vector() %>% table() %>% min()
  if(min_sample_size < 3){
    # 当样本量小于3时，使用log2 fold-change进行缩放
    result <- DESeq2::lfcShrink(dds, coef = pairname, type = "apeglm") # coef 指定了比较组
    ### The contens of the LFC dataframe contain the log2 fold-change, as well as the p-value and adjusted p-value
  }else{
    result <- results(dds, name = pairname)
  }
  
  
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
  write.csv(as.data.frame(result), file.path(outdir,paste0(pairname, "_DEseq2_results.csv")))
  
  
  ### 绘制热图 ==========================================
  sig_genes <- subset(result, padj < pval & abs(log2FoldChange) > lfc) %>% rownames()
  groupA <- meta_table %>% filter(GroupID == t) %>% rownames()
  groupB <- meta_table %>% filter(GroupID == untreated) %>% rownames() 
  #pheatmap( assay(dds)[sig_genes,c(groupA, groupB)], scale = "row" )
  p.heatmap <- pheatmap::pheatmap(assay(vst(dds))[sig_genes,c(groupA, groupB)], scale = "row", show_rownames = F)
  filename <- paste0(pairname, ".Heatmap.pdf")
  ggsave(plot = p.heatmap, filename = file.path(outdir,filename), width = 6, height = 6, create.dir = TRUE )

  
  ### 将差异基因都整理到一个表格中，便于控制输入和输出文件
  diffgenes <- rbind(diffgenes, result %>% data.frame %>% filter(abs(log2FoldChange) > pval, padj < pval) %>% mutate(compare = pairname))

}

diffgenes <- diffgenes %>% tibble::rownames_to_column(var = "genename")
write.table(x = diffgenes, file = file.path(outdir, "diffgenes.tsv"), 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)





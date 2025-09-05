setwd("~/Desktop/04.湘湖实验室/李攀RNA-seq/")
library(rtracklayer)
library(dplyr)
library(tximport)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(apeglm)
library(VennDiagram)



# ====================================================================
#  参数设置
# ====================================================================


#Set your log-fold-change and p-value thresholds
lfc = 1
pval = 0.05
untreated <- "mKO"


# ====================================================================
#  导入数据
# ====================================================================

### 读取元数据 ==================================
meta_table <- read.table("~/Desktop/04.湘湖实验室/李攀RNA-seq/metafile", sep = ",", header = T, row.names = 1)
meta_table$GroupID <- as.factor(meta_table$GroupID)


treated <- setdiff(meta_table$GroupID %>% as.vector() %>% unique(), untreated)

### 从gtf文件中提取gene_id和transcript_id信息 =====
gtf <- import("~/Desktop/04.湘湖实验室/李攀RNA-seq/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gtf")
tx2gene <- gtf %>%
  as.data.frame() %>%
  select(transcript_id,gene_id) %>% # 顺序必须是先转录本ID，再基因ID
  filter(!is.na(transcript_id)) %>%
  unique()

genenamemapfile <- gtf %>% 
  as.data.frame() %>% 
  select(gene_id, gene_name ) %>% 
  unique() 
rownames(genenamemapfile) <- genenamemapfile$gene_id


### 读取RNA表达数量信息 ==============================
samples = read.table("~/Desktop/04.湘湖实验室/李攀RNA-seq/samples.txt")

### 导入转录本水平数据汇总到基因水平
counts.imported <- tximport::tximport(files = as.character(samples[,2]), type = 'salmon', tx2gene = tx2gene)

###创建DESeq2对象 =====================================
dds <- DESeqDataSetFromTximport(counts.imported, colData = meta_table, design = ~GroupID )



# ====================================================================
#  DESeq2分析
# ====================================================================

# 制定比较组，将mKO指定为case
# 不指定的话会根据字符顺序比较
dds$GroupID <- relevel(dds$GroupID, ref = untreated)


# 运行DESeq2
dds <- DESeq(dds)  # 这一步自动进行了标准化

# plotDispEsts(dds)
# resultsNames(dds)

# ====================================================================
# PCA
# ====================================================================

### 样本量比较少的时候和比较多的时候所使用的方法不同

### 样本量比较少时使用以下方法
rld <- DESeq2::rlog(dds, blind = TRUE) #This function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size.
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))
pca <- pca$x %>% as.data.frame() 
pca$GroupID <- meta_table[rownames(pca),"GroupID"]
p <- ggplot(data = pca, aes(x = PC1, y = PC2, colour = GroupID)) + geom_point()
ggsave(filename = "PCA.pdf", plot = p, width = 6, height = 5)



# ====================================================================
#  每个比较对信息提取
# ====================================================================


for (t in treated) {
  
  ### 从dds文件中提取结果 ===============================
  pairname <- paste0("GroupID","_", t, "_vs_", untreated)
  res <-  results(dds,  contrast = c("GroupID", t, untreated))   

  ### 对不可靠的倍数变化估计进行贝叶斯收缩 ==============
  ### 当样本量较少时这种收缩是很有必要的
  LFC <- DESeq2::lfcShrink(dds, coef = pairname, type = "apeglm") %>% as.data.frame()
  LFC$gene_name <- genenamemapfile[rownames(LFC), "gene_name"]

  
  ###绘制火山图 ========================================
  p <- EnhancedVolcano(LFC,
                  lab = LFC$gene_name,
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  title = pairname,
                  pCutoff = pval,
                  FCcutoff = lfc,
                  labSize = 3,
                  gridlines.major = FALSE,
                  gridlines.minor = FALSE,
                  legendPosition = "right"
                  )
  filename <- paste0(pairname, ".Volcano.pdf")
  ggsave(filename = filename, plot = p, width = 10, height = 7)
  
  
  ### 输出差异分析表格 ==================================
  LFC <- LFC[order(LFC$pvalue),]
  write.csv(as.data.frame(LFC), paste0(pairname, "_DEseq2_results.csv"))
  
  
}













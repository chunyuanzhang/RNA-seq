setwd("~/Desktop/04.湘湖实验室/李攀RNA-seq/")
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
suppressMessages(library(tximport))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(apeglm))
suppressMessages(library(VennDiagram))
suppressMessages(library(ggrepel))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(optparse))
suppressMessages(library(org.Gg.eg.db))
library(clusterProfiler)
library(enrichplot)

# ====================================================================
#  参数设置
# ====================================================================


#Set your log-fold-change and p-value thresholds
lfc = 1
pval = 0.05
untreated <- "mKO"
metafile <- "~/Desktop/04.湘湖实验室/李攀RNA-seq/metafile"
gtf <- "~/Desktop/04.湘湖实验室/李攀RNA-seq/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gtf"
metafile <- "~/Desktop/04.湘湖实验室/李攀RNA-seq/metafile"
mappedfiles <- c("result/576-F4/quant.sf",  "result/576-F9/quant.sf","result/576-W5/quant.sf", "result/W2-576/quant.sf", "result/589-1-GF/quant.sf", "result/589-4-GF/quant.sf" )

# ====================================================================
#  命令行参数传递
# ====================================================================

option_list <- list(
  make_option("--metafile", type="character", default=NULL, help="差异倍数"),
  make_option("--species", type="character", default=NULL, help="物种，与注释时使用的背景数据有关系"),
  make_option("--lfc", type="double", default=1, help="差异倍数"),
  make_option("--pval", type="double", default=0.05, help="p值阈值"),
  make_option("--gtf", type="character", default=NULL, help="gtf文件"),
  make_option("--untreated", type="character", default=NULL, help="指定对照组"),
  make_option("--mappedfiles", type="character", default=NULL, help="RNAseq数据比对文件")
)

args <- parse_args(OptionParser(option_list=option_list))

metafile <- args$metafile
lfc <- args$lfc
pval <- args$pval
gtf <- args$gtf
out <- args$out
untreated <- args$untreated
mappedfiles <- args$mappedfiles
mappedfiles <- mappedfiles %>% strsplit(",") %>% unlist()
species <- args$species

if(stringr::str_to_lower(species) == "chicken"){
  backgroupDB = org.Gg.eg.db
  organism = "gga"
}

if(stringr::str_to_lower(species) == "duck"){
  backgroupDB = 
  organism 
}


# ====================================================================
#  函数
# ====================================================================

GO_enrichment <- function(diff_genes, backgroupDB){
  
  diff_genes <- diff_genes
  backgroupDB <- backgroupDB
  
  ego <- enrichGO(gene = diff_genes,
                  OrgDb = backgroupDB,  # org.Gg.eg.db 是鸡的数据包
                  keyType = "SYMBOL",  # 直接指定ENSEMBL类型
                  ont = "ALL",         # 三种过程同时计算
                  pvalueCutoff = 1,    # 输出所有计算结果
                  qvalueCutoff = 1     # 输出所有计算结果
  )
  
  p.GO <- dotplot(ego, 
                  showCategory=10,      # 显示top10富集项
                  color="p.adjust",     # 颜色映射校正p值
                  size="Count",         # 点大小映射基因数量
                  split="ONTOLOGY") +   # 按GO类别分组
    facet_grid(ONTOLOGY~., scale="free") + # 自由缩放y轴
    scale_color_gradient(low="red", high="blue") # 自定义颜色
  
  
  return(list(ego, p.GO))
  
}

KEGG_enrchment <- function(diff_genes, backgroupDB){
  
  diff_genes <- diff_genes
  backgroupDB <- backgroupDB
  
  gene_df <- bitr(diff_genes, 
                  fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = backgroupDB)
  
  #### KEGG使用的是ENTREZID，需要重新映射一下啊
  kk <- enrichKEGG(gene = gene_df$ENTREZID,
                   organism = organism,
                   keyType = "kegg",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   pAdjustMethod = "BH"
  )
  
  #### 输出结果需要把 ENTREZID 重新映射为 Symbol
  kk@result$symbols <- apply(X = kk@result, MARGIN = 1, function(row){
    geneIDs <- row["geneID"] %>% as.character() %>% strsplit("/") %>% unlist() 
    symbols <- bitr(geneID = geneIDs, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = backgroupDB)$SYMBOL
    symbols <- paste0(symbols, collapse = ",")
    return(symbols)
  })
  
  p.kegg <- dotplot(kk, 
          showCategory = 10,  # 显示top10通路
          color = "p.adjust", # 颜色表示显著性
          size = "Count",     # 大小表示基因数量
          title = "KEGG Pathway Enrichment") +
    scale_color_gradient(low="red", high="blue") # 自定义颜色
  
  return(list(kk, p.kegg))
  
}









# ====================================================================
#  导入数据
# ====================================================================

### 读取元数据 ==================================
meta_table <- read.table(metafile, sep = ",", header = T, row.names = 1)
meta_table$GroupID <- as.factor(meta_table$GroupID)

treated <- setdiff(meta_table$GroupID %>% as.vector() %>% unique(), untreated)

### 从gtf文件中提取gene_id和transcript_id信息 =====
gtf <- import(gtf)

tx2gene <- gtf %>%
  as.data.frame() %>%
  dplyr::select(transcript_id, gene_id) %>% # 顺序必须是先转录本ID，再基因ID
  dplyr::filter(!is.na(transcript_id)) %>%
  unique()

genenamemapfile <- gtf %>% 
  as.data.frame() %>% 
  dplyr::select(gene_id, gene_name ) %>% 
  unique() 
rownames(genenamemapfile) <- genenamemapfile$gene_id


### 读取RNA表达数量信息 ==============================

### 导入转录本水平数据汇总到基因水平
counts.imported <- tximport::tximport(files = as.character(mappedfiles), type = 'salmon', tx2gene = tx2gene)

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
  # res <-  results(dds,  contrast = c("GroupID", t, untreated))   

  ### 对不可靠的倍数变化估计进行贝叶斯收缩 ==============
  ### 当样本量较少时这种收缩是很有必要的
  LFC <- DESeq2::lfcShrink(dds, coef = pairname, type = "apeglm") # coef 指定了比较组
  ### The contens of the LFC dataframe contain the log2 fold-change, as well as the p-value and adjusted p-value
  LFC.result <- as.data.frame(LFC)
  LFC.result$gene_name <- genenamemapfile[rownames(LFC), "gene_name"]

  
  ###绘制火山图 ========================================
  p <- EnhancedVolcano(LFC.result,
                  lab = LFC.result$gene_name,
                  x = 'log2FoldChange',
                  y = 'pvalue',
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
  ggsave(filename = filename, plot = p, width = 10, height = 7)
  
  
  ### 输出差异分析表格 ==================================
  LFC.result <- LFC.result[order(LFC.result$pvalue),]
  write.csv(as.data.frame(LFC.result), paste0(pairname, "_DEseq2_results.csv"))
  
  
  ### 绘制热图 ==========================================
  sig_genes <- subset(LFC, pvalue < pval & abs(log2FoldChange) > lfc) %>% rownames()
  groupA <- meta_table %>% filter(GroupID == t) %>% rownames()
  groupB <- meta_table %>% filter(GroupID == untreated) %>% rownames() 
  #pheatmap( assay(dds)[sig_genes,c(groupA, groupB)], scale = "row" )
  p.heatmap <- pheatmap(assay(vst(dds))[sig_genes,c(groupA, groupB)], scale = "row", show_rownames = F)
  filename <- paste0(pairname, ".Heatmap.pdf")
  ggsave(plot = p.heatmap, filename = "~/Desktop/04.湘湖实验室/李攀RNA-seq/GroupID_wko_vs_mKO.heatmap.pdf", width = 6, height = 6  )

  
  ### GO注释 和 KEGG注释 =================================
  
  #### 所有差异基因
  sig_genes.Symbol <- LFC.result %>% 
    filter(pvalue < pval & abs(log2FoldChange) > lfc) %>% 
    filter(!is.na(gene_name)) %>% 
    pull(gene_name)
  go.result <- GO_enrichment(sig_genes.Symbol, backgroupDB = backgroupDB)
  filename <- paste0(pairname, ".all.GO.csv")
  write.csv(go.result[[1]]@result, file = filename)
  filename <- paste0(pairname, ".all.GO.pdf")
  ggsave(filename = filename, plot = go.result[[2]], width = 10, height = 10)
  
  kegg.result <- KEGG_enrchment(sig_genes.Symbol, backgroupDB = backgroupDB)
  filename <- paste0(pairname, ".all.KEGG.csv")
  write.csv(kegg.result[[1]]@result, file = filename)
  filename <- paste0(pairname, ".all.KEGG.pdf")
  ggsave(filename = filename, plot = kegg.result[[2]], width = 10, height = 10)
  
  
  
  
  #### case 高表达的基因，log2FoldChange > 0的基因
  sig_genes.Symbol <- LFC.result %>% 
    filter(pvalue < pval & log2FoldChange > lfc) %>% 
    filter(!is.na(gene_name)) %>% 
    pull(gene_name)
  go.result <- GO_enrichment(sig_genes.Symbol, backgroupDB = backgroupDB)
  filename <- paste0(pairname,".", t ,".GO.csv")
  write.csv(go.result[[1]]@result, file = filename)
  filename <- paste0(pairname,".", t, ".GO.pdf")
  ggsave(filename = filename, plot = go.result[[2]], width = 10, height = 10)
  
  kegg.result <- KEGG_enrchment(sig_genes.Symbol, backgroupDB = backgroupDB)
  filename <- paste0(pairname,".", t ,".KEGG.csv")
  write.csv(kegg.result[[1]]@result, file = filename)
  filename <- paste0(pairname,".", t, ".KEGG.pdf")
  ggsave(filename = filename, plot = kegg.result[[2]], width = 10, height = 10)
  
  #### control 高表达的基因，log2FoldChange < 0的基因
  sig_genes.Symbol <- LFC.result %>% 
    filter(pvalue < pval & log2FoldChange < -lfc) %>% 
    filter(!is.na(gene_name)) %>% 
    pull(gene_name)
  go.result <- GO_enrichment(sig_genes.Symbol, backgroupDB = backgroupDB)
  filename <- paste0(pairname,".", untreated ,".GO.csv")
  write.csv(go.result[[1]]@result, file = filename)
  filename <- paste0(pairname,".", untreated, ".GO.pdf")
  ggsave(filename = filename, plot = go.result[[2]], width = 10, height = 10)
  
  kegg.result <- KEGG_enrchment(sig_genes.Symbol, backgroupDB = backgroupDB)
  filename <- paste0(pairname,".", untreated ,".KEGG.csv")
  write.csv(kegg.result[[1]]@result, file = filename)
  filename <- paste0(pairname,".", untreated, ".KEGG.pdf")
  ggsave(filename = filename, plot = kegg.result[[2]], width = 10, height = 10)

  
  
  
  
}







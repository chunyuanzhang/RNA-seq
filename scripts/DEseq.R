#setwd("/home/zhangchunyuan/zhangchunyuan/RNA-seq")
#.libPaths("/home/zhangchunyuan/tools/rstudio/lib/R/library")


suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))

option_list <- list(
  make_option("--infotable", type="character", default=NULL, help="差异倍数"),
  make_option("--lfc", type="double", default=1, help="差异倍数"),
  make_option("--pval", type="double", default=0.05, help="显著性阈值，转录组一般使用padj"),
  make_option("--design", type="character", default=NULL, help="指定比较组"),
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
design <- args$design

if(design == "Species"){
  orthologgenes <- args$orthologgenes
  cat("与种内差异不同的是，中间差异需要先提取1:1同源的基因列表，再进行标准话，随后进行差异分析\n")
}

untreated <- args$untreated
CountingFiles.isoforms <- args$CountingFiles.isoforms
CountingFiles.isoforms <- CountingFiles.isoforms %>% strsplit(",") %>% unlist()
CountingMethod <- args$CountingMethod
outdir <- args$outdir


# lfc = 1
# pval = 0.05
# untreated <- "Chicken"
# infotable <- "/home/zhangchunyuan/zhangchunyuan/RNA-seq/infotable.csv"
# CountingMethod <- "rsem"
# outdir <- "result/DEseq"
# CountingFiles.isoforms <- list.files(path = "result/RSEM", pattern = ".isoforms.results", full.names = T)
# orthologgenes <- "/home/zhangchunyuan/zhangchunyuan/Orthologgenes/Chicken_vs_Duck.one_to_one.orthologenes.txt"



#====================================================================
# 读取元信息
#====================================================================

meta_table <- read.table(infotable, sep = ",", header = T)
rownames(meta_table) <- meta_table$SampleID
treated <- setdiff(meta_table[[design]] %>% as.vector() %>% unique(), untreated)


#################
#################

if(design == "Species"){

  #====================================================================
  # 读取1:1的同源基因，并给名称不同的同源基因重命名，用于后续提取基因
  #====================================================================

  meta_table_list <- split(x = meta_table, f = ~Species)
  species <- names(meta_table_list)

  orthologgenes <- read.delim(file = orthologgenes)
  orthologgenes <- orthologgenes %>% mutate(
    genename =  ifelse(!!as.name(species[1]) == !!as.name(species[2]),!!as.name(species[1])  ,  paste0(!!as.name(species[1]), "_vs_", !!as.name(species[2])))
  ) 
  
  #====================================================================
  # 1、先在种内分别归一化、然后提取同源基因的子集、再归一化、差异分析
  #====================================================================

  dds_list <- list()
  size_factors_combined <- c()
  merged_raw_counts <- c()

  for(s in names(meta_table_list)){
    
    sub_table <- meta_table_list[[s]]
    pattern <- rownames(sub_table) %>% paste0(collapse = "|")
    sub_isoforms <- CountingFiles.isoforms[grep(pattern = pattern, x = CountingFiles.isoforms)] # 按照meta_table的顺序进行调整
    cat("正在读取物种 ",s," 对应样本 ", pattern, "并按照同源基因提取子集")
    tx2gene <- read.delim(file = sub_isoforms[1], header = T) %>% select(transcript_id, gene_id)
    counts.genes <- tximport::tximport(files = as.character(sub_isoforms), 
                                      type = CountingMethod, 
                                      tx2gene = tx2gene,
                                      countsFromAbundance = "lengthScaledTPM" )  # 将种内各样本的表达数据导入
    
    ### 统计size factor，并保存在 size_factors_combined 中，后续需要继续使用
    dds_list[[s]] <- DESeq2::DESeqDataSetFromTximport(counts.genes, colData = sub_table, design = ~1)
    dds_list[[s]] <- DESeq2::estimateSizeFactors(dds_list[[s]])
    size_factors_combined <- c(size_factors_combined, sizeFactors(dds_list[[s]]) )
    ### 提取同源基因的原始counts
    raw_counts_array <- counts(dds_list[[s]], normalized = FALSE)[orthologgenes[[s]],]
    rownames(raw_counts_array) <- orthologgenes[["genename"]]
    
    ### 合并矩阵
    merged_raw_counts <- cbind(merged_raw_counts, raw_counts_array)

  }
  ### 重新给 meta_table 排序，保证矩阵中样本的顺序和表格中保持一致
  meta_table <- meta_table[colnames(merged_raw_counts), ]

  ### 将合并后的数据重新创建DEseq对象
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = merged_raw_counts,
    colData = meta_table,
    design = ~ Species  # 现在可以指定species作为变量
  )

  dds$Species <- relevel(dds$Species, ref = untreated)

  ### 将种内各自估计的size factor赋给dds
  sizeFactors(dds) <- size_factors_combined
  ### 该步骤会使用已经存在的size factor进行分析
  dds <- DESeq2::DESeq(dds)

} else {

  tx2gene <- read.delim(file = CountingFiles.isoforms[1], header = T) %>% dplyr::select(transcript_id, gene_id)
  counts.genes <- tximport::tximport(files = as.character(CountingFiles.isoforms), type = CountingMethod, tx2gene = tx2gene )

  ###创建DESeq2对象 
  design_formula <- as.formula(paste("~", design))
  dds <- DESeq2::DESeqDataSetFromTximport(counts.genes, colData = meta_table, design = design_formula )
  # 不指定的话会根据字符顺序比较
  dds[[design]] <- relevel(dds[[design]], ref = untreated)

  # 运行DESeq2
  dds <- DESeq2::DESeq(dds)  # 这一步自动进行了标准化

}


#################
#################

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
  pairname <- paste0(design, "_", t, "_vs_", untreated)
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
  
  
  ### 输出差异分析总表格 ==================================
  result <- result[order(result$padj),]
  result_dataframe <- as.data.frame(result) %>% 
    tibble::rownames_to_column(var = "genename") %>%
    mutate(
      Up = case_when(
        padj < pval & log2FoldChange > lfc ~ t,
        padj < pval & log2FoldChange < -lfc ~ untreated,
        TRUE ~ "NotSig"
      )
    ) 
  # ========
  if(design == "Species"){
    result_dataframe <- result_dataframe %>%
      left_join(., orthologgenes, by = "genename")
  }
  # ========

  write.csv(result_dataframe, file.path(outdir,paste0(pairname, ".all.csv")),
    quote = F,  row.names = F)
  
  
  ### 绘制热图 ==========================================
  sig_genes <- subset(result, padj < pval & abs(log2FoldChange) > lfc) %>% rownames()
  groupA <- meta_table %>% filter(!!as.name(design) == t) %>% rownames()
  groupB <- meta_table %>% filter(!!as.name(design) == untreated) %>% rownames() 
  #pheatmap( assay(dds)[sig_genes,c(groupA, groupB)], scale = "row" )
  p.heatmap <- pheatmap::pheatmap(assay(vst(dds))[sig_genes,c(groupA, groupB)], scale = "row", show_rownames = F)
  filename <- paste0(pairname, ".Heatmap.pdf")
  ggsave(plot = p.heatmap, filename = file.path(outdir,filename), width = 6, height = 6, create.dir = TRUE )
  
}



# 这个脚本只在种间比较的时候使用

source("scripts/functions.R")

suppressMessages({
  library(dplyr)
  library(optparse)
  library(tximport)
  library(DESeq2)
  library(ggplot2)
})


#-------------------------------------------------------------------------------
# 参数传递
#-------------------------------------------------------------------------------

option_list <- list(
  make_option("--infotable", type="character", default=NULL, help="差异倍数"),
  make_option("--lfc", type="double", default=1, help="差异倍数"),
  make_option("--padj", type="double", default=0.05, help="显著性阈值，转录组一般使用padj"),
  make_option("--design", type="double", default=0.05, help="在谁之间比对"),
  make_option("--untreated", type="character", default=NULL, help="指定对照组"),
  make_option("--confoundingvariable", type="character", default=NULL, help="混淆变量，多个变量用逗号分隔"),
  make_option("--Orthologgenes", type="character", default=NULL, help="1:1同源基因列表"),
  make_option("--indir", type = "character", default = "result/Salmon/", help = "Samlon输出的主文件夹"),
  make_option("--outdir", type = "character", default = "result/Salmon/", help = "差异分析结果存储路径")
)

args <- parse_args(OptionParser(option_list=option_list))

infotable <- args$infotable
lfc <- args$lfc
padj <- args$padj
design <- args$design
untreated <- args$untreated
confoundingvariable <- args$confoundingvariable
confoundingvariable <- confoundingvariable %>% strsplit(split = ",") %>% unlist() %>% trimws()
indir <- args$indir
indir <- check_path(indir)
outdir <- args$outdir
outdir <- check_path(outdir)


#-------------------------------------------------------------------------------
# 读取元信息
#-------------------------------------------------------------------------------

# indir <- "Salmon"
# untreated <- "Chicken"
# infotable <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/MyAnalyses/infotable.csv"
# Orthologgenes <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/同源基因/Chicken_vs_Duck.Ensembl_plus_Orthofinder.one_to_one.orthologenes.tsv"
# lfc <- 1
# padj <- 0.05

sampletable <- read.csv(infotable, check.names = F, comment.char = "#", blank.lines.skip = T) %>% filter(!is.na(Transcriptome)) 
rownames(sampletable) <- sampletable$SampleID

treated <- setdiff(sampletable[[design]] %>% as.vector() %>% unique(), untreated)
# 如果存在多组比较，则直接根据treated循环比较
sampletable[design] <- factor(x = sampletable[[design]], levels = c(untreated, treated))

#print(head(sampletable))


#-------------------------------------------------------------------------------
# 导入salmon结果
#-------------------------------------------------------------------------------


if(design == "Species"){
  
  print("种间差异分析，读取表达数据")
  # 读取同源基因
  Orthologgenes <- args$Orthologgenes
  Orthologgenes <- read.delim(file = Orthologgenes, header = T) %>%
    mutate(
      CombindName = paste0(!!as.name(untreated), "_vs_", !!as.name(treated))
    )
  
  # 拆分sample table
  sampletable_list <- split(x = sampletable, f = sampletable[[design]])
  species <- names(sampletable_list)
  
  # 1. 各物种导入(保留 length,不要 lengthScaledTPM)
  txi_list <- list()
  for (s in species) {
    tx2gene <- read.table(file = paste0(indir, s, "_trans2symbol.tsv"), header = T)
    
    files <- file.path(indir, rownames(sampletable_list[[s]]), "quant.sf")
    names(files) <- rownames(sampletable_list[[s]])
    txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
    txi_list[[s]] <- txi.salmon
    remove(tx2gene, txi.salmon, files)
  }
  
  # 2. 按同源基因抽 counts + length,合并成一个 txi
  keep <- Orthologgenes[[untreated]] %in% rownames(txi_list[[untreated]]$counts) &
    Orthologgenes[[treated]] %in% rownames(txi_list[[treated]]$counts)
  ort <- Orthologgenes[keep, ]
  
  pull <- function(txi, g)
    lapply(txi[c("counts","length","abundance")], function(m){
      m <- m[g, , drop = FALSE]; rownames(m) <- ort$CombindName; m })
  a <- pull(txi_list[[untreated]], ort[[untreated]])
  b <- pull(txi_list[[treated]], ort[[treated]])
  
  txi <- list(
    counts    = cbind(a$counts,    b$counts),
    length    = cbind(a$length,    b$length),
    abundance = cbind(a$abundance, b$abundance),
    countsFromAbundance = "no"
  )
  
  sampletable <- sampletable[colnames(txi$counts), ] 
  
  print("完成种间数据合并")
  
  
} else {
  
  print("种内差异分析，读取表达数据")
  s <- sampletable$Species %>% unique() %>% as.character()
  tx2gene <- read.table(file = paste0(indir, s, "_trans2symbol.tsv"), header = T)
  files <- file.path(indir, rownames(sampletable), "quant.sf")
  names(files) <- rownames(sampletable)
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
  remove(tx2gene,  files)
  sampletable <- sampletable[colnames(txi$counts), ] 

  print("完成Samlon数据读取")
  
}


#-------------------------------------------------------------------------------
# 创建DEseq对象，并进行差异表达分析
#-------------------------------------------------------------------------------


# 拼接所有项：协变量 + 主效应
allvariable <- c(confoundingvariable, design)


dds <- DESeqDataSetFromTximport(
  txi,
  sampletable,
  design = reformulate(allvariable)
)

dds[[design]] <- relevel(dds[[design]], ref = untreated)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)


for (t in treated) {
  # 4. 差异基因
  res <- results(dds, contrast = c(design, t, untreated)) %>% as.data.frame()
  # res <- results(dds, lfcThreshold = lfc, altHypothesis = "greaterAbs", contrast = c(design, t, untreated)) %>% as.data.frame()   # 检验 H₀「|LFC| ≤ 1」, 把变化倍数显著的提取出来
  res[[untreated]] <- extract_gene_from_gene_vs_gene_pair(genes = rownames(res), side = 1)
  res[[t]] <- extract_gene_from_gene_vs_gene_pair(genes = rownames(res), side = 2)
  res <- res[order(-abs(res$stat)), ]
  
  outfile <- paste0(outdir, paste0(design, "_", t, "_vs_", untreated, ".all.tsv"))
  write.table(x = res, file = outfile, sep = "\t", row.names = F, quote = F)
  
  deg <- res %>% filter(abs(log2FoldChange) > lfc, padj < !!padj, baseMean > 10, lfcSE < 1)
  deg <- deg[order(-abs(deg$stat)), ]
  
  outfile <- paste0(outdir, paste0(design, "_", t, "_vs_", untreated, ".diffexp.tsv"))
  write.table(x = deg, file = outfile, sep = "\t", row.names = F, quote = F)
  

  remove(res, deg, outfile)
  
}

# log2FoldChange < 0 是在鸡[untreated]中差异高表达
# log2FoldChange > 0 是在鸭[treated]中差异高表达

print("完成差异分析")
print("命名方式为 untreated_vs_treated, 例如 Chicken_vs_Duck, log2FoldChange < 0 是在鸡中差异高表达, log2FoldChange > 0 是在鸭中差异高表达")



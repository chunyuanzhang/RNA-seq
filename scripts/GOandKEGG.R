.libPaths("/home/zhangchunyuan/tools/rstudio/lib/R/library")

suppressMessages(library(dplyr))
suppressMessages(library(VennDiagram))
suppressMessages(library(optparse))

# suppressMessages(library(org.Gg.eg.db))
# library(clusterProfiler)
# library(enrichplot)
# library(org.AnasPlatyrhynchos.domestica.eg.sqlite)


# ====================================================================
#  注释背景文件选择
# ====================================================================


if(stringr::str_to_lower(species) == "chicken"){
  backgroupDB = org.Gg.eg.db
  organism = "gga"
}

if(stringr::str_to_lower(species) == "duck"){
  backgroupDB = org.AnasPlatyrhynchos.domestica.eg.sqlite 
  organism 
}


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
  make_option("--CountingMethod", type="character", default="rsem", help="表达量统计工具，影响数据导入的格式"),
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
CountingMethod <- args$CountingMethod
mappedfiles <- mappedfiles %>% strsplit(",") %>% unlist()

species <- args$species


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
  ggsave(filename = file.path(outdir, filename), plot = kegg.result[[2]], width = 10, height = 10)
  
  
  
  
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

  
  












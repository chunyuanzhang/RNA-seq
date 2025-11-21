
suppressMessages(library(clusterProfiler))
# suppressMessages(library(enrichplot))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(GO.db))

cat("使用emapper获取注释信息，然后使用enrichKEGG进行富集分析；\n在进行物种间的富集分析时，只对物种内差异高表达的基因进行注释;\n当前物种间注释时只使用了表达的同源基因，非同源基因目前并不包含在本分析中，若要加进来，只需要更改backgroundgenes即可\n")


source("scripts/functions.R")

option_list <- list(
  make_option("--target", type="character", default = NULL, help = "希望进行KEGG分析的物种或分组名称"),
  make_option("--diffgenes",type="character",  default=NULL, help="差异分析提取的差异基因，用于进行注释" ),
  make_option("--backgroundgenes",type="character",  default=NULL, help="所有参与表达的基因作为背景基因，这里暂时不包含非同源的基因" ),
  make_option("--pval", type="double",  default=0.05, help = "富集分析的显著性阈值"),
  make_option("--emapperannotations", type="character", default=NULL, help="emapper比对结果文件，在没有现成的注释数据库时使用"),
  make_option("--gtf",type="character", default = NULL, help = "gtf文件"),
  make_option("--outdir", type = "character", default = "result/DEseq/", help = "差异分析结果存储路径")
)

args <- parse_args(OptionParser(option_list=option_list))


target <- args$target
diffgenes <- args$diffgenes
backgroundgenes <- args$backgroundgenes
pval <- args$pval
emapperannotations <- args$emapperannotations
gtf <- args$gtf
outdir <- args$outdir

# target <- "Chicken"
# diffgenes <- "result/DEseq/Species_Duck_vs_Chicken.all.csv"
# backgroundgenes <- "result/DEseq/Species_Duck_vs_Chicken.all.csv"
# pval <- 0.05
# emapperannotations <- "/home/zhangchunyuan/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b/chicken.emapper.annotations"
# gtf <- "/home/zhangchunyuan/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf"
# outdir <- "result/DEseq"

pairname <- basename(diffgenes) %>% strsplit("\\.") %>% unlist() %>% head(n=1)

# ====================================================================
# 读取GTF文件，方便把蛋白质ID和基因ID对应
# 读取 eggNOG-mapper 输出的结果
# 非模式物种的注释信息偏少，需要通过比对获取注释信息；
# 本流程采用 eggNOG-mapper 线上流程进行分析，直接把参考基因组的蛋白质序列文件上传进行比对即可，然后保存输出文件为emapperannotations
# ====================================================================

diffgenes <- read.csv(diffgenes)
backgroundgenes <- read.csv(backgroundgenes)

# ====================================================================
# 从差异分析的结果中读取目标基因和背景基因
# ====================================================================

### 读取目标差异的基因
if(target %in% colnames(diffgenes)){
  side <- ifelse(target == "Chicken", yes = 1, no = 2)
  diffgenes <-  diffgenes  %>% filter(Up == target) %>% pull(target) # 这里已经考虑物种内自己的基因名了
  backgroundgenes <- backgroundgenes %>% pull(target)  # 所有表达的基因作为“背景基因”，并不是基因组上的所有基因，注意有些基因并不表达
}else{
  ### 如果是物种内的比较，直接提取genename即可
  diffgenes <-  diffgenes %>% filter(Up == target) %>% pull(genename)
  backgroundgenes <- backgroundgenes %>% pull(genename)
}

diffgenes <- extract_gene_from_gene_vs_gene_pair(genes = diffgenes, side = side)
backgroundgenes <- extract_gene_from_gene_vs_gene_pair(genes = backgroundgenes, side = side )

# ==========================================================
# GO富集分析
# ==========================================================

enrichment_GO <- GO_enrichment(diffgenes =  diffgenes, backgroundgenes = backgroundgenes, GTF = gtf, emapperannotations = emapperannotations)

### 输出富集表格
write.csv(x = enrichment_GO@result, file = file.path(outdir, paste0(pairname, ".", target, ".GO.csv")), quote = F,  row.names = F)


# ==============================================================
# KEGG富集分析
# ==============================================================

enrichment_kegg <- KEGG_enrichment(diffgenes =  diffgenes, backgroundgenes = backgroundgenes, GTF = gtf, emapperannotations = emapperannotations)

### 输出结果
write.csv(x = enrichment_kegg@result, file = file.path(outdir, paste0(pairname, ".", target, ".KEGG.csv")), quote = F, row.names = F)




suppressMessages(library(clusterProfiler))
# suppressMessages(library(enrichplot))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))


option_list <- list(
  make_option("--infotable", type="character", default=NULL, help="差异倍数"),
  make_option("--lfc", type="double", default=1, help="差异倍数"),
  make_option("--pval", type="double", default=0.05, help="显著性阈值，转录组一般使用padj"),
  make_option("--untreated", type="character", default=NULL, help="指定对照组"),
  make_option("--emapperannotations", type="character", default=NULL, help="emapper比对结果文件，在没有现成的注释数据库时使用"),
  make_option("--diffgenes",type="character",  default="result/DEseq/diffgenes.tsv", help="差异分析提取的差异基因，用于进行注释" ),
  make_option("--gtf",type="character", default = NULL, help = "gtf文件"),
  make_option("--target", type="character", default = NULL, help = "希望进行KEGG分析的物种或分组名称"),
  make_option("--outdir", type = "character", default = "result/DEseq/", help = "差异分析结果存储路径")
)

args <- parse_args(OptionParser(option_list=option_list))


infotable <- args$infotable
lfc <- args$lfc
pval <- args$pval
untreated <- args$untreated
target <- args$target
emapperannotations <- args$emapperannotations
diffgenes <- args$diffgenes
gtf <- args$gtf
outdir <- args$outdir



# infotable <- "infotable.csv"
# lfc <- 1
# pval <- 0.05
# untreated <- "Chicken"
# target <- "Chicken"
# emapperannotations <- "/home/zhangchunyuan/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b/chicken.emapper.annotations"
# gtf <- "/home/zhangchunyuan/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf"
# diffgenes <- "result/DEseq/diffgenes.csv"
# outdir <- "result/DEseq"

# ====================================================================
#   读取元信息
# ====================================================================

meta_table <- read.table(infotable, sep = ",", header = T, row.names = 1)
meta_table$Species <- as.factor(meta_table$Species)
treated <- setdiff(meta_table$Species %>% as.vector() %>% unique(), untreated)


# ====================================================================
# 读取GTF文件，方便把蛋白质ID和基因ID对应
# ====================================================================


GTF <- rtracklayer::import(gtf)
GTF <- as.data.frame(GTF)
gene_to_protein_mapping <- GTF %>% select(gene_id, protein_id) %>% filter(!is.na(protein_id)) %>% unique() 


# ====================================================================
# 读取 eggNOG-mapper 输出的结果
# ====================================================================


header <- readLines(emapperannotations, n = 5)
header <- header[5] %>% gsub("#", "", .) %>% strsplit(split = "\t") %>% unlist()
data <- read.delim(emapperannotations, comment.char = '#', header = F)
colnames(data) <- header
remove(header)

### 提取KEGG相关注释信息
eggNOG_kegg <- data %>% select(KEGG_ko, query) %>% 
  tidyr::separate_rows(KEGG_ko, sep = ",") %>%
  filter(KEGG_ko != "-") %>%
  mutate(KEGG_ko = gsub(pattern = "ko:", replacement = "", KEGG_ko)) %>%
  unique() %>% 
  tidyr::drop_na() %>% 
  dplyr::rename("protein_id" = "query") %>%
  merge(., gene_to_protein_mapping, by = "protein_id") %>%
  select(KEGG_ko, gene_id) %>%
  unique()


### 提取GO相关注释信息




# ====================================================================
# 
# ====================================================================


keggpathways <- c()

for (t in treated) {

  ### 读取背景基因
  pairname <- paste0("Species","_", t, "_vs_", untreated)
  alldifffile <- file.path(dirname(diffgenes),  paste0(pairname, "_DEseq2_results.csv"))
  background_genes <- read.csv(alldifffile) %>% pull(target)  # 背景基因

  ### 将背景基因的KEGG号提取出来
  background_kegg <- eggNOG_kegg %>%
    dplyr::filter(gene_id %in% background_genes) %>%
    unlist() %>%
    as.vector()
  
  ### 读取目标差异的基因
  interesting_set <-  read.csv(diffgenes)  %>% filter(Up == target) %>% pull(target)

  ### 把差异基因对应的KEGG号提取出来
  interesting_set_kegg <- eggNOG_kegg %>%
    dplyr::filter(gene_id %in% interesting_set) %>%
    unlist() %>%
    as.vector()
  
  
  enrichment_kegg <- enrichKEGG(interesting_set_kegg,
                                organism = "ko",
                                keyType = "kegg",
                                pvalueCutoff = pval,
                                qvalueCutoff = pval,
                                pAdjustMethod = "BH",
                                universe = background_kegg,
                                use_internal_data = FALSE)
  
  result <- enrichment_kegg@result
  
  result$gene_ids <- apply(result, 1, function(row){
    genekeggs <-  row["geneID"] %>% strsplit("/") %>% unlist()
    gene_ids <- eggNOG_kegg[match(genekeggs, eggNOG_kegg$KEGG_ko),]$gene_id %>% paste0(collapse = "/")
    return(gene_ids)
  })
  
  write.csv(x = result, file = file.path(outdir, paste0(pairname, ".", target, ".KEGG.csv")), quote = F, col.names = T, row.names = F)
  
  keggpathways <- rbind(keggpathways, result %>% filter(pvalue < 0.05) %>% mutate(compare = pairname) )
  
}


write.csv(x = keggpathways, file = file.path(outdir, paste0(target, ".KEGG.csv")), quote = F, row.names = F)












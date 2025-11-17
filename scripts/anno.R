
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

GTF <- rtracklayer::import(gtf)
GTF <- as.data.frame(GTF)
emapperannotations_data <- get_eggNOGmapperData(emapperannotations = emapperannotations)
diffgenes <- read.csv(diffgenes)
backgroundgenes <- read.csv(backgroundgenes)

# ====================================================================
# 将基因和蛋白对应起来
# ====================================================================

gene_to_protein_mapping <- get_Protein2Gene(GTF)

# ====================================================================
# 读取目标基因和背景基因调整为KO号，用于KEGG注释
# ====================================================================

### 提取KEGG相关注释信息
eggNOG_kegg2gene <- get_KO2Gene(emapperannotations_data = emapperannotations_data, gene_to_protein = gene_to_protein_mapping)


### 读取目标差异的基因
if(target %in% colnames(diffgenes)){

  ### 如果是物种间比较，diffgenes的格式稍有不同
  diff_genes <-  diffgenes  %>% filter(Up == target) %>% pull(target) # 这里已经考虑物种内自己的基因名了
  ### 把差异基因对应的KEGG号提取出来
  diff_genes_kegg <- eggNOG_kegg2gene %>% dplyr::filter(gene_id %in% diff_genes) %>% pull(KEGG_ko) 
  background_genes <- backgroundgenes %>% pull(target)  # 所有表达的基因作为“背景基因”，并不是基因组上的所有基因，注意有些基因并不表达
  ### 将背景基因的KEGG号提取出来
  background_genes_kegg <- eggNOG_kegg2gene %>%
      dplyr::filter(gene_id %in% background_genes) %>%
      pull(KEGG_ko) 

}else{

  ### 如果是物种内的比较，直接提取genename即可
  diff_genes <-  diffgenes %>% filter(Up == target) %>% pull(genename)
  diff_genes_kegg <- eggNOG_kegg2gene %>% dplyr::filter(gene_id %in% diff_genes) %>% pull(KEGG_ko) 
  background_genes <- backgroundgenes %>% pull(genename)
  background_genes_kegg <- eggNOG_kegg2gene %>%
    dplyr::filter(gene_id %in% background_genes) %>%
    pull(KEGG_ko) 

}


# ==============================================================
# KEGG富集分析
# ==============================================================


    
### 富集分析
enrichment_kegg <- enrichKEGG(diff_genes_kegg,
                              organism = "ko",
                              keyType = "kegg",
                              pvalueCutoff = pval,
                              qvalueCutoff = pval,
                              pAdjustMethod = "BH",
                              universe = background_genes_kegg,
                              use_internal_data = FALSE)


### 给富集分析的结果添加必要的注释信息
enrichment_kegg@result <- enrichment_kegg@result %>%
  rowwise() %>%
  mutate(
    # 添加 gene_ids
    gene_ids = {
      genekeggs <- strsplit(geneID, "/")[[1]]
      gene_ids <- eggNOG_kegg2gene[eggNOG_kegg2gene$KEGG_ko %in% genekeggs, ]$gene_id %>% 
        unique()
      gene_ids[gene_ids %in% diff_genes] %>% 
        paste0(collapse = "/")
    },
    # 添加 CLASS
    CLASS = {
      tryCatch({
        keggfile <- KEGGREST::keggGet(ID)
        CLASS <- keggfile[[1]]$CLASS
        if(is.null(CLASS)) NA_character_ else CLASS
      }, error = function(e) NA_character_)
    }
  ) %>%
  ungroup() %>%
  # 拆分 CLASS 为 L1 和 L2
  tidyr::separate(CLASS, into = c("L1", "L2"), sep = "; ", remove = FALSE, fill = "right")


### 输出结果
write.csv(x = enrichment_kegg@result, file = file.path(outdir, paste0(pairname, ".", target, ".KEGG.csv")), quote = F, row.names = F)
  

# ==========================================================
# GO富集分析
# ==========================================================


### 提取GO相关注释信息
eggNOG_goterm2gene <- get_GOterm2Gene(emapperannotations_data = emapperannotations_data, gene_to_protein = gene_to_protein_mapping)

goterm2name <- AnnotationDbi::select(GO.db::GO.db, keys=unique(eggNOG_goterm2gene$term), 
                                   columns=c("GOID", "TERM"), 
                                   keytype="GOID") 
colnames(goterm2name) <- c("term", "name")



### GO富集, 使用基因symbol
enrichment_GO <- enricher(gene = diff_genes, 
           TERM2GENE = eggNOG_goterm2gene,
           TERM2NAME = goterm2name,
           pvalueCutoff = 0.05, 
           qvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           universe = background_genes,
           )

enrichment_GO@result <- enrichment_GO@result %>%
  left_join(
    select(GO.db::GO.db,
           keys = .$ID,
           columns = c("GOID", "ONTOLOGY"),
           keytype = "GOID"),
    by = c("ID" = "GOID")
  )


### 输出富集表格
write.csv(x = enrichment_GO@result, file = file.path(outdir, paste0(pairname, ".", target, ".GO.csv")), quote = F,  row.names = F)
  
# ### 可视化
# # 按ONTOLOGY分面的条形图
# top_terms <- enrichment_GO@result %>%
#     group_by(ONTOLOGY) %>%
#     slice_min(order_by = pvalue, n = 10) %>%
#     ungroup()
# p.GO <- ggplot(top_terms, aes(x = reorder(Description, -log10(pvalue)), 
#                         y = -log10(pvalue), 
#                         fill = ONTOLOGY)) +
#     geom_bar(stat = "identity") +
#     coord_flip() +
#     facet_wrap(~ONTOLOGY, scales = "free_y", ncol = 1) +
#     labs(x = "GO Term", y = "-log10(P-value)", 
#          title = "GO Enrichment Analysis") +
#     theme_minimal() +
#     theme(legend.position = "none")
#   
# ggsave(filename = file.path(outdir, paste0(pairname, ".", target, ".GO.pdf")), plot = p.GO, width = 5, height = 6)
# 


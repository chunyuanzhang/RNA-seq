
suppressMessages(library(clusterProfiler))
# suppressMessages(library(enrichplot))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(GO.db))

cat("使用emapper获取注释信息，然后使用enrichKEGG进行富集分析；\n在进行物种间的富集分析时，只对物种内差异高表达的基因进行注释;\n当前物种间注释时只使用了表达的同源基因，非同源基因目前并不包含在本分析中，若要加进来，只需要更改backgroundgenes即可")


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
# ====================================================================

GTF <- rtracklayer::import(gtf)
GTF <- as.data.frame(GTF)
gene_to_protein_mapping <- GTF %>% dplyr::select(gene_id, protein_id) %>% filter(!is.na(protein_id)) %>% unique() 

# ====================================================================
# 读取 eggNOG-mapper 输出的结果
# 每次注释都使用该物种对应的emapperannotations文件，从而实现物种特异的注释
# ====================================================================

# 非模式物种的注释信息偏少，需要通过比对获取注释信息；
# 本流程采用 eggNOG-mapper 线上流程进行分析，直接把参考基因组的蛋白质序列文件上传进行比对即可，然后保存输出文件为emapperannotations

header <- readLines(emapperannotations, n = 5)
header <- header[5] %>% gsub("#", "", .) %>% strsplit(split = "\t") %>% unlist()
data <- read.delim(emapperannotations, comment.char = '#', header = F)
colnames(data) <- header
remove(header)

### 提取KEGG相关注释信息
eggNOG_kegg <- data %>% dplyr::select(KEGG_ko, query) %>% 
  tidyr::separate_rows(KEGG_ko, sep = ",") %>%
  filter(KEGG_ko != "-") %>%
  mutate(KEGG_ko = gsub(pattern = "ko:", replacement = "", KEGG_ko)) %>%
  unique() %>% 
  tidyr::drop_na() %>% 
  dplyr::rename("protein_id" = "query") %>%
  merge(., gene_to_protein_mapping, by = "protein_id") %>%
  dplyr::select(KEGG_ko, gene_id) %>%
  unique()


### 提取GO相关注释信息
eggNOG_term2gene <- data %>% dplyr::select(GOs, query) %>%
  tidyr::separate_rows(GOs, sep = ",") %>%
  filter(GOs != "-") %>%
  tidyr::drop_na() %>% 
  dplyr::rename("protein_id" = "query") %>%
  merge(., gene_to_protein_mapping, by = "protein_id") %>%
  dplyr::select(GOs, gene_id) %>%
  dplyr::rename("term" = "GOs", "gene" = "gene_id") %>%
  unique() 

term2name <- select(GO.db::GO.db, keys=unique(eggNOG_term2gene$term), 
                    columns=c("GOID", "TERM"), 
                    keytype="GOID") 
colnames(term2name) <- c("term", "name")


# ====================================================================
# 读取感兴趣的差异基因，以及背景基因【当前使用的是所有表达的同源基因，非同源基因没有参与分析】
# ====================================================================

### 读取背景基因
background_genes <- read.csv(backgroundgenes) %>% pull(target)  # 所有表达的基因作为“背景基因”，并不是基因组上的所有基因，注意有些基因并不表达

### 将背景基因的KEGG号提取出来
background_kegg <- eggNOG_kegg %>%
    dplyr::filter(gene_id %in% background_genes) %>%
    unlist() %>%
    as.vector()
  
### 读取目标差异的基因
interesting_set <-  read.csv(diffgenes)  %>% filter(Up == target) %>% pull(target) # 这里已经考虑物种内自己的基因名了

# ==============================================================
# KEGG富集分析
# ==============================================================
  
### 把差异基因对应的KEGG号提取出来
interesting_set_kegg <- eggNOG_kegg %>%
    dplyr::filter(gene_id %in% interesting_set) %>%
    unlist() %>%
    as.vector()
### 富集分析
enrichment_kegg <- enrichKEGG(interesting_set_kegg,
                              organism = "ko",
                              keyType = "kegg",
                              pvalueCutoff = pval,
                              qvalueCutoff = pval,
                              pAdjustMethod = "BH",
                              universe = background_kegg,
                              use_internal_data = FALSE)
  
### 把基因ID添加进来
enrichment_kegg@result$gene_ids <- apply(enrichment_kegg@result, 1, function(row){
    genekeggs <-  row["geneID"] %>% strsplit("/") %>% unlist()
    gene_ids <- eggNOG_kegg[match(genekeggs, eggNOG_kegg$KEGG_ko),]$gene_id %>% paste0(collapse = "/")
    return(gene_ids)
})

### 输出结果
write.csv(x = enrichment_kegg@result, file = file.path(outdir, paste0(pairname, ".", target, ".KEGG.csv")), quote = F, col.names = T, row.names = F)
  

  
# ==========================================================
# GO富集分析
# ==========================================================
  
### 把差异基因对应的GO号提取出来
interesting_set_GO <- eggNOG_term2gene %>%
    dplyr::filter(gene %in% interesting_set) %>%
    unlist() %>%
    as.vector()
  
### GO富集
enrichment_GO <- enricher(gene = interesting_set, 
           TERM2GENE = eggNOG_term2gene,
           TERM2NAME = term2name,
           pvalueCutoff = 0.05, 
           qvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           universe = background_genes,
           )
  
### 添加 ONTOLOGY （"BP" "CC" NA   "MF"） 信息
go_ontology <- select(GO.db::GO.db,
                    keys = enrichment_GO@result$ID,
                    columns = c("GOID", "ONTOLOGY"),
                    keytype = "GOID") 
enrichment_GO@result <- left_join(x = enrichment_GO@result, y = go_ontology, by = c("ID" = "GOID")) 
### 输出富集表格
write.csv(x = enrichment_GO@result, file = file.path(outdir, paste0(pairname, ".", target, ".GO.csv")), quote = F,  row.names = F)
  
### 可视化
# 按ONTOLOGY分面的条形图
top_terms <- enrichment_GO@result %>%
    group_by(ONTOLOGY) %>%
    slice_min(order_by = pvalue, n = 10) %>%
    ungroup()
p.GO <- ggplot(top_terms, aes(x = reorder(Description, -log10(pvalue)), 
                        y = -log10(pvalue), 
                        fill = ONTOLOGY)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~ONTOLOGY, scales = "free_y", ncol = 1) +
    labs(x = "GO Term", y = "-log10(P-value)", 
         title = "GO Enrichment Analysis") +
    theme_minimal() +
    theme(legend.position = "none")
  
ggsave(filename = file.path(outdir, paste0(pairname, ".", target, ".GO.pdf")), plot = p.GO, width = 5, height = 6)



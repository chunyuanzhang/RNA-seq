
suppressMessages(library(clusterProfiler))
# suppressMessages(library(enrichplot))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(GO.db))

option_list <- list(
    make_option("--diffgenes",type="character",  default=NULL, help="目标基因，文件只有一列，不要文件名，一行一个基因" ),
    make_option("--backgroundgenes",type="character",  default=NULL, help="所有参与表达的基因作为背景基因，这里暂时不包含非同源的基因，不要文件名，一行一个基因" ),
    make_option("--gtf", type="character", default=NULL, help="gtf文件"),
    make_option("--emapperannotations", type="character", default=NULL, help="emapper比对结果文件，在没有现成的注释数据库时使用，【注意】此时决定了物种"),
    make_option("--outdir", type = "character", default = "./", help = "输出文件路径")
)

args <- parse_args(OptionParser(option_list=option_list))

diffgenes <- args$diffgenes
backgroundgenes <- args$backgroundgenes
emapperannotations <- args$emapperannotations
gtf <- args$gtf
outdir <- args$outdir


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

term2name <- AnnotationDbi::select(GO.db::GO.db, keys=unique(eggNOG_term2gene$term), 
                    columns=c("GOID", "TERM"), 
                    keytype="GOID") 
colnames(term2name) <- c("term", "name")



# ====================================================================
# 读取感兴趣的差异基因，以及背景基因【当前使用的是所有表达的同源基因，非同源基因没有参与分析】
# ====================================================================

diffgenes <- read.csv(diffgenes, header = F)
backgroundgenes <- read.csv(backgroundgenes)






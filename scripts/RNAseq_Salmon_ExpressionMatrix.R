


suppressMessages({
  library(tximport)
  library(optparse)
  library(dplyr)
})
#-------------------------------------------------------------------------------
# 传递参数
#-------------------------------------------------------------------------------


check_path <- function(path) {
  sub("/?$", "/", path)
}

option_list <- list(
  make_option("--infotable", type="character", default=NULL, help="差异倍数"),
  make_option("--Species", type="character", default=NULL, help="物种，如果在种内分析，则不添加该参数"),
  make_option("--gff", type="character", default=NULL, help="参考基因组GFF文件，注意不是GTF"),
  make_option("--indir", type = "character", default = "result/Salmon/", help = "Salmon结果存储路径"),
  make_option("--outdir", type = "character", default = "result/Salmon/", help = "表达结果存储路径")
)


args <- parse_args(OptionParser(option_list=option_list))

infotable <- args$infotable
Species <- args$Species
gff <- args$gff
indir <- args$indir
outdir <- args$outdir


#infotable <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/MyAnalyses/sampleinfo.csv"
#Species <- "Duck"
#gff <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/MyAnalyses/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.gff"
#indir <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/MyAnalyses/RNA-seq/Salmon/"
#outdir <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/MyAnalyses/RNA-seq/Salmon/"


indir <- check_path(indir)
outdir <- check_path(outdir)

#-------------------------------------------------------------------------------
# 读取infotable，并按照物种信息提取样本
#-------------------------------------------------------------------------------

sampletable <- read.csv(infotable, check.names = F, comment.char = "#", blank.lines.skip = T) %>% filter(!is.na(Transcriptome)) %>% filter(Species == !!Species)
rownames(sampletable) <- sampletable$SampleID

#-------------------------------------------------------------------------------
# 读取gff文件
#-------------------------------------------------------------------------------

# library(txdbmaker) # 处理GFF文件
# txdb <- makeTxDbFromGFF(gff, format = "gff3")
# k    <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
# t2g  <- AnnotationDbi::select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")
# t2g <- na.omit(t2g)

library(stringr)
gff_lines <- readLines(gff)
gff_lines <- gff_lines[!startsWith(gff_lines, "#")]

# 只取有转录本特征的行 (mRNA / transcript / *_RNA 等)
df <- read.delim(text = gff_lines, header = FALSE, sep = "\t", quote = "",
                 col.names = c("seqid","source","type","start","end",
                               "score","strand","phase","attr"))
tx_types <- c("mRNA","transcript","lnc_RNA","ncRNA","rRNA","tRNA",
              "snRNA","snoRNA","miRNA","guide_RNA","antisense_RNA","V_gene_segment","C_gene_segment")
tx <- df %>% filter(type %in% tx_types)

# 从 attr 抠出 转录本ID 和 GeneID
t2g <- tx %>%
  mutate(
    TXNAME = str_match(attr, "(?:^|;)ID=(?:rna-)?([^;]+)")[,2],   # 裸转录本ID
    # GENEID = str_match(attr, "GeneID:([0-9]+)")[,2],              # 稳定数字ID
    SYMBOL = str_match(attr, "(?:^|;)gene=([^;]+)")[,2]           # 顺带留 symbol
  ) %>%
  filter(!is.na(TXNAME), !is.na(SYMBOL)) %>%
  select(TXNAME, SYMBOL) %>%
  distinct()

write.table(t2g, file.path(outdir, paste0(Species, "_trans2symbol.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

#-------------------------------------------------------------------------------
# 读取Salmon输出的quant.sf文件，并输出基因表达矩阵
#-------------------------------------------------------------------------------

files <- paste0(indir, sampletable$SampleID, "/quant.sf")
names(files) <- sampletable$SampleID

txi <- tximport(files, type = "salmon", txOut = TRUE)

# 处理 rna- 前缀: salmon ID 带 rna-, txdb 是裸 ID
rownames(txi$abundance) <- sub("^(rna-|gene-)", "", rownames(txi$abundance))
rownames(txi$counts)    <- sub("^(rna-|gene-)", "", rownames(txi$counts))
rownames(txi$length)    <- sub("^(rna-|gene-)", "", rownames(txi$length))


# 汇到基因层，先看基因表达是否有差异
txi_gene <- summarizeToGene(txi, tx2gene = t2g[, c("TXNAME", "SYMBOL")], countsFromAbundance = "lengthScaledTPM")
txi_gene_raw <- summarizeToGene(txi, tx2gene = t2g[, c("TXNAME", "SYMBOL")], countsFromAbundance = "no")

# ---- 导出 gene-level 三个矩阵 ----
out_counts <- data.frame(gene_id = rownames(txi_gene_raw$counts),    txi_gene_raw$counts,    check.names = FALSE)
out_tpm    <- data.frame(gene_id = rownames(txi_gene$abundance), txi_gene$abundance, check.names = FALSE)
out_len    <- data.frame(gene_id = rownames(txi_gene$length),    txi_gene$length,    check.names = FALSE)


write.table(out_counts, file.path(outdir,  paste0(Species,"_gene_counts.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(out_tpm, file.path(outdir,   paste0(Species,"_gene_lengthScaledTPM.tsv")),sep = "\t", quote = FALSE, row.names = FALSE)
write.table(out_len, file.path(outdir,  paste0(Species,"_gene_length.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)







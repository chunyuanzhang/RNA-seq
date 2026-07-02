
suppressMessages({
  library(stringr)
  library(optparse)
  library(dplyr)
})

source("scripts/functions.R")


#-------------------------------------------------------------------------------
# 参数传递
#-------------------------------------------------------------------------------

option_list <- list(
  make_option("--gfffile", type="character", default=NULL, help="gff文件，从中提取信息"),
  make_option("--Species", type="character", default=NULL, help="输出文件，用于后续使用"),
  make_option("--MappingTool", type="character", default=NULL, help="比对工具，影响t2g提取转录本的名字"),
  make_option("--outdir", type="character", default="result/RSEM/", help="输出路径")
)

args <- parse_args(OptionParser(option_list=option_list))

gfffile <- args$gfffile
Species <- args$Species
MappingTool <- args$MappingTool
outdir <- args$outdir
outdir <- check_path(outdir)

# gfffile <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/MyAnalyses/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff"
# Species <- "Chicken"
# outdir <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/MyAnalyses/RNA-seq/RSEM/"



#-------------------------------------------------------------------------------
# 读取gff文件，并提取t2g
#-------------------------------------------------------------------------------


gff_lines <- readLines(gfffile)
gff_lines <- gff_lines[!startsWith(gff_lines, "#")]

# 只取有转录本特征的行 (mRNA / transcript / *_RNA 等)
df <- read.delim(text = gff_lines, header = FALSE, sep = "\t", quote = "",
                 col.names = c("seqid","source","type","start","end",
                               "score","strand","phase","attr"))

tx_types <- c("mRNA","transcript","lnc_RNA","ncRNA","rRNA","tRNA",
              "snRNA","snoRNA","miRNA","guide_RNA","antisense_RNA","V_gene_segment","C_gene_segment")
tx <- df %>% filter(type %in% tx_types)

tx <- tx[grep(tx$attr, pattern = "transcript_id="),]


# 从 attr 抠出 转录本ID 和 GeneID

if(MappingTool == "salmon"){
  t2g <- tx %>%
    mutate(
      TXNAME = str_match(attr, "^ID=([^;]+)")[,2],   # 裸转录本ID
      # GENEID = str_match(attr, "GeneID:([0-9]+)")[,2],              # 稳定数字ID
      SYMBOL = str_match(attr, "(?:^|;)gene=([^;]+)")[,2]           # 顺带留 symbol
    ) %>%
    filter(!is.na(TXNAME), !is.na(SYMBOL)) %>%
    dplyr::select(TXNAME, SYMBOL) %>%
    distinct()
} else if(MappingTool == "rsem"){
  t2g <- tx %>%
    mutate(
      TXNAME = str_match(attr, "Name=([^;]+)")[,2],   # 裸转录本ID
      # GENEID = str_match(attr, "GeneID:([0-9]+)")[,2],              # 稳定数字ID
      SYMBOL = str_match(attr, "(?:^|;)gene=([^;]+)")[,2]           # 顺带留 symbol
    ) %>%
    filter(!is.na(TXNAME), !is.na(SYMBOL)) %>%
    dplyr::select(TXNAME, SYMBOL) %>%
    distinct()
}


#-------------------------------------------------------------------------------
# 输出结果
#-------------------------------------------------------------------------------


outfile <- paste0(outdir, Species, "_trans2symbol.tsv")
write.table(x = t2g, file = outfile, quote = F, row.names = F, sep = "\t")




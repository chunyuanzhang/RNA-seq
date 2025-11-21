
################################################################################
# 拆分用_vs_ 链接的同源基因
################################################################################

extract_gene_from_gene_vs_gene_pair <- function(genes, side){
  
  genes <- genes
  side <- side
  
  sapply(genes, function(gene) {
    if (grepl("_vs_", gene)) {
      # 如果包含 _vs_，则拆分
      parts <- strsplit(gene, "_vs_")[[1]]
      if (side == 1) {
        return(parts[1])
      } else if (side == 2) {
        return(parts[2])
      } else {
        stop("side 只能是 1 或 2, 1 表示提取左边的基因，2 表示提取右边的基因")
      }
    } else {
      # 如果不包含 _vs_，直接返回
      return(gene)
    }
  }, USE.NAMES = FALSE)
  

}


################################################################################
# KEGG and GO
################################################################################

# 由GTF文件中提取基因和蛋白质的对应关系 ========================================

get_Protein2Gene <- function(GTF){
  GTF <- GTF
  
  cat("从GTF中提取蛋白和基因的对应关系\n")

  # 由于emapperannotations_data的注释信息是由蛋白质比对得到的，因此我们需要从蛋白质对应到基因
  GTF <- rtracklayer::import(GTF) %>% as.data.frame()
  gene_to_protein_mapping <- GTF %>% dplyr::select(gene_id, protein_id) %>% filter(!is.na(protein_id)) %>% unique() 
  
  return(gene_to_protein_mapping)
}


# 提取 eggNOG-mapper 比对数据到内存，方便后续从中提取数据

get_eggNOGmapperData <- function(emapperannotations){
  
  emapperannotations <- emapperannotations # gtf比对后从网站上下载的结果
  
  header <- readLines(emapperannotations, n = 5)
  header <- header[5] %>% gsub("#", "", .) %>% strsplit(split = "\t") %>% unlist()
  emapperannotations_data <- read.delim(emapperannotations, comment.char = '#', header = F)
  colnames(emapperannotations_data) <- header
  remove(header)
  
  return(emapperannotations_data)
  
}



# 由emapper 的注释文件中提取 KO 和 gene 的对应关系 =============================

get_KO2Gene <- function(emapperannotations_data, gene_to_protein){
  
  cat("将基因映射为 KO id\n")

  emapperannotations_data <- emapperannotations_data
  gene_to_protein <- gene_to_protein
  
  # KEGG 注释时需要把基因都转换成KO，再使用enricher进行富集分析
  eggNOG_KOtoGene <- emapperannotations_data %>% dplyr::select(KEGG_ko, query) %>% 
    tidyr::separate_rows(KEGG_ko, sep = ",") %>%
    filter(KEGG_ko != "-") %>%
    mutate(KEGG_ko = gsub(pattern = "ko:", replacement = "", KEGG_ko)) %>%
    unique() %>% 
    tidyr::drop_na() %>% 
    dplyr::rename("protein_id" = "query") %>%
    merge(., gene_to_protein, by = "protein_id") %>%
    unique()
  
  return(eggNOG_KOtoGene)
  
}

# 由 emapper 的注释文件提取 GO 和 gene 的对应关系 ==============================

get_GOterm2Gene <- function(emapperannotations_data, gene_to_protein){
  
  cat("将基因映射为 GO term\n")

  emapperannotations_data <- emapperannotations_data
  gene_to_protein<- gene_to_protein
  # GO 注释时
  
  eggNOG_TermtoGene <- emapperannotations_data %>% 
    dplyr::select(GOs, query) %>%
    tidyr::separate_rows(GOs, sep = ",") %>%
    filter(GOs != "-") %>%
    tidyr::drop_na() %>% 
    dplyr::rename("protein_id" = "query") %>%
    merge(., gene_to_protein, by = "protein_id") %>%
    dplyr::select(-protein_id) %>%
    unique() %>%
    dplyr::rename("term" = "GOs", "gene" = "gene_id") %>%
    unique() 
  
  return(eggNOG_TermtoGene)
  
}

# 输入目标差异基因、背景基因、GTF文件和emapperannotations文件，直接得到GO注释结果==

GO_enrichment <- function(diffgenes, backgroundgenes, GTF, emapperannotations){
 
  cat("执行GO 注释 \n")

  GTF <- GTF
  emapperannotations <- emapperannotations
  diffgenes <- diffgenes
  backgroundgenes <- backgroundgenes
  
  #############
  
  gene_to_protein_mapping <- get_Protein2Gene(GTF = GTF)
  emapperannotations_data <- get_eggNOGmapperData(emapperannotations = emapperannotations)
  
  eggNOG_goterm2gene <- get_GOterm2Gene(emapperannotations_data = emapperannotations_data, gene_to_protein = gene_to_protein_mapping)
  
  goterm2name <- AnnotationDbi::select(GO.db::GO.db, keys=unique(eggNOG_goterm2gene$term), 
                                       columns=c("GOID", "TERM"), 
                                       keytype="GOID") 
  colnames(goterm2name) <- c("term", "name")
  
  cat("\tenricher\n")
  enrichment_GO <- clusterProfiler::enricher(gene = diffgenes, 
                                             TERM2GENE = eggNOG_goterm2gene,
                                             TERM2NAME = goterm2name,
                                             pvalueCutoff = 0.05, 
                                             qvalueCutoff = 0.05,
                                             pAdjustMethod = "BH",
                                             universe = backgroundgenes,
  )
  
  enrichment_GO@result <- enrichment_GO@result %>%
    left_join(
      AnnotationDbi::select(GO.db::GO.db,
             keys = .$ID,
             columns = c("GOID", "ONTOLOGY"),
             keytype = "GOID"),
      by = c("ID" = "GOID")
    )
  
  return(enrichment_GO)
}


# 输入目标差异基因、背景基因、GTF文件和emapperannotations文件，直接得到KEGG注释结果=

KEGG_enrichment <- function(diffgenes, backgroundgenes, GTF, emapperannotations){
  
  cat("执行KEGG注释\n")

  options(timeout = 300) # 增加时长，避免因网速不够导致的错误推出

  diffgenes <- diffgenes
  backgroundgenes <- backgroundgenes
  GTF <- GTF
  emapperannotations <- emapperannotations
  
  gene_to_protein_mapping <- get_Protein2Gene(GTF = GTF)
  emapperannotations_data <- get_eggNOGmapperData(emapperannotations = emapperannotations)
  
  eggNOG_kegg2gene <- get_KO2Gene(emapperannotations_data = emapperannotations_data, gene_to_protein = gene_to_protein_mapping)
   
  diffgenes2kegg <- eggNOG_kegg2gene %>% dplyr::filter(gene_id %in% diffgenes) %>% pull(KEGG_ko) 
  backgroundgenes2kegg <- eggNOG_kegg2gene %>%  dplyr::filter(gene_id %in% backgroundgenes) %>% pull(KEGG_ko) 
 
  cat("\tenricherKEGG\n")
  enrichment_kegg <- clusterProfiler::enrichKEGG(diffgenes2kegg,
                                organism = "ko",
                                keyType = "kegg",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                pAdjustMethod = "BH",
                                universe = backgroundgenes2kegg,
                                use_internal_data = FALSE)
  
  enrichment_kegg@result <- enrichment_kegg@result %>%
    rowwise() %>%
    mutate(
      # 添加 gene_ids
      gene_ids = {
        genekeggs <- strsplit(geneID, "/")[[1]]
        gene_ids <- eggNOG_kegg2gene[eggNOG_kegg2gene$KEGG_ko %in% genekeggs, ]$gene_id %>% 
          unique()
        gene_ids[gene_ids %in% diffgenes] %>% 
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
  
   return(enrichment_kegg)

}

# 绘制点图

anno_plot <- function(object, title){
  
  object <- object
  title <- title
  
  p <- enrichplot::dotplot(object, showCategory=20, font.size=10, label_format=70) + 
    scale_size_continuous(range = c(1,7)) + 
    theme_minimal() + 
    ggtitle(title)
  
  return(p)
  
}

################################################################################
# WGCNA
################################################################################


# 从 dds 文件开始，准备 WGCNA 分析的输入数据 ===================================


# prepare_data_for_wgcna <- function(dds, min_count = 10, min_samples = 0.1) {
  
#   dds <- dds
#   min_count <- min_count
#   min_samples <- min_samples
  
  
#   # 获取标准化的表达数据
#   # 使用rlog或vst转换以获得同方差稳定的数据
#   if (ncol(dds) < 30) {
#     expr_data <- rlog(dds, blind = FALSE)
#     cat("样本数 < 30，使用rlog转换\n")
#   } else {
#     expr_data <- vst(dds, blind = FALSE)
#     cat("样本数 >= 30，使用VST转换\n")
#   }
  
#   # 提取表达矩阵
#   expr_matrix <- assay(expr_data)
  
#   # 过滤低表达基因
#   # 保留在至少10%样本中表达量 > min_count的基因
#   min_samples_count <- ceiling(ncol(expr_matrix) * min_samples)
#   keep_genes <- rowSums(expr_matrix > min_count) >= min_samples_count
#   expr_matrix <- expr_matrix[keep_genes, ]
  
#   cat(sprintf("过滤后保留 %d 个基因\n", nrow(expr_matrix)))
  
#   # 转置矩阵（WGCNA需要样本为行，基因为列）
#   datExpr <- t(expr_matrix)
  
#   # 检查数据质量
#   gsg <- goodSamplesGenes(datExpr, verbose = 3)
#   if (!gsg$allOK) {
#     if (sum(!gsg$goodGenes) > 0) {
#       cat(sprintf("移除 %d 个质量不佳的基因\n", sum(!gsg$goodGenes)))
#     }
#     if (sum(!gsg$goodSamples) > 0) {
#       cat(sprintf("移除 %d 个质量不佳的样本\n", sum(!gsg$goodSamples)))
#     }
#     datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
#   }
  
#   return(list(
#     datExpr = datExpr,
#     sample_info = colData(dds)[rownames(datExpr), ],
#     original_dds = dds
#   ))
# }

prepare_data_for_wgcna <- function(dds, quantile, min_count = 10, min_samples = 0.1 ) {
  
  dds <- dds
  quantile <- quantile
  min_count <- min_count
  min_samples <- min_samples
  
  ### vsd处理数据
  vsd_data <- DESeq2::getVarianceStabilizedData(dds)
  
  # 保留在至少 10 %样本中表达量 > min_count的 基因
  min_samples_count <- ceiling(ncol(vsd_data) * min_samples)
  keep_genes <- rowSums(vsd_data > min_count) >= min_samples_count
  vsd_data <- vsd_data[keep_genes, ]
  
  
  ### 过滤掉低方差的基因
  rv <- rowVars(vsd_data)
  quantiles <- quantile(rv, quantile)
  filtered_data <- vsd_data[rv > quantiles, ]
  
  cat("剩余 ", dim(filtered_data)[1], " 个基因和 ", dim(filtered_data)[2], "个样本，样本示例： ", colnames(filtered_data)[1:5]  )

  return(filtered_data)
  
}


#  自动选择软阈值 ==============================================================


auto_select_soft_Power <- function(datExpr){
  datExpr <- datExpr
  
  # 计算不同阈值下的网络拓扑特性
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
  sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  
  # 绘制软阈值选择图
  # pdf("soft_threshold_selection.pdf", width = 12, height = 5)
  par(mfrow = c(1, 2))
  
  # 无标度拟合指数
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit,signed R^2",
       type = "n", main = "Scale independence")
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = 0.9, col = "red")
  abline(h = 0.90, col = "red")
  
  # 平均连接度
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity", type = "n",
       main = "Mean connectivity")
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = "red")
  # dev.off()
  
  # 自动选择软阈值
  if (is.na(sft$powerEstimate)) {
    softPower <- 6  # 默认值
    cat("无法自动确定软阈值，使用默认值 6\n")
  } else {
    softPower <- sft$powerEstimate
    cat(sprintf("选择的软阈值: %d\n", softPower))
  }
  
  return(softPower)
  
  
}


# WGCNA网络分析，并绘制图形

net_and_plot <- function(datExpr, softPower, labeltext = "Module colors"){
  
  datExpr <- datExpr
  softPower <- softPower
  
  net <- WGCNA::blockwiseModules(t(datExpr), 
                               power = softPower,
                               TOMType = "unsigned", 
                               minModuleSize = 30,
                               reassignThreshold = 1e-6,   # 初步未分配的模块是否重新分配，0表示不分配, 1e-6 默认值，几乎所有基因都会被分配
                               mergeCutHeight = 0.25,
                               numericLabels = TRUE, 
                               pamRespectsDendro = FALSE,
                               saveTOMs = TRUE,
                               saveTOMFileBase = "TOM",
                               verbose = 3
  )
  
  
  modulecolors <- WGCNA::labels2colors(net$colors)
  
  plotDendroAndColors(net$dendrograms[[1]], 
                      modulecolors[net$blockGenes[[1]]],
                      labeltext,
                      dendroLabels = FALSE, 
                      hang = 0.03,
                      addGuide = TRUE, 
                      guideHang = 0.05
                      )

  return(net)
}








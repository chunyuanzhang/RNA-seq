


################################################################################
# KEGG and GO
################################################################################

# 由GTF文件中提取基因和蛋白质的对应关系 ========================================

get_Protein2Gene <- function(GTF){
  GTF <- GTF
  
  # 由于emapperannotations_data的注释信息是由蛋白质比对得到的，因此我们需要从蛋白质对应到基因
  
  gene_to_protein_mapping <- GTF %>% dplyr::select(gene_id, protein_id) %>% filter(!is.na(protein_id)) %>% unique() 
  
  return(gene_to_protein_mapping)
}


# 由emapper 的注释文件中提取 KO 和 gene 的对应关系 =============================

get_KO2Gene <- function(emapperannotations_data, gene_to_protein){
  
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
    dplyr::rename("term" = "GOs", "gene" = "gene_id") %>%
    unique() 
  
  return(eggNOG_TermtoGene)
  
}









################################################################################
# WGCNA
################################################################################


# 从 dds 文件开始，准备 WGCNA 分析的输入数据 ===================================


prepare_data_for_wgcna <- function(dds, min_count = 10, min_samples = 0.1) {
  
  dds <- dds
  min_count <- min_count
  min_samples <- min_samples
  
  
  # 获取标准化的表达数据
  # 使用rlog或vst转换以获得同方差稳定的数据
  if (ncol(dds) < 30) {
    expr_data <- rlog(dds, blind = FALSE)
    cat("样本数 < 30，使用rlog转换\n")
  } else {
    expr_data <- vst(dds, blind = FALSE)
    cat("样本数 >= 30，使用VST转换\n")
  }
  
  # 提取表达矩阵
  expr_matrix <- assay(expr_data)
  
  # 过滤低表达基因
  # 保留在至少10%样本中表达量 > min_count的基因
  min_samples_count <- ceiling(ncol(expr_matrix) * min_samples)
  keep_genes <- rowSums(expr_matrix > min_count) >= min_samples_count
  expr_matrix <- expr_matrix[keep_genes, ]
  
  cat(sprintf("过滤后保留 %d 个基因\n", nrow(expr_matrix)))
  
  # 转置矩阵（WGCNA需要样本为行，基因为列）
  datExpr <- t(expr_matrix)
  
  # 检查数据质量
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0) {
      cat(sprintf("移除 %d 个质量不佳的基因\n", sum(!gsg$goodGenes)))
    }
    if (sum(!gsg$goodSamples) > 0) {
      cat(sprintf("移除 %d 个质量不佳的样本\n", sum(!gsg$goodSamples)))
    }
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  }
  
  return(list(
    datExpr = datExpr,
    sample_info = colData(dds)[rownames(datExpr), ],
    original_dds = dds
  ))
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









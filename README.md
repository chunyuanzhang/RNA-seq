本流程参考资料

>https://github.com/CebolaLab/RNA-seq?tab=readme-ov-file
https://github.com/hbctraining/DGE_workshop_salmon/blob/master/lessons/02_DGE_count_normalization.md
https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html



>可通过网页使用服务器上配的 R studio
>


当前只计算基因的差异表达，并没有统计转录本的差异表达，还需要把转录本的差异表达也写出来



# 转录组分析

## 转录组分析步骤：
1. 质量控制，获得CleanData
1. 数据比对，并对比对结果进行质量评估
1. 比对结果标准化/归一化【样本内不同基因的比较，种内群体间的比较，种间比较】
1. 组间差异分析

## 标准化/归一化如何选择



# 种内和种间的差异分析
因 DEseq2 通过 design 设置差异比较分组，因此我们通过设置 design 来控制比较组

infotable.csv中包含Species列，如果 design 选择了 Species，流程就会进入物种间的分析，否则就是在种内分析

## 种内的 RNA-seq 分析
- RSEM统计表达量
- 导入到DESeq中，按照默认参数进行差异分析


## 跨物种的 RNA-seq 分析

目前的思路：
- 种内标准化【通过分别计算 effecter size 实现种内标准化】
- 提取1:1的基因【】
- 差异分析
- KEGG和GO注释
    - 使用差异的1:1同源基因作为感兴趣的基因
    - 使用所有表达的1:1同源基因作为背景基因【是否有必要更换为所有表达的基因，而不局限在同源基因？】


目前已经完成使用自定义的KEGG注释流程，还需要增加GO注释
目前多组比较是通过在R内循环实现的，最好改成snakemake自动控制，不再在R中循环，不好控制输出结果
尽量让种间比较和种内比较保持一致，目前种内和种间的代码是有差别的
比较数据库和使用数据库直接注释的KEGG结果存在差异，看起来自定义的结果更多，目前保留

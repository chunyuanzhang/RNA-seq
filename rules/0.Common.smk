
import pandas as pd
import numpy as np
import yaml
import json
import os
import glob

# === 工具或分析流程选择 ===============================
pipeline = "onestep" # 使用rsem工具一步完成比对和统计

print(pipeline)
# === 指定分析步骤 ====================================
step = "map"


# === 元信息 =========================================
infotable = "metafile"
untreated = "mKO"
lfc = 1
pval = 0.05
# infotable文件内容示意，为了代码重复使用方便，infotable不要列
# 第一列是样本名，第二列是组名
# White_1W_F1.sra,White
# White_1W_F2.sra,White
# White_1W_F3.sra,White
# Wild_1W_F1.sra,Wild
# Wild_1W_F2.sra,Wild
# Wild_1W_F3.sra,Wild


# === 指定参考基因组 ==================================

## 参考基因组建立索引
### 在基因组所在文件夹执行下行命令 双端150bp测序时, --sjdbOverhang 149
### STAR --runThreadN 64 --runMode genomeGenerate --genomeDir ./STARindex --genomeFastaFiles ./STARindex/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.fna --sjdbGTFfile  ./GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.gtf --sjdbOverhang 149 --genomeSAindexNbases 13
## rsem也需要创建参考基因组

### ~/tools/rsem/bin/rsem-prepare-reference --gtf GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.gtf --star GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.fna RSEMindex/IASCAAS_PekinDuck_T2T -p 64 

## 创建transcript.fa文件
### 本流程使用*_Aligned.toTranscriptome.out.bam文件提取表达量
### 由于该文件是使用转录本进行定位的，因此需要根据参考基因组和gtf文件制作转录本构成的新参考基因组
### gffread -w GCF_047663525.1_IASCAAS_PekinDuck_T2T_transcripts.fna -g GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.fna GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.gff 
### sed 's/>rna-/>/'  GCF_047663525.1_IASCAAS_PekinDuck_T2T_transcripts.fna -i 

# === 外部参数处理 ===================================

if os.path.exists(infotable):
    metainfo = pd.read_csv(infotable, sep = ",", dtype=str, index_col = False, header=None, skiprows=1)
    metainfo.dropna(how='all', inplace=True) # 删除可能存在的空行
    metainfo.columns = [f'V{i+1}' for i in range(len(metainfo.columns))]
    metainfo = metainfo[ ~np.array([v.startswith("#") for v in metainfo.V1.to_list()])]  # 删除被#注释掉的行
    # metainfo.index = metainfo.V1
    samples = metainfo.V1.tolist()
    # 
    metainfo_dict_Group = dict(zip(metainfo.V1, metainfo.V2)) 
    # 
    metainfo_dict_Genome = dict(zip(metainfo.V1, metainfo.V3)) 
    referenceDir = "/home/zhangchunyuan/zhangchunyuan/reference/"

    

    # 设置wildcard约束
    wildcard_constraints:
        samples = "|".join(samples)
    




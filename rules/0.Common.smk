
import pandas as pd
import numpy as np
import yaml
import json
import os


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
## 本流程使用STAR进行比对，注意在比对前需要先手动对参考基因组建立索引
genomeDir = "/home/zhangchunyuan/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b_Ensembl/STARindex"
ref_fa = "/home/zhangchunyuan/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b_Ensembl/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa"
ref_gtf = "/home/zhangchunyuan/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b_Ensembl/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gtf"
ref_transscripts = "/home/zhangchunyuan/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b_Ensembl/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.transcripts.fa"

## 在基因组所在文件夹执行下行命令 双端150bp测序时, --sjdbOverhang 149
## STAR --runThreadN 64 --runMode genomeGenerate --genomeDir ./STARindex --genomeFastaFiles ./STARindex/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.fna --sjdbGTFfile  ./GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.gtf --sjdbOverhang 149 --genomeSAindexNbases 13



# === 外部参数处理 ===================================

if os.path.exists(infotable):
    metainfo = pd.read_csv(infotable, sep = ",", dtype=str, index_col = False, header=None, skiprows=1)
    metainfo.dropna(how='all', inplace=True) # 删除可能存在的空行
    metainfo.columns = [f'V{i+1}' for i in range(len(metainfo.columns))]
    metainfo = metainfo[ ~np.array([v.startswith("#") for v in metainfo.V1.to_list()])]  # 删除被#注释掉的行
    # metainfo.index = metainfo.V1

    metainfo_dict = dict(zip(metainfo.V1, metainfo.V2)) 
    
    samples = metainfo.V1.tolist()

    # 设置wildcard约束
    wildcard_constraints:
        samples = "|".join(samples)
    




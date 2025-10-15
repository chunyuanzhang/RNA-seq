
import pandas as pd
import numpy as np
import yaml
import json
import os
import glob


# === 工具或分析流程选择 ===============================
pipeline = "onestep" # 使用rsem工具一步完成比对和统计

# === 指定分析步骤 ====================================
step = "all"


# === 元信息 =========================================
infotable = "infotable.csv"
untreated = "Chicken"
lfc = 1
pval = 0.05
Orthologgenes = "/home/zhangchunyuan/zhangchunyuan/Orthologgenes/Chicken_vs_Duck.one_to_one.orthologenes.txt"

# infotable文件内容示意，为了代码重复使用方便，infotable不要列
# 第一列是样本名，第二列是组名，第三列为参考基因组名称【文件夹名称】
# 若还有其他信息，可以从第四列开始继续添加
# White_1W_F1.sra,White,IASCAAS_PekinDuck_T2T
# White_1W_F2.sra,White,IASCAAS_PekinDuck_T2T
# White_1W_F3.sra,White,IASCAAS_PekinDuck_T2T
# Wild_1W_F1.sra,Wild,IASCAAS_PekinDuck_T2T
# Wild_1W_F2.sra,Wild,IASCAAS_PekinDuck_T2T
# Wild_1W_F3.sra,Wild,IASCAAS_PekinDuck_T2T


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

# === 指定参考基因组所在路径 ===========================
referenceDir = "/home/zhangchunyuan/zhangchunyuan/reference/"
### 我们所有的参考基因组都在该路径下，只需要指定参考基因组的名字即可

# === 外部参数处理 ===================================

if os.path.exists(infotable):
    metainfo = pd.read_csv(infotable, sep = ",", dtype=str, index_col = False, header=None, skiprows=1)
    metainfo.dropna(how='all', inplace=True) # 删除可能存在的空行
    metainfo.columns = [f'V{i+1}' for i in range(len(metainfo.columns))]
    metainfo = metainfo[ ~np.array([v.startswith("#") for v in metainfo.V1.to_list()])]  # 删除被#注释掉的行
    # metainfo.index = metainfo.V1
    samples = metainfo.V1.tolist()
    # 第一列存储样本名称，第二列存储分组名称，第三列存储使用的参考基因组
    metainfo_dict_Group = dict(zip(metainfo.V1, metainfo.V2)) 
    metainfo_dict_Genome = dict(zip(metainfo.V1, metainfo.V3)) 

    # infotable不一定有V4、V5这些列，需要进行判断，如果有，则存储到字典中
    if len(metainfo.columns) > 3:
        additional_dicts = {}
        for v in range(3, len(metainfo.columns)):
            col_name = metainfo.columns[v]
            dict_name = f"metainfo_dict_{col_name}"
            additional_dicts[dict_name] = dict(zip(metainfo.V1, metainfo[col_name]))
    
    # infotable 的第三列是参考基因组，如果参考基因组不相同，那么归一化过程就需要在同一个参考基因组内部进行，归一化后再进行组间的差异比较
    references = set(metainfo_dict_Genome.values())
    if len(references) > 1:
        interspecies = True
        species = metainfo.V4.drop_duplicates().to_list()
        print(species)
        # 第四列存储物种信息
        metainfo_dict_Species_to_Genome = dict(zip(metainfo.V4, metainfo.V3)) 
    else:
        interspecies = False


    # 设置wildcard约束
    wildcard_constraints:
        samples = "|".join(samples)
    





import pandas as pd
import numpy as np
import yaml
import json
import os
import glob

### load config file======================================================================================================
configfile: "config/config.yaml"

# === 工具或分析流程选择 ===============================
 # 使用rsem工具一步完成比对和统计
pipeline = config["pipeline"]

# === 指定分析步骤 ====================================
step = config["step"]

# === 元信息 =========================================
infotable = config["infotable"]
design = config["design"]
untreated = config["untreated"]
lfc = config["lfc"]
pval = config["pval"]


# infotable文件内容示意，列名是固定的，不可更改
# 第一列是样本名，第二列是组名，第三列为参考基因组名称【文件夹名称】
# SampleID,GroupID,GenomeName,Species
# FHSY-C-31,FHSY,IASCAAS_PekinDuck_T2T,Duck
# FHSY-C-32,FHSY,IASCAAS_PekinDuck_T2T,Duck
# FHSY-C-33,FHSY,IASCAAS_PekinDuck_T2T,Duck
# FHSY-C-34,FHSY,IASCAAS_PekinDuck_T2T,Duck
# FHSY-C-35,FHSY,IASCAAS_PekinDuck_T2T,Duck
# XSJ-C-31,XSJ,bGalGal1_mat_broiler_GRCg7b,Chicken
# XSJ-C-32,XSJ,bGalGal1_mat_broiler_GRCg7b,Chicken
# XSJ-C-33,XSJ,bGalGal1_mat_broiler_GRCg7b,Chicken
# XSJ-C-34,XSJ,bGalGal1_mat_broiler_GRCg7b,Chicken
# XSJ-C-35,XSJ,bGalGal1_mat_broiler_GRCg7b,Chicken

# === 指定参考基因组所在路径 ===========================
referenceDir = config["referenceDir"]
### 我们所有的参考基因组都在该路径下，只需要指定参考基因组的名字即可

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
    metainfo = pd.read_csv(infotable, sep = ",", dtype=str)
    metainfo.dropna(how='all', inplace=True) # 删除可能存在的空行
    # metainfo.columns = [f'V{i+1}' for i in range(len(metainfo.columns))]
    metainfo = metainfo[ ~np.array([v.startswith("#") for v in metainfo.SampleID.to_list()])]  # 删除被#注释掉的行
    # metainfo.index = metainfo.V1
    samples = metainfo.SampleID.tolist()
    # 第一列存储样本名称，第二列存储分组名称，第三列存储使用的参考基因组，第四列指定物种，将所有信息分别存储在字典中，方便随时调用
    metainfo_dict_Group = dict(zip(metainfo.SampleID, metainfo.GroupID)) 
    metainfo_dict_Genome = dict(zip(metainfo.SampleID, metainfo.GenomeName)) 
    metainfo_dict_Species = dict(zip(metainfo.SampleID, metainfo.Species)) 
    
    # 为了让我们的代码可以适应多组间比较的情况
    loops = metainfo[design].values
    loops = list(set(loops))
    treated = [x for x in loops if x not in untreated]
    pairnames = [f"{design}_{t}_vs_{untreated}" for t in treated]
    # print(pairnames)
    
    # 如果在物种间进行比较，需要把物种和参考基因组去冗余，便于控制命令重复次
    if design == "Species":
        # 将物种和参考基因组的对应关系放在字典中
        metainfo_dict_Species_to_Genome = dict(zip(metainfo.Species, metainfo.GenomeName)) 
        # 物种间差异分析需要提供1:1的同源基因
        Orthologgenes = config["Orthologgenes"]

    # 设置wildcard约束
    wildcard_constraints:
        samples = "|".join(samples)
    

    


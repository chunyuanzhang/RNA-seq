
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
    metainfo = metainfo[ ~np.array([v.startswith("#") for v in metainfo.SampleID.to_list()])]  # 删除被#注释掉的行

    # 验证 design 列是否存在
    if design not in metainfo.columns:
        raise ValueError(f"Design column '{design}' not found in {infotable}. Available columns: {metainfo.columns.tolist()}")

    samples = metainfo.SampleID.tolist()
    
    metainfo_dict_Group = dict(zip(metainfo.SampleID, metainfo.GroupID)) 
    metainfo_dict_Genome = dict(zip(metainfo.SampleID, metainfo.GenomeName)) 
    metainfo_dict_Species = dict(zip(metainfo.SampleID, metainfo.Species)) 
    # metainfo_dict_Species_to_Genome = dict(zip(metainfo.Species, metainfo.GenomeName)) 
    metainfo_dict_Design_to_Genome = dict(zip(metainfo[design], metainfo.GenomeName))

    # === 成对比较逻辑（根据 design 列进行分组）==========
    all_design_values = metainfo[design].unique().tolist()

    treated = [x for x in all_design_values if x != untreated]
    pairnames = [f"{design}_{t}_vs_{untreated}" for t in treated]
    
    # 为每个 pairname 定义其包含的两个分组
    pairname_groups = {}
    for t in treated:
        pairname = f"{design}_{t}_vs_{untreated}"
        pairname_groups[pairname] = [t, untreated]

    # 生成用于 expand 的组合列表
    gokegg_outputs = []
    for pairname in pairnames:
        for group in pairname_groups[pairname]:
            gokegg_outputs.append(f"result/DEseq/{pairname}.{group}.KEGG.csv")
            gokegg_outputs.append(f"result/DEseq/{pairname}.{group}.GO.csv")

    #print(gokegg_outputs)

    # 如果在物种间进行比较，需要把物种和参考基因组去冗余，便于控制命令重复次
    if design == "Species":
        # 物种间差异分析需要提供1:1的同源基因
        Orthologgenes = config["Orthologgenes"]

    # 设置wildcard约束
    wildcard_constraints:
        samples = "|".join(samples)


else:
    raise FileNotFoundError(f"Infotable file not found: {infotable}")


# === 辅助函数 =========================================

def get_gtf_by_design_value(wildcards):
    """
    根据 design 列的值（可能是 GroupID 或 Species）获取对应的 GTF 文件
    wildcards.loop 的值来自 design 列
    """
    genome = metainfo_dict_Design_to_Genome.get(wildcards.loop)
    gtf_files = glob.glob(f"{referenceDir}{genome}/*.gtf")
    if not gtf_files:
        raise FileNotFoundError(f"No GTF file found in {referenceDir}{genome}/")
    return gtf_files[0]


def get_emapper_by_design_value(wildcards):
    """
    根据 design 列的值（可能是 GroupID 或 Species）获取对应的 emapper.annotations 文件
    """
    genome = metainfo_dict_Design_to_Genome.get(wildcards.loop)
    emapper_files = glob.glob(f"{referenceDir}{genome}/*.emapper.annotations")
    if not emapper_files:
        raise FileNotFoundError(f"No emapper.annotations file found in {referenceDir}{genome}/")
    return emapper_files[0]

def get_genome_gtf(wildcards):
    """根据样本获取对应的基因组 GTF 文件"""
    genome = metainfo_dict_Genome.get(wildcards.sample)
    gtf_files = glob.glob(f"{referenceDir}{genome}/*.gtf")
    if not gtf_files:
        raise FileNotFoundError(f"No GTF file found in {referenceDir}{genome}/")
    return gtf_files[0]




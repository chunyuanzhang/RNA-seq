#!/usr/bin/env bash
# ============================================================
# build_salmon_index.sh
# 一次性建立 decoy-aware salmon 索引 (鸡 + 鸭)
# 方案 A: selective alignment, 基因组作 decoy
#
# 用法:
#   bash build_salmon_index.sh
#
# 前提:
#   - salmon 已安装且版本 >= 0.14 (推荐 >= 1.5), 在 PATH 中或改下方 SALMON 变量
#   - 鸡鸭的 转录组 FASTA + 基因组 FASTA 均已具备
# ============================================================

set -euo pipefail   # 出错即停 / 未定义变量报错 / 管道错误传播

# ------------------------------------------------------------
# 0. 全局参数 (按需修改)
# ------------------------------------------------------------
SALMON="salmon"          # salmon 可执行文件; 若用 conda 环境, 先 conda activate 再跑
THREADS=16               # 建索引线程数
KMER=31                  # k-mer 大小 (reads >= 75bp 用 31 即可)

# 索引输出根目录 (每物种一个子目录), 与 config.yaml 的 salmon_index_root 保持一致
INDEX_ROOT="SALMONIndex"

# ------------------------------------------------------------
# 1. 每物种的参考文件路径
#    transcripts: gffread 生成的转录组 FASTA (cDNA)
#    genome:      基因组 FASTA (用作 decoy)
#    gtf:         基因组 GTF 
# ------------------------------------------------------------

# ----- 鸭 (Pekin Duck T2T) -----
GENOME="GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.fna"
GFF="GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.gff"
TRANSCRIPTS="GCF_047663525.1_IASCAAS_PekinDuck_T2T_forSalmon_transcripts.fa"

echo "生成 ${TRANSCRIPTS} 文件"
module load gffread
gffread -w "${TRANSCRIPTS}" -g "${GENOME}" "${GFF}"

# ============================================================
# 以下逻辑无需修改
# ============================================================

# ------------------------------------------------------------
# 建索引函数: build_one <物种名> <转录组FASTA> <基因组FASTA>
# ------------------------------------------------------------
build_one() {
    local sp="$1"
    local transcripts="$2"
    local genome="$3"
    local idxdir="${INDEX_ROOT}/${sp}"

    echo ""
    echo "==================================================="
    echo "==== 物种: ${sp}"
    echo "==================================================="

    # --- 检查输入文件存在 ---
    if [[ ! -f "${transcripts}" ]]; then
        echo "[错误] 转录组 FASTA 不存在: ${transcripts}" >&2; exit 1
    fi
    if [[ ! -f "${genome}" ]]; then
        echo "[错误] 基因组 FASTA 不存在: ${genome}" >&2; exit 1
    fi

    mkdir -p "${idxdir}"

    # --- 若已建好则跳过 (seq.bin 是完成标志之一) ---
    if [[ -f "${idxdir}/seq.bin" ]]; then
        echo "[跳过] ${sp} 索引已存在: ${idxdir} (如需重建请先删除该目录)"
        return 0
    fi

    # --- Step 1: decoy 文件 (基因组序列名) ---
    echo "[1/3] 生成 decoys.txt ..."
    grep '^>' "${genome}" | sed 's/^>//; s/[[:space:]].*//' > "${idxdir}/decoys.txt"
    echo "      decoy 数量: $(wc -l < "${idxdir}/decoys.txt")"

    # --- Step 2: gentrome (转录组在前, 基因组在后), 自动兼容 .gz ---
    echo "[2/3] 拼接 gentrome.fna.gz ..."
    cat_any() { case "$1" in *.gz) zcat "$1";; *) cat "$1";; esac; }
    ( cat_any "${transcripts}"; cat_any "${genome}" ) | gzip -c > "${idxdir}/gentrome.fna.gz"

    # --- Step 3: salmon index ---
    echo "[3/3] 构建 salmon 索引 (decoy-aware) ..."
    "${SALMON}" index \
        -t "${idxdir}/gentrome.fna.gz" \
        -d "${idxdir}/decoys.txt" \
        -i "${idxdir}" \
        -k "${KMER}" \
        -p "${THREADS}"

    echo "[完成] ${sp} 索引: ${idxdir}"
}

# ------------------------------------------------------------
# 主流程
# ------------------------------------------------------------
module load salmon
echo "salmon 版本:"
"${SALMON}" --version || { echo "[错误] 找不到 salmon, 请检查 SALMON 变量或 conda 环境" >&2; exit 1; }

mkdir -p "${INDEX_ROOT}"

build_one "duck"    "${TRANSCRIPTS}"    "${GENOME}"

echo ""
echo "==================================================="
echo "全部索引构建完成。"
echo "  ${INDEX_ROOT}/duck"
echo ""
echo "==================================================="

# eccToolkit 快速入门指南

## 安装

```bash
# 1. 解压并进入目录
tar -xzf eccToolkit-0.6.0.tar.gz
cd eccToolkit-0.6.0

# 2. 创建 conda 环境 (推荐)
conda create -n ecctoolkit python=3.11
conda activate ecctoolkit

# 3. 安装依赖
pip install -e .

# 4. 验证安装
ecc --version    # 应显示: eccToolkit, version 0.6.0
ecc --help       # 显示所有可用命令
```

## 依赖要求

- Python >= 3.8
- minimap2 (用于 sim-region)
- bedtools (用于富集分析，可选)

```bash
# 安装 minimap2
conda install -c bioconda minimap2

# 安装 bedtools (可选)
conda install -c bioconda bedtools
```

## 核心入口：sim-region

`sim-region` 生成标准 eccDNA 区域输出（`.all.bed`/`.all.fa`），这是后续
`sim-reads`、检测与分析的统一入口。

典型链路：

```
sim-region → sim-reads → detection/enrichment/hotspot/TE
```

## 快速测试

```bash
# 运行单元测试
pytest tests/ -v -m "not slow"

# 测试 sim-region (需要 minimap2)
ecc sim-region \
    -r tests/data/test_reference.fa \
    -o test_output \
    -u 50 -t 4 --seed 42

# 测试 sim-reads
ecc sim-reads \
    -i test_output/test_output.all.fa \
    -o test_reads \
    --cov-ngs 500 \
    --skip-hifi --skip-ont
```

## 配置模板

- `configs/sim-region.yaml`: sim-region 参数清单（CLI 对照）
- `configs/sim-reads.yaml`: sim-reads YAML 配置模板（可直接 `--config` 使用）

常用高级参数（详见 `ecc sim-region --help`）：
`--min-coverage`, `--no-hit-policy`, `--split-length`, `--multiplier-u`, `--multiplier-m`。

示例：

```bash
ecc sim-reads \
    -i test_output/test_output.all.fa \
    -o test_reads \
    --config configs/sim-reads.yaml
```

## 完整工作流示例

### 1. 模拟 eccDNA 区域

```bash
ecc sim-region \
    -r /path/to/reference.fa \
    -o my_simulation \
    -u 1000 \          # 1000 Unique eccDNA
    -m 100 \           # 100 Multi-mapped eccDNA
    -c 50 \            # 50 Chimeric eccDNA
    -t 8 \             # 8 线程
    --seed 42          # 随机种子 (可重复)

# 输出文件:
# my_simulation/my_simulation.all.bed     - 所有 eccDNA 区域
# my_simulation/my_simulation.all.fa      - 所有 eccDNA 序列
# my_simulation/my_simulation.unique.fa   - 仅 Unique eccDNA
# my_simulation/my_simulation.multi.fa    - 仅 Multi-mapped eccDNA
# my_simulation/my_simulation.chimeric.fa - 仅 Chimeric eccDNA
# my_simulation/my_simulation.qc.log      - QC 报告
```

### 2. 模拟测序读段 (RCA)

```bash
# 按读段数量指定
ecc sim-reads \
    -i my_simulation/my_simulation.all.fa \
    -o reads_output \
    --config configs/sim-reads.yaml \
    --cov-ngs 10000 \    # 10000 条 NGS reads
    --cov-hifi 1000 \    # 1000 条 HiFi reads
    --cov-ont 1000       # 1000 条 ONT reads

# 按覆盖度指定
ecc sim-reads \
    -i my_simulation/my_simulation.all.fa \
    -o reads_output \
    --output-mode coverage \
    --cov-ngs 50 \       # 50x NGS 覆盖度
    --cov-hifi 30        # 30x HiFi 覆盖度

# 输出文件:
# reads_output/reads_ngs_R1.fastq  - NGS 双端 R1
# reads_output/reads_ngs_R2.fastq  - NGS 双端 R2
# reads_output/reads_hifi.fastq    - PacBio HiFi
# reads_output/reads_ont.fastq     - Oxford Nanopore
# reads_output/truth.tsv           - Ground truth
```

### 3. 数据处理

```bash
# 合并多个 CSV 文件
ecc merge -i results_dir/ -o merged.csv

# 过滤低质量数据
ecc filter -i merged.csv -o filtered.csv --min-percent 80

# 生成汇总报告
ecc report -i filtered.csv -o summary
```

## 可用命令列表

| 类别 | 命令 | 状态 | 说明 |
|------|------|------|------|
| **模拟** | `sim-region` | ✅ 核心入口 | 模拟 eccDNA 区域 |
| | `sim-reads` | ✅ 可用 | 模拟测序读段 |
| **处理** | `merge` | ✅ 可用 | 合并 CSV 文件 |
| | `parse` | ✅ 可用 | 解析 eccDNA CSV |
| | `filter` | ✅ 可用 | 过滤数据 |
| | `report` | ✅ 可用 | 生成报告 |
| **检测** | `circlemap` | 🔧 待迁移 | Circle-Map 流程 |
| | `validate` | 🔧 待迁移 | 验证 eccDNA |
| **富集** | `enrich-cnv` | 🔧 待迁移 | CNV 富集分析 |
| | `enrich-tad` | 🔧 待迁移 | TAD 边界富集 |
| **热点** | `hotspot-detect` | 🔧 待迁移 | 热点检测 |
| **TE** | `te-analyze` | 🔧 待迁移 | TE 组成分析 |

## 获取帮助

```bash
# 查看所有命令
ecc --help

# 查看特定命令帮助
ecc sim-region --help
ecc sim-reads --help
```

## 问题反馈

GitHub Issues: https://github.com/yourusername/eccToolkit/issues

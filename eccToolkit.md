# eccToolkit: 一套用于eccDNA分析的综合工具集



## 项目简介 (Project Overview)

**eccToolkit** 是一系列为了研究环状DNA (eccDNA) 而开发的生物信息学分析脚本。本工具集涵盖了从处理原始测序数据、鉴定eccDNA、寻找断点重复序列、进行饱和度分析到最终结果的统计与整合等多个方面。

这些脚本旨在提供一套可重复、可定制的分析流程，以支持eccDNA相关的科学研究。



## 目录 (Table of Contents)

- [安装指南](https://www.google.com/search?q=%23安装指南)
- [工具脚本详解](https://www.google.com/search?q=%23工具脚本详解)
  - [`run_circlemap_pipeline.py`](https://www.google.com/search?q=%23run_circlemap_pipelinepy-circle-map-eccdna-鉴定工作流)
  - [`validate_eccdna_amplicons.py`](https://www.google.com/search?q=%23validate_eccdna_amplicons_py-靶向测序eccdna验证流程)
  - [`find_eccdna_terminal_repeats.py`](https://www.google.com/search?q=%23find_eccdna_terminal_repeatspy-eccdna-末端重复序列鉴定流程)
  - [`satCurveSeeker.py`](https://www.google.com/search?q=%23satcurveseekerpy-饱和度曲线数据生成工作流)
  - [`combine_fled_circlemap.py`](https://www.google.com/search?q=%23combine_fled_circlemappy-fled-与-circle-map-结果整合脚本)
  - [`generate_summary_reports.py`](https://www.google.com/search?q=%23generate_summary_reportspy-结果统计与摘要报告生成脚本)
- [如何引用](https://www.google.com/search?q=%23如何引用)
- [许可证](https://www.google.com/search?q=%23许可证)



## 安装指南 (Installation)

#### 1. 克隆本仓库

```bash
git clone https://github.com/your-username/eccToolkit.git
cd eccToolkit
```



#### 2. 安装 Python 依赖库

本工具集依赖于多个 Python 库。我们建议您在虚拟环境中进行安装。

```bash
# 创建并激活虚拟环境 (推荐)
python3 -m venv venv
source venv/bin/activate

# 安装所有必要的 Python 库
pip install -r requirements.txt
```

**注意**: `requirements.txt` 文件需要您手动创建。它应该至少包含以下内容：

Plaintext

```
pandas
pysam
```

您可以通过运行 `pip freeze > requirements.txt` 来生成包含您环境中所有包的完整列表。



#### 3. 安装外部命令行工具



本工具集中的部分脚本需要调用外部的生物信息学软件。请确保以下工具已经安装，并且它们的可执行文件位于您的系统路径 (`PATH`) 中：

- `fastp`
- `bwa`
- `samtools` (>= 1.10)
- `bedtools`
- `blastn` (来自 NCBI BLAST+ 套件)
- `seqkit`
- `Circle-Map`
- `CircleSeeker` (如果相关脚本被使用)
- `FLED` (如果相关脚本被使用)

------



## 工具脚本详解 (Toolkit Scripts)



本章节详细介绍了 `code/` 目录下的每一个核心脚本。



### `run_circlemap_pipeline.py`: Circle-Map eccDNA 鉴定工作流





#### 功能简介



该脚本是一个自动化的生物信息学流程，用于从双端测序数据 (paired-end FASTQ) 中鉴定环状DNA (eccDNA)。它整合了多个生信领域的标准工具，实现从原始数据质控、序列比对到最终使用 Circle-Map 鉴定 eccDNA 的完整分析流程。



#### 依赖项



- **命令行工具**: `fastp`, `bwa`, `samtools`, `Circle-Map`



#### 使用方法



Bash

```
python code/run_circlemap_pipeline.py \
    -t 16 \
    -1 /path/to/your/reads_1.fastq.gz \
    -2 /path/to/your/reads_2.fastq.gz \
    -r /path/to/your/reference.fasta \
    -o /path/to/your/output_directory \
    -s MySample_01
```



#### 参数详解



| 参数                  | 类型  | 是否必须 | 说明                                       |
| --------------------- | ----- | -------- | ------------------------------------------ |
| `-t`, `--threads`     | `int` | 是       | 用于所有分析步骤的线程数。                 |
| `-1`, `--fastq1`      | `str` | 是       | 输入的 FASTQ 文件路径 (Read 1)。           |
| `-2`, `--fastq2`      | `str` | 是       | 输入的 FASTQ 文件路径 (Read 2)。           |
| `-r`, `--reference`   | `str` | 是       | 参考基因组的 FASTA 文件路径。              |
| `-o`, `--output_dir`  | `str` | 是       | 用于存储所有输出文件的目录。               |
| `-s`, `--sample_name` | `str` | 是       | 样本的唯一名称，将用作所有输出文件的前缀。 |



### `validate_eccdna_amplicons.py`: 靶向测序eccDNA验证流程





#### 方法概述



该脚本实现了一种新颖的、用于从靶向扩增子测序数据中**验证**候选eccDNA的方法。其核心思想在于为每个候选eccDNA构建一个“连接点中心”的人工参考序列，然后通过分析比对到该序列上的普通跨越reads和分割跨越reads，来为eccDNA的存在提供高置信度证据。



#### 依赖项



- **命令行工具**: `bwa`, `samtools`
- **Python 库**: `pysam`, `pandas`
- **输入文件格式**: 一个包含 `Chr`, `Start`, `End` 列的TSV文件。



#### 使用方法



Bash

```
python code/validate_eccdna_amplicons.py \
    -i candidates.tsv \
    -1 your_amplicon_reads_R1.fastq.gz \
    -2 your_amplicon_reads_R2.fastq.gz \
    -r /path/to/your/genome.fasta \
    -o MyProject_validation_results \
    -t 16
```

*(详细参数和输出说明请参考该脚本的独立文档部分)*



### `find_eccdna_terminal_repeats.py`: eccDNA 末端重复序列鉴定流程





#### 方法概述



该脚本实现了一种高灵敏度的、用于鉴定候选eccDNA断点处是否存在末端重复序列（正向或反向）的方法。它通过提取每个eccDNA断点两侧的旁翼序列，使用`blastn`进行自我比对，并根据比对结果与断点的邻近度来筛选真实的末端重复。该过程使用多进程并行处理以提高效率。



#### 依赖项



- **命令行工具**: `bedtools`, `blastn`
- **Python 库**: `pysam`, `pandas`
- **输入文件格式**: 一个TSV文件，**关键**：坐标系统需为 **0-based, inclusive**。



#### 使用方法



Bash

```
python code/find_eccdna_terminal_repeats.py \
    -i your_candidates.tsv \
    -g /path/to/your/genome.fasta \
    -o MyProject_repeat_results.csv \
    -t 16 \
    --breakpoint_tolerance 5
```

*(详细参数和输出说明请参考该脚本的独立文档部分)*



### `satCurveSeeker.py`: 饱和度曲线数据生成工作流





#### 功能简介



该脚本通过对输入数据进行系统性的梯度下采样（例如10%, 20%, ..., 90%），并对每个数据子集运行核心分析工具 `CircleSeeker`，来自动化地生成用于绘制饱和度曲线所需的数据。



#### 依赖项



- **命令行工具**: `seqkit`, `CircleSeeker`



#### 使用方法



Bash

```
python code/satCurveSeeker.py \
    -i your_assembled_reads.fasta \
    -r /path/to/your/reference.fasta \
    -p MySample \
    -o /path/to/saturation_curve_analysis \
    -t 16
```

*(详细参数和输出说明请参考该脚本的独立文档部分)*



### `combine_fled_circlemap.py`: FLED 与 Circle-Map 结果整合脚本





#### 功能简介



该脚本用于处理、筛选并整合来自两种不同 eccDNA 鉴定工具（FLED 和 Circle-Map）的输出文件，最终生成筛选后的独立结果以及合并后的总览表格。



#### 依赖项



- **Python 库**: `pandas`



#### 使用方法



Bash

```
python code/combine_fled_circlemap.py \
    --fled-dir /path/to/your/fled_outputs \
    --circlemap-dir /path/to/your/circlemap_outputs \
    --output-dir /path/to/your/analysis_results
```

*(详细参数和输出说明请参考该脚本的独立文档部分)*



### `generate_summary_reports.py`: 结果统计与摘要报告生成脚本





#### 功能简介



该脚本用于对合并后的结果文件进行统计分析，主要计算样本计数和长度分布两个核心指标，并生成易于解读的CSV报告。



#### 依赖项



- **Python 库**: `pandas`



#### 使用方法



Bash

```
python code/generate_summary_reports.py \
    -i /path/to/your/FLED_and_CircleMap_combined.csv \
    -o MyProject_Analysis \
    --sep ","
```

*(详细参数和输出说明请参考该脚本的独立文档部分)*

------



## 如何引用 (Citation)



如果您在您的研究中使用了 **eccToolkit**，请引用我们的论文：

> [请在这里替换成您的论文引用信息]
>
> Zhang, Y., et al. (2025). *Title of Your Amazing Paper*. Journal of Scientific Excellence, 1(1), 1-10.



## 许可证 (License)



本项目采用 [MIT 许可证](https://opensource.org/licenses/MIT)。
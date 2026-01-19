# combine_FLED_CircleMap.py

# eccDNA 数据整合工具使用说明

## 功能概述

`combine_FLED_CircleMap.py` 是一个用于整合和分析环状染色体外DNA（eccDNA）检测数据的Python脚本，专门用于处理和合并来自两种不同检测工具（FLED和Circle-Map）的结果。该脚本能够读取多个样本的数据，应用质量过滤，并生成统一格式的综合报告。

## 主要功能

- 读取和处理FLED工具的OnesegJunction.out输出文件
- 读取和处理Circle-Map工具的BED格式输出文件
- 对两种数据源应用特定的质量过滤标准
- 计算每个eccDNA区域的长度并进行标准化处理
- 整合两种工具检测到的eccDNA位点信息
- 生成CSV和TSV格式的完整和简化结果表格

## 数据处理流程

1. **FLED数据处理**:
   - 从文件名提取样本信息
   - 清理和格式化eccDNA读段ID
   - 创建唯一名称（染色体-起始-结束）
   - 计算eccDNA长度
   - 过滤掉没有完整通过的记录（Nfullpass = 0）
   - 对浮点数值进行四舍五入处理
2. **Circle-Map数据处理**:
   - 应用严格的质量过滤标准：Circle_score > 50, Split_reads > 2, Discordants > 2, Coverage_increase_start > 0.33, Coverage_increase_end > 0.33, 长度 < 10Mb
   - 标准化浮点数和命名格式
   - 添加样本标识和长度计算
3. **结果整合**:
   - 创建一个包含两种工具关键信息的简化合并表格
   - 保持数据格式一致性
   - 按染色体和位置排序

## 输入文件

脚本中预设了以下文件作为输入：

- **FLED文件**:

  ```
  3SEP_LR_25_1.DiGraph.OnesegJunction.out
  3SEP_LR_25_2.DiGraph.OnesegJunction.out
  3SEP_LR_25_3.DiGraph.OnesegJunction.out
  Circel-Seq_LR_25_1.DiGraph.OnesegJunction.out
  Circel-Seq_LR_25_2.DiGraph.OnesegJunction.out
  Circel-Seq_LR_25_3.DiGraph.OnesegJunction.out
  ```

- **Circle-Map文件**:

  ```
  eccDNA_Circel-Seq_NGS_25_1_CM.bed
  eccDNA_Circel-Seq_NGS_25_2_CM.bed
  eccDNA_Circel-Seq_NGS_25_3_CM.bed
  ```

## 输出文件

脚本生成以下输出文件：

- **FLED_Merged_Results.csv/tsv**: 合并所有FLED样本的详细结果
- **Circle_Map_Merged_Results.csv/tsv**: 合并所有Circle-Map样本的详细结果
- **Combined_Basic_Results.csv/tsv**: 包含两种工具基本信息的综合表格

## 使用方法

```bash
python combine_FLED_CircleMap.py
```

运行脚本时，确保当前工作目录中存在所有需要处理的输入文件。脚本会自动处理文件并生成结果。

## 处理逻辑和过滤标准

1. **FLED数据过滤**:
   - 仅保留Nfullpass > 0的记录
2. **Circle-Map数据过滤**:
   - Circle_score > 50
   - Split_reads > 2
   - Discordants > 2
   - Coverage_increase_start > 0.33
   - Coverage_increase_end > 0.33
   - eccDNA长度 < 10,000,000 bp

通过这些严格的过滤标准，该脚本能够提高eccDNA检测结果的可靠性，并为后续分析提供高质量的数据集。

## 优势特点

- 自动处理多个样本数据
- 应用严格的质量控制标准
- 生成标准化且易于分析的输出格式
- 提供详细的处理日志和统计信息
- 支持不同工具结果的整合比较
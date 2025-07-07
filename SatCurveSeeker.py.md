# SatCurveSeeker.py

# eccDNA 饱和度曲线分析工具使用说明

## 功能概述

`satCurveSeeker.py` 是一个用于进行环状染色体外DNA（eccDNA）测序饱和度分析的Python脚本。该工具通过对输入的FASTA文件进行不同比例的下采样，并对每个采样文件运行CircleSeeker分析，以评估测序深度与eccDNA检测效果之间的关系。这种分析有助于确定最优测序深度，避免数据不足或过度测序。

## 主要功能

1. **序列下采样**：使用seqkit工具对输入的FASTA文件进行多个比例的下采样（10%至90%，每10%一档）
2. **eccDNA检测**：对每个采样文件执行CircleSeeker分析
3. **串行处理**：采用严格串行执行模式，确保资源合理利用
4. **完整日志**：记录每个步骤的执行时间和结果

## 工作流程

1. 脚本接收输入FASTA文件和参考基因组等参数
2. 创建输出目录（如果不存在）
3. 按10%、20%...90%的比例，依次对输入文件进行下采样
4. 对每个采样文件运行CircleSeeker（等待一个完成后再执行下一个）
5. 记录详细的执行日志和时间消耗

## 参数说明

- `-i, --input_fasta`：输入的FASTA格式文件，通常是长读测序数据
- `-r, --reference`：参考基因组FASTA文件
- `-p, --prefix`：输出文件的样本前缀
- `-o, --outdir`：结果输出目录
- `-t, --threads`：分析过程使用的线程数

## 使用示例

```bash
python /data9/home/yxzhang/Tools/eccToolkit/satCurveSeeker.py \
  -i GBM_JBH.fasta \
  -r /data9/home/yxzhang/Genome_Ref/Human/hg38.p13/GRCh38.p13.fna \
  -p GBM_P1 \
  -o GBM_P1_saturation \
  -t 32 &
```

这个命令在后台运行饱和度曲线分析，使用GBM_JBH.fasta作为输入文件，hg38参考基因组，将结果保存在GBM_P1_saturation目录，使用32个线程进行分析，并为所有输出文件添加GBM_P1前缀。

## 输出文件

脚本会在指定的输出目录中生成多个文件：

1. 下采样后的FASTA文件：如`GBM_P1_10Per.fasta`、`GBM_P1_20Per.fasta`等
2. 每个采样比例对应的CircleSeeker结果文件
3. 详细的执行日志信息

## 注意事项

- 确保已安装seqkit和CircleSeeker工具，并且它们在系统PATH中可用
- 分析过程为严格串行执行，完整分析可能需要较长时间
- 为确保结果可重复性，下采样使用固定的随机种子（42）
- 该工具支持X染色体分析（通过`--enable_X`参数传递给CircleSeeker）

通过分析不同测序深度下的eccDNA检测效果，研究人员可以确定最佳的测序深度，优化实验设计和资源利用。
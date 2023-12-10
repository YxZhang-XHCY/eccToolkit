#!/bin/bash

# 默认覆盖度阈值
coverage_threshold=0.85

# 默认目录为当前目录
directory="."

# 使用 getopts 解析命令行参数
while getopts "t:d:" opt; do
  case $opt in
    t) coverage_threshold=$OPTARG ;;
    d) directory=$OPTARG ;;
    \?) echo "无效选项: -$OPTARG" >&2; exit 1 ;;
  esac
done

# 遍历指定目录下所有以 _CM.bed 结尾的文件
for bed_file in $directory/eccDNA_*_CM.bed; do
    # 检查文件是否存在
    if [ -f "$bed_file" ]; then
        # 提取文件名中的基本部分，去掉 eccDNA_ 和 _CM.bed 后缀
        base_name=$(basename $bed_file | sed 's/eccDNA_//' | sed 's/_CM.bed$//')

        # 应用 awk 命令进行过滤，并输出到新文件
        filt_file="${directory}/${base_name}.filt.csv"
        awk '{if (($4 >= 1) && ($5 >= 2) && ($6 >=50) && ($9 >=0.33) && ($10 >= 0.33) ){print $0; }}' "$bed_file" > "$filt_file"

        # 定义与 filt 文件对应的 ef 文件
        ef_file="${directory}/${base_name}.csv"

        # 检查对应的 ef 文件是否存在
        if [ -f "$ef_file" ]; then
            # 构建输出文件名
            output_file="${directory}/${base_name}.merge.${coverage_threshold}.csv"

            # 执行 bedtools intersect 命令
            nohup bedtools intersect -nonamecheck -a "$filt_file" -b "$ef_file" -f $coverage_threshold -wa -wb > "$output_file" &

            # 等待 bedtools intersect 命令完成
            wait

            # 处理 merge 文件
            final_output="${directory}/${base_name}.final_output.${coverage_threshold}.csv"
            # 提取前11列并去除重复行
            awk '!seen[$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11]++' "$output_file" | cut -f1-11 > temp_file.csv

            # 生成新的 eccDNA_Name 列和 eccDNA_length 列，并重新排列列
            awk -F'\t' 'BEGIN{OFS="\t"} {print "'$base_name'", $1,$2,$3,$1"_"$2"_"$3,$3-$2+1,$4,$5,$6,$7,$8,$9,$10,$11}' temp_file.csv > temp_file_12col.csv

            # 添加行名
            {
                echo -e "Sample\tChr_C\tStart_C\tEnd_C\teccDNA_Name\teccDNA_length\tDiscordants\tSplit_Reads\tCircle_Score\tMean_Cov\tStd_Dev\tStart_Cov_Inc\tEnd_Cov_Inc\tCov_Cont"
                cat temp_file_12col.csv
            } > "$final_output"

            # 清理临时文件
            rm temp_file.csv temp_file_12col.csv

        else
            echo "未找到与 $filt_file 对应的 ef 文件: $ef_file"
        fi
    else
        echo "未找到文件: $bed_file"
    fi
done

#!/bin/bash

# Pipeline 1: Software Installation
# ----------------------------------
echo "Starting software installation..."
mamba install -c bioconda seqkit=2.3.1
mamba install -c agbiome bbtools
echo "Software installation completed."

# Pipeline 2: Generate Subsets based on Declining Coverage
# --------------------------------------------------------
echo "Generating subsets based on declining coverage..."
input_r1="/data5/home/yxzhang2/EccDNA/simulation/GE-2_combined.R1.fastq.gz"
input_r2="/data5/home/yxzhang2/EccDNA/simulation/GE-2_combined.R2.fastq.gz"
initial_depth=40

for i in {1..8}; do
  subset_coverage=$(($initial_depth - 4 * $i))
  retain_probability=$(echo "scale=2; $subset_coverage / $initial_depth" | bc)
  output_r1="GE_fastq_R1_${subset_coverage}X.fastq.gz"
  output_r2="GE_fastq_R2_${subset_coverage}X.fastq.gz"

  echo "Generating subset with coverage ${subset_coverage}X"
  seqtk sample -s100 $input_r1 $retain_probability | gzip > $output_r1 
  seqtk sample -s100 $input_r2 $retain_probability | gzip > $output_r2
done
echo "Subsetting based on declining coverage completed."

# Pipeline 3: Parallel Sampling of Data
# -------------------------------------
echo "Starting parallel sampling of data..."
declare -A samples=(
  ["/data5/home/yxzhang2/All_data/circle-seq/GE-1_2nd.circle-seq.R1.fastq.gz"]=0.51
  ["/data5/home/yxzhang2/All_data/circle-seq/GE-1_2nd.circle-seq.R2.fastq.gz"]=0.51
  # ... Add all other samples and their associated probabilities here ...
  ["/data5/home/yxzhang2/All_data/circle-seq/stem-2.circle-seq.R1.fastq.gz"]=0.63
  ["/data5/home/yxzhang2/All_data/circle-seq/stem-2.circle-seq.R2.fastq.gz"]=0.63
)

for sample in "${!samples[@]}"; do
  base_name=$(basename $sample)
  output_file="${base_name%.*}_12X.fastq.gz"
  nohup seqtk sample -s100 $sample ${samples[$sample]} > $output_file &
done

wait
echo "Parallel sampling completed."

# End of the script
echo "All pipelines executed successfully!"


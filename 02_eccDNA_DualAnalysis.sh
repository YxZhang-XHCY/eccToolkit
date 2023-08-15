#!/bin/bash
set -e
trap 'echo Pipeline failed at: $BASH_COMMAND' ERR

while read line; do

  arr=($line)

  ref=${arr[0]}
  intdir=${arr[1]}
  outdir=${arr[2]}
  oldnam=${arr[3]}
  newnam=${arr[4]}
  Tissue=${arr[5]}
  Thread=${arr[6]}

  # Ecc finder
  cd $intdir
  mkdir output1
  cd output1
  mkdir $newnam
  
  cd /data5/home/yzhang/test_zyx/ecc_finder
  ecc_finder_run=0
  nohup python ecc_finder.py map-sr $ref \
    $intdir/$oldnam.R1.fastq.gz \
    $intdir/$oldnam.R2.fastq.gz \
    -t $Thread \
    -r $ref \
    -o $intdir/output1/$newnam \
    -x $newnam > nohup.$newnam.ef.out || ecc_finder_run=$?
  
  if [ $ecc_finder_run -ne 0 ]; then
    echo "Ecc finder failed for sample $newnam" >&2
    exit 1
  fi

  cd $intdir/output1/$newnam
  cp $newnam.csv $outdir/eccDNA_$newnam.ef.csv

  # Circle map 
  cd $intdir/output1/$newnam/align_files
  rm *fastq.gz $newnam.bam

  nohup samtools sort -@ $Thread -n -o qname_eccDNA_$newnam.bam $newnam.sam
  nohup samtools sort -@ $Thread -o sorted_eccDNA_$newnam.bam $newnam.sam
  nohup samtools index -@ $Thread sorted_eccDNA_$newnam.bam

  nohup qualimap bamqc -bam sorted_eccDNA_$newnam.bam -outfile $newnam.pdf \
    -outformat PDF --java-mem-size=100G

  cd $outdir
  mkdir qualimap
  cd $intdir/output1/$newnam/align_files/sorted_eccDNA_$newnam\_stats
  mv *pdf $outdir/qualimap

  cd $intdir/output1/$newnam/align_files
  rm $newnam.sam

  nohup /data6/tool/anaconda3-2019.07/bin/Circle-Map ReadExtractor \
    -i qname_eccDNA_$newnam.bam -o eccDNA_$newnam\_candidates.bam

  nohup samtools sort -@ $Thread eccDNA_$newnam\_candidates.bam -o sort_eccDNA_$newnam\_candidates.bam

  rm eccDNA_$newnam\_candidates.bam

  nohup samtools index -@ $Thread sort_eccDNA_$newnam\_candidates.bam

  nohup /data6/tool/anaconda3-2019.07/bin/Circle-Map Realign -t $Thread \
    -i sort_eccDNA_$newnam\_candidates.bam -qbam qname_eccDNA_$newnam.bam \
    -sbam sorted_eccDNA_$newnam.bam \
    -fasta $ref -o eccDNA_$newnam\_CM.bed > nohup_$newnam\_CM

  awk '{if (($4 >= 1) && ($5 >= 2) && ($6 >=50) && ($9 >=0.33) && ($10 >= 0.33) ){print $0; }}' \
    eccDNA_$newnam\_CM.bed > eccDNA_$newnam\_CM.filt.csv

  cp eccDNA_$newnam\_CM.bed eccDNA_$newnam\_CM.filt.csv $outdir

  cd $outdir
  mv *.bam *.bai $outdir/bam

  cd $outdir
  nohup bedtools intersect -a eccDNA_$newnam\_CM.filt.csv -b eccDNA_$newnam.ef.csv -f 0.50 \
    -wa -wb > $newnam.EccDNA.merge.csv

  Rscript /data5/home/yzhang/test_zyx/Merge.R $newnam $Tissue
  
done < $1

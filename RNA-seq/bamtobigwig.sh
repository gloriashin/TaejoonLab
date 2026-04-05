#!/bin/bash


for bam in $(ls *.bam)
do
  OUT=${bam/.bam/}
  
  OUT_bigwig=$OUT"_bam.bigwig"
  


  echo "##########"
  echo $bam
  bamCoverage -b $bam -o $OUT_bigwig

done
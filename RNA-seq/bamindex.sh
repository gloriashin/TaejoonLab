#!/bin/bash


for bam in $(ls *.bam)
do
  OUT=${bam/.bam/}
  
  
  


  echo "##########"
  echo $bam
  samtools index $bam

done
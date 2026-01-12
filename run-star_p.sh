#!/bin/bash
NUM_THREADS=6

DB="/data/home/ssy0213/plouhinec2017/DB/db.STAR/xenLae10"
GTF="/data/home/ssy0213/plouhinec2017/DB/xenLae10.2025_02.gtf"

for FQ1 in $(ls /data/project.xenopus/xenopus.tx.XENLA/Plouhinec2017_EctoMap.PRJNA400602_XENLAtx/fastq/*NBa*_R1.trim.fastq.gz)
do
  FQ2=${FQ1/_R1/_R2}
  SAM=${FQ1/_R1.trim.fastq.gz}".xenLae10.STAR.sam"
  SAM=$(basename $SAM)
  OUT="Start aligning "${FQ1/_R1}
   

  BAM=${SAM/.sam/}"unsorted.bam"


  if [ ! -e $BAM ]; then
    echo $OUT
    STAR --genomeDir $DB --runThreadN $NUM_THREADS --readFilesCommand gunzip -c --readFilesIn $FQ1 $FQ2 --sjdbGTFfile $GTF --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix $BAM

    
  fi
done

for FQ1 in $(ls /data/project.xenopus/xenopus.tx.XENLA/Plouhinec2017_EctoMap.PRJNA400602_XENLAtx/fastq/*NBp*_R1.trim.fastq.gz)
do
  FQ2=${FQ1/_R1/_R2}
  SAM=${FQ1/_R1.trim.fastq.gz}".xenLae10.STAR.sam"
  SAM=$(basename $SAM)
  OUT="Start aligning "${FQ1/_R1}
  

  BAM=${SAM/.sam/}".bam"


  if [ ! -e $BAM ]; then
    echo $OUT
    STAR --genomeDir $DB --runThreadN $NUM_THREADS --readFilesCommand gunzip -c --readFilesIn $FQ1 $FQ2 --sjdbGTFfile $GTF --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix $BAM

  fi
done

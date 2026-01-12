#!/bin/bash
NUM_THREADS=6

DB="$HOME/project.xenopus/db/db.hisat2/xenLae10"

for FQ1 in $(ls ../fastq/*_R1.trim.fastq.gz)
do
  FQ2=${FQ1/_R1/_R2}
  SAM=${FQ1/_R1.trim.fastq.gz}".xenLae10.hisat2.sam"
  SAM=$(basename $SAM)
  echo $OUT

  BAM=${SAM/.sam/}".bam"
  SORT_N=${SAM/.sam/}".sortN.bam"
  FIXMATE=${SAM/.sam/}".fixmate.bam"
  SORT_C=${SAM/.sam/}".sortC.bam"
  GTF=${SAM/.sam/}".stringtie.gtf"

  if [ ! -e $GTF ]; then
    hisat2 -p $NUM_THREADS -x $DB -1 $FQ1 -2 $FQ2 -S $SAM
    samtools view -@ $NUM_THREADS -bS -o $BAM $SAM
    rm $SAM

    samtools sort -@ $NUM_THREADS -n -o $SORT_N $BAM
    samtools fixmate -@ $NUM_THREADS $SORT_N $FIXMATE
    samtools sort -@ $NUM_THREADS -o $SORT_C $FIXMATE
    rm $FIXMATE $SORT_N

    stringtie $SORT_C -p $NUM_THREADS -o $GTF 
  fi
done

REF_GTF="$HOME/project.xenopus/xenopus.tx.XENLA/xenLae10.2025_01.gtf"

for SORT_C in $(ls *.sortC.bam)
do
  GTF_OUT=${SORT_C/.sortC.bam/}".stringtie_xb2025-01.gtf"
  TPM_OUT=${SORT_C/.sortC.bam/}".stringtie_xb2025-01.tpm.tsv"

  stringtie $SORT_C -G $REF_GTF -A $TPM_OUT -p $NUM_THREADS -o $GTF_OUT --ref $REF_FA
done

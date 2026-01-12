#!/bin/bash
NUM_THREADS=4

for IN_1 in $(ls *_R1.raw.fastq.gz)
do
  IN_2=${IN_1/_R1/_R2}
  OUT=${IN_1/_R1.raw.fastq/}
  OUT=${OUT/.gz/}
  OUT_HTML=$OUT".fastp_log.html"
  OUT_JSON=$OUT".fastp_log.json"

  OUT_1=$OUT"_R1.trim.fastq.gz"
  OUT_2=$OUT"_R2.trim.fastq.gz"

  echo "##########"
  echo $IN_1 $IN_2
  fastp --thread $NUM_THREADS --in1 $IN_1 --in2 $IN_2 --out1 $OUT_1 --out2 $OUT_2 \
    -l 40 -h $OUT_HTML -j $OUT_JSON 
    # --unpaired1 $UN_1 --unpaired2 $UN_2 
    # --reads_to_trim 1000000
done

# PacBio onso reference
# cutadapt --cores=$NUM_THREADS --max-n 1 --pair-filter=any \
#   -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
#  --minimum-length 50 \
#  -o $FQ1_OUT -p $FQ2_OUT \
#  $FQ1 $FQ2

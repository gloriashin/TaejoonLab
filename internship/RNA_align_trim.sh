#!/usr/bin/env bash

TRIMMOMATIC="$HOME/miniconda3/bin/trimmomatic"
BWA="$HOME/miniconda3/bin/bwa"
STAR="$HOME/miniconda3/bin/STAR"

THREAD_NUM=6

BWA_DB="/work/SeqRef/bwaRef/RAT/rat98_longest_bwa"
BWA_DB_NAME=$(basename $BWA_DB)
STAR_DB_DIR="/home/ssy0213/flyRNAseq/dmel-all-chromosome-r6.46"
STAR_DB_NAME="droso_6.46"

for FASTQ in $(ls *_R1.fastq.gz)
do
		FASTQ_NM=$(basename $FASTQ _R1.fastq.gz)
		FASTQ2=$FASTQ_NM"_R2.fastq.gz"
		TRIM_FASTQ1=$FASTQ_NM".R1.trimmed.fastq.gz"
		unTRIM_FASTQ1=$FASTQ_NM".R1un.trimmed.fastq.gz"
		TRIM_FASTQ2=$FASTQ_NM".R2.trimmed.fastq.gz"
		unTRIM_FASTQ2=$FASTQ_NM".R2un.trimmed.fastq.gz"
		TRIM_LOG=$FASTQ_NM".trimmed.log"
		if [ ! -e $TRIM_FASTQ1 ]; then
			
			$TRIMMOMATIC PE -threads $THREAD_NUM -trimlog $TRIM_LOG $FASTQ $FASTQ2 $TRIM_FASTQ1 $unTRIM_FASTQ1 $TRIM_FASTQ2 $unTRIM_FASTQ2 ILLUMINACLIP:NEB_adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

		fi
	
		SAM_NAME1=$(basename $TRIM_FASTQ1 .fastq.gz)
		SAM_NAME2=$(basename $TRIM_FASTQ2 .fastq.gz)
        BWA_SAM=$SAM_NAME"."$BWA_DB_NAME".bwa_mem.sam"

#        if [ -e $TRIM_FASTQ ] && [! -e $BWA_SAM ]; then
#
#                echo "Make $SAM by BWA" 
#                $BWA mem -t $THREAD_NUM $BWA_DB $TRIM_FASTQ > $BWA_SAM 2>log
#
#        fi
	
		STAR_SAM=$FASTQ_NM"."$STAR_DB_NAME".STAR."
		STAR_SAM_DIR="$FASTQ_NM"
		if [ -e $TRIM_FASTQ ] && [ ! -e $STAR_SAM ]; then
		
				echo "Make $SAM by STAR"
                if [ ! -d $STAR_SAM_DIR ]; then
	        
						mkdir $STAR_SAM_DIR

	        	fi		
				$STAR --genomeDir $STAR_DB_DIR --runThreadN $THREAD_NUM --readFilesIn $TRIM_FASTQ1 $TRIM_FASTQ2 --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts \
						--outFileNamePrefix ./$STAR_SAM_DIR/$STAR_SAM  
		
		
		fi
done





#!/bin/bash

#https://docs.nvidia.com/clara/parabricks/4.0.1/documentation/tooldocs/compatiblecpusoftwareversions.html

REF=/scratch/user/u.gz28467/input_files/TAIR10_chr_all.fasta

TDR_1=/scratch/user/u.gz28467/input_files/TDr-7_10M_1.fastq
TDR_2=/scratch/user/u.gz28467/input_files/TDr-7_10M_2.fastq

# Number of threads for GATK to use
OMP_NUM_THREADS=8

function fq2bam_pe {
        # https://docs.nvidia.com/clara/parabricks/4.0.1/documentation/tooldocs/man_fq2bam.html#man-fq2bam
        REF=$1; FQ1=$2; FQ2=$3; PREFIX=$4
        # Align and sort reads
        bwa mem -t 32 -K 10000000 \
                -R '@RG\tID:sample_rg1\tLB:lib1\tPL:bar\tSM:sample\tPU:sample_rg1' \
                $REF $FQ1 $FQ2 | gatk SortSam \
                        --java-options -Xmx30g --MAX_RECORDS_IN_RAM 5000000 \
                        -I /dev/stdin -O ${PREFIX}_sorted.bam --SORT_ORDER coordinate
        
        # Mark duplicates
        gatk MarkDuplicates --java-options -Xmx30g \
                -I ${PREFIX}_sorted.bam \
                -O ${PREFIX}_sorted_marked.bam \
                -M ${PREFIX}_metrics.txt
        rm ${PREFIX}_sorted.bam
        # Index the bam file
        samtools index ${PREFIX}_sorted.bam
}

# Align Bay-0 reads
time fq2bam_pe $REF $TDR_1 $TDR_2 TDr-7_10M_cpu

# Each took about 8 minutes to run on 16 CPUs
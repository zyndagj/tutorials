#!/bin/bash

#https://docs.nvidia.com/clara/parabricks/4.0.1/documentation/tooldocs/compatiblecpusoftwareversions.html

REF=/scratch/user/u.gz28467/input_files/TAIR10_chr_all.fasta

BAY_FQ=/scratch/user/u.gz28467/input_files/Bay-0_20M.fastq
SHA_FQ=/scratch/user/u.gz28467/input_files/Shadara_20M.fastq

function fq2bam_se {
        # https://docs.nvidia.com/clara/parabricks/4.0.1/documentation/tooldocs/man_fq2bam.html#man-fq2bam
        REF=$1; FQ=$2; PREFIX=$3
        # Align and sort reads
        bwa mem -t 32 -K 10000000 \
                -R '@RG\tID:sample_rg1\tLB:lib1\tPL:bar\tSM:sample\tPU:sample_rg1' \
                $REF $FQ | ${GATK} SortSam \
                        --java-options -Xmx30g --MAX_RECORDS_IN_RAM 5000000 \
                        -I /dev/stdin -O ${PREFIX}_sorted.bam --SORT_ORDER coordinate
        
        # Mark duplicates
        gatk MarkDuplicates --java-options -Xmx30g \
                -I ${PREFIX}_sorted.bam \
                -O ${PREFIX}_sorted_marked.bam \
                -M ${PREFIX}_metrics.txt
        rm ${PREFIX}_sorted.bam
}

# Align Bay-0 reads
time fq2bam_se $REF $BAY_FQ Bay-0_20M_cpu

# Align Shadara reads
time fq2bam_se $REF $SHA_FQ Shadara_20M_cpu

# Each took about 8 minutes to run on 16 CPUs
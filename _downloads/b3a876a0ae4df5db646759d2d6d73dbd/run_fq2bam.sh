#!/bin/bash

REF=/scratch/user/u.gz28467/input_files/TAIR10_chr_all.fasta

TDR_1=/scratch/user/u.gz28467/input_files/TDr-7_10M_1.fastq
TDR_2=/scratch/user/u.gz28467/input_files/TDr-7_10M_2.fastq

# Align TDr-7 reads
time pbrun fq2bam --ref $REF --in-fq $TDR_1 $TDR_2 --out-bam TDr-7_10M_pb.bam --num-gpus 1

#!/bin/bash

REF=input_files/TAIR10_chr_all.fasta

TDR_1=input_files/TDr-7_10M_R1.fastq
TDR_2=input_files/TDr-7_10M_R2.fastq

# Align TDr-7 reads
time pbrun fq2bam --ref $REF --in-fq $TDR_1 $TDR_2 --out-bam TDr-7_10M_pb.bam --num-gpus 1
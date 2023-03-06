#!/bin/bash

REF=/scratch/user/u.gz28467/input_files/TAIR10_chr_all.fasta

BAY_FQ=/scratch/user/u.gz28467/input_files/Bay-0_20M.fastq
SHA_FQ=/scratch/user/u.gz28467/input_files/Shadara_20M.fastq

# Align Bay-0 reads
time pbrun fq2bam --ref $REF --in-se-fq $BAY_FQ --out-bam Bay-0_20M_pb.bam --num-gpus 1

# Align Shadara reads
time pbrun fq2bam --ref $REF --in-se-fq $SHA_FQ --out-bam Shadara_20M_pb.bam --num-gpus 1
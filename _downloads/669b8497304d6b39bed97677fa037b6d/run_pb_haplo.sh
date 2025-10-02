#!/bin/bash

REF=input_files/TAIR10_chr_all.fasta

TDR_BAM=TDr-7_10M_pb.bam

# Call variants on TDr-7 bam file
time pbrun haplotypecaller --ref $REF --in-bam $TDR_BAM --out-variants TDr-7_10M_pb.vcf
#!/bin/bash

REF=/scratch/user/u.gz28467/input_files/TAIR10_chr_all.fasta

BAY_BAM=Bay-0_20M_pb.bam
SHA_BAM=Shadara_20M_pb.bam

# Call variants on Bay-0 bam file
time pbrun haplotypecaller --ref $REF --in-bam $BAY_BAM --out-variants Bay-0_20M_pb.vcf

# Call variants on Shadara bam file
time pbrun haplotypecaller --ref $REF --in-bam $SHA_BAM --out-variants Shadara_20M_pb.vcf
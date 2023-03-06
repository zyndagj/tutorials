#!/bin/bash

REF=/scratch/user/u.gz28467/input_files/TAIR10_chr_all.fasta

BAY_BAM=Bay-0_20M_cpu_sorted_marked.bam
SHA_BAM=Shadara_20M_cpu_sorted_marked.bam

function haplo_se {
	REF=$1; BAM=$2; PREFIX=$3
	gatk HaplotypeCaller \
		--java-options -Xmx30g \
		--input $BAM \
		--output ${PREFIX}_cpu.vcf \
		--reference ${REF} \
		--native-pair-hmm-threads 16
}

# Index the inputs
samtools index $BAY_BAM
samtools index $SHA_BAM

# Call variants on Bay-0 bam file (11min)
time haplo_se $REF $BAY_BAM Bay-0_20M

# Call variants on Shadra bam file (13min)
time haplo_se $REF $SHA_BAM Shadara_20M
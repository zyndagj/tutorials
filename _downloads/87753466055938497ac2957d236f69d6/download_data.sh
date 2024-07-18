#!/bin/bash

[ -e input_data ] || mkdir input_data
cd input_data

# Download A. thaliana SNPs
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
tabix -p vcf 1001genomes_snp-short-indel_only_ACGTN.vcf.gz

# Download and downsample TDr-7 data
# https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP012869
curl -o - "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR519/SRR519591/SRR519591_1.fastq.gz" | zcat | seqtk sample - 10000000 > TDr-7_10M_R1.fastq
curl -o - "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR519/SRR519591/SRR519591_2.fastq.gz" | zcat | seqtk sample - 10000000 > TDr-7_10M_R2.fastq

# Download the TAIR10 reference and build indices
curl https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz | zcat > TAIR10_chr_all.fasta
samtools faidx TAIR10_chr_all.fasta
bwa index TAIR10_chr_all.fasta
gatk CreateSequenceDictionary -R TAIR10_chr_all.fasta

# SNPmatch database for 1001 Genomes
wget https://figshare.com/ndownloader/files/9547051 && tar -xf 9547051 && rm 9547051
#!/bin/bash

# Download A. thaliana SNPs
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
tabix -p vcf 1001genomes_snp-short-indel_only_ACGTN.vcf.gz

# Download and downsample Shadara
wget https://1001genomes.org/data/JGI/JGIHeazlewood2008/releases/2011_11_15/TAIR10/strains/Sha/Shadara.snp.aa.txt.gz
curl https://1001genomes.org/data/JGI/JGIHeazlewood2008/releases/2011_11_15/TAIR10/strains/Sha/Shadara.bam | samtools fastq - | singularity exec seqtk_v1.3-1-deb_cv1.sif seqtk sample - 20000000 > Shadara_20M.fastq

# Download and downsample Bay-0
wget https://1001genomes.org/data/JGI/JGIHeazlewood2008/releases/2011_11_15/TAIR10/strains/Bay-0/Bay-0.snp.aa.txt.gz
curl https://1001genomes.org/data/JGI/JGIHeazlewood2008/releases/2011_11_15/TAIR10/strains/Bay-0/Bay-0.bam | samtools fastq - | singularity exec seqtk_v1.3-1-deb_cv1.sif seqtk sample - 20000000 > Bay-0_20M.fastq

# Download the TAIR10 reference and build indices
curl https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz | zcat - > TAIR10_chr_all.fasta
singularity exec bwa_v0.7.15_cv4.sif bwa index TAIR10_chr_all.fasta
singularity exec gatk_4.2.0.0.sif gatk CreateSequenceDictionary -R TAIR10_chr_all.fasta
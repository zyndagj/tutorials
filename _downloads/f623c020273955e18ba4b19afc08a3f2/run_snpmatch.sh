#!/bin/bash

VCF=TDr-7_10M_pb.vcf

RM=/scratch/user/u.gz28467/input_files/1001genomes_snp-short-indel_only_ACGTN.BIALLELIC.hdf5
CM=/scratch/user/u.gz28467/input_files/1001genomes_snp-short-indel_only_ACGTN.BIALLELIC.acc.hdf5

snpmatch inbred -i $VCF -d $RM -e $CM -o TDr-7_10M_pb_snpmatch
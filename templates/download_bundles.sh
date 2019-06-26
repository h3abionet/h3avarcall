#!/bin/bash

## Download data
# lftp -e "pget -n 20 ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz; bye"
# lftp -e "pget -n 20 ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz; bye"

## Download all the files in the list
while read file
do
    lftp -e "pget -n 20 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/$file; bye"
done < !{b37_list}

## Unzip all the files
for file in *.gz;
do
    gunzip $file
done

## Generate BWA indexe
bwa index human_g1k_v37_decoy.fasta



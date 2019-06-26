#!/bin/bash

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



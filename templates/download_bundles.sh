#!/bin/bash

while read file
do
    lftp -e "pget -n 20 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/$file; bye"
done < b37_files_minimal.txt

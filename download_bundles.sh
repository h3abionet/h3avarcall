#!/bin/bash

# curl --list-only ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/ > b37_files.txt

while read file
do
    lftp -e "pget -n 20 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/$file; bye"
done < b37_files.txt

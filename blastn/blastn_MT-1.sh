#!/bin/bash

# blastn: 2.10.1+
# Package: blast 2.10.1, build May 12 2020 13:06:02
#
# bin:   ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
# nt.gz: ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/
#
# create a blast database:  makeblastdb -in nt -input_type fasta -dbtype nucl
#
# run example:  blastn -db nt -query test.fasta -out o.results -num_threads 16

files=$(ls *fasta)

for FILE in $(echo $files) ; do

 head -n250  ${FILE} > temp_${FILE} 

 blastn -db nt -query  temp_${FILE}  -out BLASTn_results_${FILE}.txt -num_threads 16 -max_target_seqs 10

 rm temp_${FILE} 

done



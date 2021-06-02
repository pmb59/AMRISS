#blastn: 2.10.1+
# Package: blast 2.10.1, build May 12 2020 13:06:02
 
files=$(ls *fasta)

for FILE in  $(echo $files) ; do

    head -n20  ${FILE} > temp_${FILE} 

    blastn -db nt -query  temp_${FILE}  -out BLASTn_results_${FILE}.txt  -remote  -max_target_seqs 10

    #rm temp_${FILE} 

done



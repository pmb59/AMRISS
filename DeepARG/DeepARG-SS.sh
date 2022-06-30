#!/bin/bash

# convert FASTQ to FASTA
cd /GLDS69
L=$(ls F*.fastq.gz)

cd /seqtk
for i in $L; do
  echo $i
  seqtk seq -a /GLDS69/${i} >/GLDS69/${i}.fa
done

#merge R1 and R2 in one FASTA file
cd /GLDS69
N=$(ls -a F*_R1-clean.fastq.gz | cut -d "R" -f 1)

for j in $N; do
  cat ${j}*.fa >${j}.fa
  cd /deeparg-ss
  python deepARG.py --align --type nucl --reads --input /GLDS69/${j}.fa --output /GLDS69/${j}.out
done

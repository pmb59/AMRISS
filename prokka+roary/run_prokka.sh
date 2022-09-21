#!/bin/bash

pushd MT-1
contigs=$(ls -a *.scaffolds.fasta | cut -d "." -f 1)
echo $contigs && popd

docker pull staphb/prokka:latest

for j in $contigs; do
    docker run -v $PWD:/data staphb/prokka:latest prokka --prefix $j --cpus 2 --outdir /data/prokka-output/prokka-output-$j /data/MT-1/${j}.scaffolds.fasta
    tar -zcvf ./prokka-output/prokka-output-${j}.tar.gz ./prokka-output/prokka-output-$j
done

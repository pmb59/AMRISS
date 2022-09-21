#!/bin/bash

cd /GLIMMER/MT-1
N=$(ls -a *.scaffolds.fasta | cut -d "." -f 1)

cd /GLIMMER/glimmer3.02/bin

for j in $N; do
    echo $j

    genome='/GLIMMER/MT-1/'${j}'.scaffolds.fasta'
    orfs='/GLIMMER/results/out.orfs'${j}
    runtrain='/GLIMMER/results/run1.train'${j}
    icm='/GLIMMER/results/run1.icm'${j}
    predictions='/GLIMMER/results/'${j}
    predictionsFASTA="/GLIMMER/results/ORF_"${j}".fa"

    long-orfs -n $genome $orfs
    extract -t $genome $orfs >$runtrain
    build-icm $icm <$runtrain
    # minimum gene length to 50 bp
    glimmer3 -g50 $genome $icm $predictions
    #get data for DeepARG
    extract -t $genome ${predictions}".predict" >$predictionsFASTA
    rm $orfs
    rm $runtrain
    rm $icm

done

#!/bin/bash

cd /GLIMMER/results/
N=$(ls -a ORF*.fa | cut -d "." -f 1)

for j in $N; do
  cd /deeparg-ss
  python deepARG.py --align --type nucl --genes --input /GLIMMER/results/${j}.fa --output /GLIMMER/deeparg/${j}.out
done

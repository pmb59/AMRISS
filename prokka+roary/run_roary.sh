#!/bin/bash

docker pull sangerpathogens/roary

# gff files are stored in $PWD/data

docker run --rm -it -v $PWD/data:/data sangerpathogens/roary roary -r -p 4 -f /data/b_cereus ./data/*.gff
 
docker run --rm -it -v $PWD/data:/data sangerpathogens/roary roary -r -p 4 -f /data/e_bugandensis ./data/*.gff

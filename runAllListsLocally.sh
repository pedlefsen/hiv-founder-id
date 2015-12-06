#!/bin/bash
for patient in `ls -c1 ${1}/*.fasta | egrep --only "_[0-9]{5,}_" | tr -d "_" | uniq`
do 
    ./runListLocally.bash ${1} $patient &
done

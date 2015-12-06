#!/bin/bash
mkdir ${2}
for patient in  `ls -c1 ${1}/*.fasta | egrep --only "_[0-9]{5,}_" | tr -d "_" | uniq`
do 
   ./runListLocally2.bash ${1} ${2} $patient &
done

#!/bin/bash
mkdir ${3}
for patient in  `ls -c1 ${1}/*.fasta | egrep --only "_[0-9]{5,}_" | tr -d "_" | uniq`
do
   ./evaluateFoundersFromInfer.bash ${1} ${2} ${3} ${patient} &
done

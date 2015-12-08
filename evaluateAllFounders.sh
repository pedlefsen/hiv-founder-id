#!/bin/bash
mkdir ${3}
for patient in  `ls -c1 ${1}/*.fasta | egrep --only "_[0-9]{5,}_" | tr -d "_" | uniq`
do
    export outputDir=${3}/${patient}
    mkdir ${outputDir}
   ./evaluateFounders.bash ${1} ${2} ${outputDir} ${patient} &
done

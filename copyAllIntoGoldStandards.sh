#!/bin/bash
mkdir ${1}/true_founders
for patient in  `ls -1 ${1}/*.fasta | egrep --only "_[0-9]{5,}_" | tr -d "_" | uniq`
do 
   ./copyIntoGoldStandards.bash ${1} $patient &
done

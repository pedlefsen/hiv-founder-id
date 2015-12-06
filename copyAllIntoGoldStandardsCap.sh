#!/bin/bash
mkdir ${1}/true_founders
for patient in  `ls -c1 ${1}/*.fasta | cut -d_ -f5 | sort -u | egrep "^[0-9]+" | uniq`
do 
   ./copyIntoGoldStandards.bash ${1} $patient &
done

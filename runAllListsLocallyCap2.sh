#!/bin/bash
mkdir ${2}
for patient in  `ls -c1 ${1}/*.fasta | cut -d_ -f5 | sort -u | egrep "^[0-9]+" | uniq`
do 
   ./runListLocally2.bash ${1} ${2} $patient &
done

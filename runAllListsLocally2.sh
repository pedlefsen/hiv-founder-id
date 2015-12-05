#!/bin/bash
for patient in  `ls -c1 ${1}/*.fasta | egrep --only "_[0-9]{5,}_" | tr -d "_"`
do 
   ./runListLocally2.bash ${1} $patient &
done

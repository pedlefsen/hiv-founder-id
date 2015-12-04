#!/bin/bash
for patient in  `ls -c1 ${1}/*.fasta | egrep --only "_[0-9]{5,}_" | tr -d "_"`
do 
   sbatch -Jhfi_${patient} -N1 -t 04:00 --mail-user=tholzman --mail-type=ALL ./runList2.bash ${1} $patient 
done

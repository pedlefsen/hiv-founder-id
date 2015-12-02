#!/bin/bash
for patient in  `ls -c1 ${1}/*.fasta | cut -d_ -f5 | sort -u | egrep "^[0-9]+"`
do 
   sbatch -Jhfi_${patient} -N1 -t 10:00 --mail-user=tholzman --mail-type=ALL ./runList2.bash ${1} $patient 
done

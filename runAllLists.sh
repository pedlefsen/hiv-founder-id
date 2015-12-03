#!/bin/bash
for patient in  `ls -c1 ${1}/*.fasta | cut -d_ -f6 | sort -u | egrep "^[0-9]+"`
do 
   sbatch -JhifCEPF${patient} -N1 -t 04:00 --mail-user=tholzman --mail-type=ALL ./runList.bash ${1} $patient 
done

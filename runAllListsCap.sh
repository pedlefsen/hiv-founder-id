#!/bin/bash
for patient in  `ls -c1 ${1}/*.fasta | cut -d_ -f5 | sort -u | egrep "^[0-9]+"`
do 
   sbatch -JhfiCEPF${patient} -N1 -t 04:00 --mail-user=tholzman --mail-type=FAIL ./runList.bash ${1} $patient 
done

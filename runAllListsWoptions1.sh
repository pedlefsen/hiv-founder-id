#!/bin/bash
export options=$2
for patient in `ls -c1 ${1}/*.fasta | egrep --only "_[0-9]{5,}_" | tr -d "_" | sort -u`
do 
    sbatch -Jhif.${patient}.${options} -N1 -t 04:00 --mail-user=tholzman --mail-type=FAIL ./runListWoptions.bash ${1} $patient $options 
done

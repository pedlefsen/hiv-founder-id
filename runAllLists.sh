#!/bin/bash
for patient in `ls -c1 ${1}/*.list  | egrep --only "[0-9]+\.list"  | tr -d "_" | sort -u`
do 
    sbatch -JhifCEPF${patient} -N1 -t 04:00 --mail-user=pedlefse --mail-type=FAIL ./runList.bash ${1} $patient 
done

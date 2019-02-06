#!/bin/bash
for patient in  `ls -c1 ${1}/*.list  | egrep --only "[0-9]+\.list" | egrep --only "[0-9]+"  | sort -u`
do 
   sbatch -Jhfi_${patient} -N1 -t 04:00 --mail-user=tholzman --mail-type=FAIL ./runList2.bash ${1} $patient 
done

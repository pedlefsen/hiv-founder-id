#!/bin/bash
mkdir ${1}/true_founders
for patient in  `ls -c1 ${1}/*.list  | egrep --only "[0-9]+\.list" | egrep --only "[0-9]+"  | sort -u`
do 
   ./copyIntoGoldStandards.bash ${1} $patient &
done

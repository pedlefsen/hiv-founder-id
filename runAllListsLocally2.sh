#!/bin/bash
echo $1
mkdir ${2}
for patient in  `ls -c1 ${1}/*.list  | egrep --only "[0-9]+\.list" | egrep --only "[0-9]+"  | sort -u`
do
    echo $patient
   ./runListLocally2.bash ${1} ${2} $patient &
done

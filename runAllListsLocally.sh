#!/bin/bash
echo $1
for patient in `ls -c1 ${1}/*.list  | egrep --only "[0-9]+\.list" | egrep --only "[0-9]+"  | sort -u`
do
    echo $patient
    ./runListLocally.bash ${1} $patient &
done

#!/bin/bash
echo $1
for patient in  `ls -c1 ${1}/*.list  | egrep --only "[0-9]+\.list" | egrep --only "[0-9]+"  | sort -u`
do
    echo $patient
   ./runListLocallyPhyML.bash ${1} $patient
done
## Now we also postprocess here because it has to wait for those to be done.
./postProcessIdentifyFounders.sh ${1} ${1} ${1}

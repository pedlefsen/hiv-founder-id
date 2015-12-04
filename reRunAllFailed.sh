#!/bin/bash
for patient in `find /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw/v3/1m -maxdepth 1 -name "*.err" -size +1000c -exec basename {} .err \;`
do 
   sbatch -JhifCEPF${patient} -N1 -t 04:00 --mail-user=tholzman --mail-type=FAIL ./runList.bash ${1} $patient 
done

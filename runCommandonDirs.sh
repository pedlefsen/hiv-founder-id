#!/bin/bash
baseDir=/fh/fast/edlefsen_p/bakeoff
for direc in ${baseDir}/analysis_sequences/raw/nflg/1m    ${baseDir}/analysis_sequences/raw/nflg/1m6m ${baseDir}/analysis_sequences/raw/nflg/6m ${baseDir}/gold_standard/caprisa_002/v3 ${baseDir}/gold_standard/rv217/nflg ${baseDir}/gold_standard/rv217/v3
do
   $1 $direc
done
#   ${baseDir}/analysis_sequences/raw/v3/1m                     
#   ${baseDir}/analysis_sequences/raw/v3/1m6m                   
#   ${baseDir}/analysis_sequences/raw/v3/6m                     


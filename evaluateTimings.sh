#!/bin/bash

### This is preparatory for evaluateTimings.R; run this on an RV217 or CAP/V3 dir to extract the dates from the file names, like so:
# ./evaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw/nflg/1m/ > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw/nflg/1m/sampleDates.tbl
# ./evaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw/v3/1m/ > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw/v3/1m/sampleDates.tbl

# Note that if you provide a second argument, it will look for "list files" called eg nnnnnnn.list instead of eg processed_nnnnnnn.list.
export mainDir=$1

if [ -z $2 ]; then
    export listPrefix="processed_"
else
    export listPrefix=""
fi

for patient in  `ls -c1 ${1}/*.fasta | egrep --only "_[0-9]{5,}_" | tr -d "_" | uniq`
do
    export listFile=${mainDir}/${listPrefix}${patient}.list
    
    export fastaFile=`head -n 1 ${listFile} | sed "s/\/home\/tholzman\/edf/\/fh\/fast\/edlefsen_p/g"`
    echo ${patient} `head -n 1 ${fastaFile} | head -n 1 | cut -d '|' -f 5`
done

#!/bin/bash

### This is preparatory for evaluateTimings.sh, evaluateFounders.bash, evaluateIsMultiple.R, etc.

# Arg one is a directory in which there are lists of "1w" fasta files
# named. nnnnnnn.list where nnnnnnn is a patient number.  Arg two is
# the analogous dir for the estimates. It expects directories within
# that directory to be called hiv_founder_id_processed_nnnnnnn unless
# a fourth argument is provided, in which case it will expect
# directories to be called hiv_founder_id_nnnnnnn. Arg three is the
# writable dir to put outputs, and the fourth argument if present
# toggles the expected dir name (as just described in the previous
# sentence).
#

# ./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw/nflg/1m/ 
export mainDir=$1
export estimateDir=$2
export outputDir=$3
export notprocessedFlag=$4

mkdir -p ${outputDir}

for patient in  `ls -c1 ${mainDir}/*.fasta | egrep --only "_[0-9]{5,}_" | tr -d "_" | uniq`
do
    if [ -z $4 ]; then 
    export inputDir=${estimateDir}/hiv_founder_id_processed_${patient}
    else
    export inputDir=${estimateDir}/hiv_founder_id_${patient}
    fi

    if [ !-e "${outputDir}/identify_founders.tab" ]
    then
        cat ${inputDir}/identify_founders.tab > ${outputDir}/identify_founders.tabatha
    else
        tail -n +2 ${inputDir}/identify_founders.tab >> ${outputDir}/identify_founders.tabatha
    fi
done

#### For now see README.postprocessing.txt

#!/bin/bash

### This is preparatory for evaluateTimings.sh, evaluateFounders.bash, evaluateIsMultiple.R, etc.
## See README.postprocessing.txt

# Arg one is a directory in which there are lists of "1w" fasta files
# named. nnnnnnn.list where nnnnnnn is a patient number.  Arg two is
# the analogous dir for the estimates. It expects directories within
# that directory to be called hiv_founder_id_processed_nnnnnnn unless
# a fourth argument is provided, in which case it will expect
# directories to be called hiv_founder_id_nnnnnnn. Arg three is the
# writable dir to put outputs, and the fourth argument if present
# toggles the expected dir name (as just described in the previous
# sentence). The fifth argument, if present, indicates that the analysis is of the partitions subdir with the fifth argument specifying the study (rv217 or caprisa002), the sixth argument specifying the time point (as 1M, 6M, or 1M6M), and the seventh argument specifying the region (only V3 is allowed fornow)
#

export mainDir=$1
export estimateDir=$2
export outputDir=$3
export notprocessedFlag=$4
export doPartitionsForStudy=$5
export doPartitionsForTime=$6
export doPartitionsForRegion=$7

mkdir -p ${outputDir}

if [ -e "${outputDir}/identify_founders.tab" ]; then
    rm "${outputDir}/identify_founders.tab"
fi

for patient in  `ls -c1 ${1}/*.list  | egrep --only "[0-9]+\.list" | egrep --only "[0-9]+"  | sort -u`
do
    echo "${patient}"
    if [ -z $5 ]; then
        export estimateDirMaybeModified=$estimateDir
    else
        export estimateDirMaybeModified="${estimateDir}/bakeoff_analysis_${doPartitionsForStudy}_${patient}_${doPartitionsForTime}_${doPartitionsForRegion}"
    fi
    if [ -z $4 ]; then 
        export inputDir=${estimateDirMaybeModified}/hiv_founder_id_processed_${patient}
    else
        export inputDir=${estimateDirMaybeModified}/hiv_founder_id_${patient}
    fi
    echo ${inputDir}

    if [ -e "${inputDir}/identify_founders.tab" ]; then
        if [ -e "${outputDir}/identify_founders.tab" ]; then
            # Skip the header since it is already there.
            echo "All but header: ${inputDir}/identify_founders.tab"
            tail -n +2 ${inputDir}/identify_founders.tab >> ${outputDir}/identify_founders.tab
        else
            echo "Including header: ${inputDir}/identify_founders.tab"
            cat ${inputDir}/identify_founders.tab > ${outputDir}/identify_founders.tab
        fi
    else
        echo "Missing input dir: ${inputDir}"
    fi
done

#### For now see README.postprocessing.txt

#!/bin/bash

## See README.postprocessing.txt

# Arg one is a directory in which there are lists of "1w" fasta files
# named. nnnnnnn.list where nnnnnnn is a patient number.  Arg two is the
# writable dir to put outputs and in which we expect to find subdirs named after the patients (nnnnn), each containing a file called evaluateFounders.tbl.
#

export mainDir=$1
export inputOutputDir=$2

if [ -e "${inputOutputDir}/evaluateFounders.tbl" ]; then
    rm "${inputOutputDir}/evaluateFounders.tbl"
fi

for patient in  `ls -c1 ${maindir}/*.list  | egrep --only "[0-9]+\.list" | egrep --only "[0-9]+"  | sort -u`
do
    export inputDir="${inputOutputDir}/${patient}"

    #echo ${patient}

    if [ -e "${inputDir}/evaluateFounders.tbl" ]; then
        #echo "file ${inputDir}/evaluateFounders.tbl exists"
        if [ -e "${inputOutputDir}/evaluateFounders.tbl" ]; then
            # Skip the header since it is already there.
            echo "All but header: ${inputDir}/evaluateFounders.tbl"
            tail -n +2 ${inputDir}/evaluateFounders.tbl >> ${inputOutputDir}/evaluateFounders.tbl
        else
            echo "Including header: ${inputDir}/evaluateFounders.tbl"
            cat ${inputDir}/evaluateFounders.tbl > ${inputOutputDir}/evaluateFounders.tbl
        fi
    else
        echo "Missing file: ${inputDir}/evaluateFounders.tbl"
    fi
done


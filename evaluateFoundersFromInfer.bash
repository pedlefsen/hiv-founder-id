#!/bin/bash
##
# Arg one is a directory in which there are lists of "1w" fasta files named. nnnnnnn.list
# where nnnnnnn is a patient number.  It expects a directories within that directory to be called 
# hiv_founder_id_processed_nnnnnnn. Arg two is the analogous dir for the estimates, arg three is the writable dir to put outputs, and the fourth argument is ptid nnnnnnn.
#
# D'OPTE 12/15
##
export mainDir=$1
export estimateDir=$2
mkdir ${3}
export outputDir=${3}/${patient}
rm -rf ${outputDir}
mkdir ${outputDir}
export patient=$4
export inputDir=${estimateDir}/founder-inference-bakeoff_${patient}
export truthDir=${mainDir}/true_founders/${patient}
export truthListFile=${mainDir}/processed_${patient}.list
export listFile=${estimateDir}/processed_${patient}.list
# Ok, so do all four kinds.  Append as we go.
export evaluateFounders_outputFilename="${outputDir}/evaluateFounders.tbl";
export evaluateFounders_append="TRUE";
for fasta_prefix in  `cat ${listFile} | rev | cut -d'/' -f-1 | rev | cut -d '.' -f 1`
do
    echo ${fasta_prefix}
    ## Fix the suffixes first. Note the fixed files are put in the _output_ directory.
    if [ -e "${inputDir}/${fasta_prefix}.outsingle.fa" ]; then
          cp "${inputDir}/${fasta_prefix}.outsingle.fa" "${outputDir}/${fasta_prefix}_outsingle.fa";
    fi
    if [ -e "${inputDir}/${fasta_prefix}.outmultiple.fa" ]; then
          cp "${inputDir}/${fasta_prefix}.outmultiple.fa" "${outputDir}/${fasta_prefix}_outmultiple.fa";
    fi
    export evaluateFounders_estimatesFilename_single="${outputDir}/${fasta_prefix}_outsingle.fa";
    export evaluateFounders_estimatesFilename_multiple="${outputDir}/${fasta_prefix}_outmultiple.fa";
    
    for truth_fasta_prefix in  `cat ${truthListFile} | rev | cut -d'/' -f-1 | rev | cut -d '.' -f 1`
    do
        echo ${truth_fasta_prefix}
        export evaluateFounders_truthsFilename_single="${truthDir}/${truth_fasta_prefix}_singlefounder.fasta";
        export evaluateFounders_truthsFilename_multiple="${truthDir}/${truth_fasta_prefix}_multifounder.fasta";
        
        export evaluateFounders_estimatesFilename=${evaluateFounders_estimatesFilename_single};
        export evaluateFounders_truthsFilename=${evaluateFounders_truthsFilename_single};
        R -f ./evaluateFounders.R --vanilla --slave
        
        export evaluateFounders_estimatesFilename=${evaluateFounders_estimatesFilename_single};
        export evaluateFounders_truthsFilename=${evaluateFounders_truthsFilename_multiple};
        R -f ./evaluateFounders.R --vanilla --slave
        
        export evaluateFounders_estimatesFilename=${evaluateFounders_estimatesFilename_multiple};
        export evaluateFounders_truthsFilename=${evaluateFounders_truthsFilename_single};
        R -f ./evaluateFounders.R --vanilla --slave
        
        export evaluateFounders_estimatesFilename=${evaluateFounders_estimatesFilename_multiple};
        export evaluateFounders_truthsFilename=${evaluateFounders_truthsFilename_multiple};
        R -f ./evaluateFounders.R --vanilla --slave
    done
done
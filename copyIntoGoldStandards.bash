#!/bin/bash
##
# Copies the identify-founders output into the true_founders directory
# (as a starting point for gold standards, possibly to be augmented by
# expert knowledge from Morgane and/or Carolyn's group).
# 
# Arg one is a directory in which there are lists of fasta files named. nnnnnnn.list
# where nnnnnnn is a patient number.  It expects a directories within that directory to be called 
# hiv_founder_id_nnnnnnn into which the results will go.  The second argument is nnnnnnn.
#
# D'OPTE 11/15
##
export mainDir=$1
export patient=$2
#export inputDir=./hiv_founder_id_processed_${patient}
export inputDir=${mainDir}/hiv_founder_id_processed_${patient}
export outputDir=${mainDir}/true_founders/${patient}
rm -rf ${outputDir}
mkdir ${outputDir}
export listFile=${mainDir}/processed_${patient}.list
# Cut out the path and cut out the .fasta suffix. NOTE THIS ASSUMES NO "/" IN PATHNAME.
for fasta_prefix in  `cat ${listFile} | rev | cut -d'/' -f-1 | rev | cut -d '.' -f 1`
do 
    cp -v ${inputDir}/${fasta_prefix}_singlefounder_cons.fasta ${outputDir}/${fasta_prefix}_singlefounder.fasta
    cp -v ${inputDir}/${fasta_prefix}_multiplefounders_cons.fasta ${outputDir}/${fasta_prefix}_multifounder.fasta 
done

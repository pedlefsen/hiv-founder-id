#!/bin/bash
##
# Submits Paul's code.  Arg one is a directory in which there are lists of fasta files named. nnnnnnn.list
# where nnnnnnn is a patient number.  It expects a directories within that directory to be called 
# hiv_founder_id_nnnnnnn into which the results will go.  The second argument is nnnnnnn.
#
# TAH 11/15
##
export mainDir=$1
export patient=$2
export outputDir=${mainDir}/hiv_founder_id_processed_${patient}
export outputFile=${outputDir}/identify_founders_PhyML.out
export errFile=${outputDir}/${patient}_PhyML.err
touch $errFile
export listFile=${mainDir}/processed_${patient}.list
perl -w ./identify_founders.pl -PTFEI -HR -t identify_founders_PhyML.tab -V -O ${outputDir}/ $listFile >$outputFile 2>>$errFile

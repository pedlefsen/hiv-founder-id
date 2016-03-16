#!/bin/bash
##
# Submits Paul's code.  Arg one is a directory in which there are lists of fasta files named. nnnnnnn.list
# where nnnnnnn is a patient number.  It expects a directories within that directory to be called 
# hiv_founder_id_nnnnnnn into which the results will go.  The second argument is nnnnnnn.
#
# The third argument is a list of options to pass to identify_founder
#
#
# TAH 3/16
##

## edit these commands to find the right libraries
export PATH=/fh/fast/edlefsen_p/tholzman/anaconda/bin:${PATH}
source activate root
export R_LIBS_USER=/home/tholzman/R/Library
export PERL5LIB=/home/tholzman/perl5/lib/perl5
##

export mainDir=$1
export patient=$2
export options=$3
export outputDir=${mainDir}/hiv_founder_id_${patient}
rm -rf ${outputDir}/*
export outputFile=${outputDir}/identify_founders.out
export errFile=${mainDir}/${patient}.err
export listFile=${mainDir}/${patient}.list
rm $errFile
touch $errFile
/usr/bin/time -a -o $errFile perl -w ./identify_founders.pl -${options} -O ${outputDir}/ $listFile >$outputFile 2>>$errFile

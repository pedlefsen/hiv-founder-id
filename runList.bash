#!/bin/bash
##
# Submits Paul's code.  Arg one is a directory in which there are lists of fasta files named. nnnnnnn.list
# where nnnnnnn is a patient number.  It expects a directories within that directory to be called 
# hiv_founder_id_nnnnnnn into which the results will go.  The second argument is nnnnnnn.
#
# TAH 11/15
##
export PATH=/fh/fast/edlefsen_p/tholzman/anaconda/bin:${PATH}
source activate root
export R_LIBS_USER=/home/tholzman/R/Library
export PERL5LIB=/home/tholzman/perl5/lib/perl5
export mainDir=$1
export patient=$2
export outputDir=${mainDir}/hiv_founder_id_${patient}
rm -rf ${outputDir}/*
export outputFile=${outputDir}/identify_founders.out
export errFile=${mainDir}/${patient}.err
export listFile=${mainDir}/${patient}.list
touch $errFile
/usr/bin/time -a -o $errFile perl -w ./identify_founders.pl  -E -P -F -O ${outputDir}/ $listFile >$outputFile 2>>$errFile

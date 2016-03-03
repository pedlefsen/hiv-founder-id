#!/bin/bash
##
# There is only one argument: a writable dir to put output files called [study]_evaluateIsMultiple.tab.
#
# D'OPTE 2/16
##
export evaluateIsMultiple_outputDir=$1

mkdir -p ${evaluateIsMultiple_outputDir}

# Ok, so do all four kinds.  Write one output file for each "study".
export evaluateIsMultiple_append="FALSE";
for study in rv217 caprisa002 rv217_v3;
do
    echo ${study};
    export evaluateIsMultiple_study="${study}"
    R -f ./evaluateIsMultiple.R --vanilla --slave
done

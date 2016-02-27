#!/bin/bash
##
# There is only one argument: a writable dir to put output files called [study]_evaluateIsMultiple.tab.
#
# D'OPTE 2/16
##
export evaluateIsMultiple_outputDir=$1
#rm -rf ${outputDir}
mkdir ${outputDir}
# Ok, so do all four kinds.  Append as we go.
export evaluateIsMultiple_append="FALSE";
for study in rv217 caprisa002;
do
    echo ${study};
    rm "${outputDir}/${study}_evaluateIsMultiple.tbl"
    export evaluateIsMultiple_study="${study}"
    R -f ./evaluateIsMultiple.R --vanilla --slave
done

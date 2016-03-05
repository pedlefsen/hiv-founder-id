#!/bin/bash
##
# Arg one is a directory (the gold standard dir) in which there are
# lists of "1w" fasta files named nnnnnnn.list where nnnnnnn is a
# patient number.  It expects a directories within that directory to
# be called hiv_founder_id_processed_nnnnnnn unless a fifth argument
# is provided, in which case it will expect directories to be called
# hiv_founder_id_nnnnnnn. Arg two is the analogous dir for the
# estimates (where unless a fifth arg is given, there should be list
# files called processed_nnnnn.list; otherwise just nnnnn.list), arg
# three is the writable dir to put outputs, and the fourth argument is
# ptid nnnnnnn.  The fifth argument if present toggles the expected
# dir name and list file name (as just described in the previous
# sentences).
#
# D'OPTE 12/15
# Modified 2/16
##
export mainDir=$1
export estimateDir=$2
export outputDir=$3
export patient=$4
if [ -z $5 ]; then 
export inputDir=${estimateDir}/hiv_founder_id_processed_${patient}
else
export inputDir=${estimateDir}/hiv_founder_id_${patient}
fi
export truthDir=${mainDir}/true_founders/${patient}
# rm -rf ${outputDir}
mkdir ${outputDir}
export truthListFile=${mainDir}/processed_${patient}.list
if [ -z $5 ]; then 
    export listFile=${estimateDir}/processed_${patient}.list
else
    export listFile=${estimateDir}/${patient}.list
fi
# For V3 we do just V3 and get output for just ENV.
export evaluateFounders_proteinsList="V3";
export evaluateFounders_genomeRegion="V3";
# Ok, so do all four kinds.  Append as we go.
export evaluateFounders_outputFilename="${outputDir}/evaluateFounders.tbl";
export evaluateFounders_append="FALSE";
rm ${evaluateFounders_outputFilename};
for fasta_prefix in  `cat ${listFile} | rev | cut -d'/' -f-1 | rev | cut -d '.' -f 1`
do
    echo ${fasta_prefix};

    export evaluateFounders_estimatesFilename_single="${inputDir}/${fasta_prefix}_singlefounder_cons.fasta";
    export evaluateFounders_estimatesFilename_multiple="${inputDir}/${fasta_prefix}_multiplefounders_cons.fasta";
    
    for truth_fasta_prefix in  `cat ${truthListFile} | rev | cut -d'/' -f-1 | rev | cut -d '.' -f 1`
    do
        echo ${truth_fasta_prefix}
        export evaluateFounders_truthsFilename_single="${truthDir}/${truth_fasta_prefix}_singlefounder.fasta";
        export evaluateFounders_truthsFilename_multiple="${truthDir}/${truth_fasta_prefix}_multifounder.fasta";
        
        export evaluateFounders_estimatesFilename=${evaluateFounders_estimatesFilename_single};
        export evaluateFounders_truthsFilename=${evaluateFounders_truthsFilename_single};
        R -f ./evaluateFounders.R --vanilla --slave
        
        export evaluateFounders_append="TRUE";
        
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

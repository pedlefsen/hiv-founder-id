#!/bin/bash
export R_LIBS_USER=~/R/Library
export workdir=$1
export recombdir=$1/recomb
for f in `find ${workdir} -maxdepth 2 -name "*.fasta"`
do
   export removeDuplicateSequencesFromAlignedFasta_inputFilename=$f
   export removeDuplicateSequencesFromAlignedFasta_outputFastaFilename=`basename $f .fasta`.dedup.fasta
   export removeDuplicateSequencesFromAlignedFasta_outputDir=$recombdir
   export removeDuplicateSequencesFromAlignedFasta_outputTableFilename=`basename $f .fasta`.dedup.tab
   echo "source('removeDuplicateSequencesFromAlignedFasta.R')" | R --vanilla --slave
done
pushd .
cd $recombdir
ls -c1 *.dedup.fasta | nl -w3 -nrz | sed 's/^\([0-9][0-9][0-9]\)\s*\(.*\)/ln -s \2 \1\.fasta/' | sh
rm `grep -c ">" [0-9][0-9][0-9].fasta | egrep ':[1,2]$' | cut -d: -f1`
popd 


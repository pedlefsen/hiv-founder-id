#!/bin/sh
##
# takes one argument.  The root of the directory tree containing fasta files.
# finds directories with fasta files.  Makes a list of fasta files and a
# directory called "hiv_founder_results
#
# TAH 11/15 
##
export root="$1"
export tmpf=$(tempfile)
echo $root
echo $tmpf
find $root -name "*.fasta" -exec dirname {} \; | sort -u >> $tmpf
while read dire 
do
   echo $dire
   find $dire -name "*.fasta" > ${dire}/fasta.list
   mkdir ${dire}/hiv_founder_results
done <$tmpf
rm $tmpf


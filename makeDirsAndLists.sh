#!/bin/sh
##
# Make a list file and a directory for each patient
#
# TAH 11/15
##

for patient in `ls -c1 ./*.fasta | egrep --only "_[0-9]{5,}_" | tr -d "_"`
do
   rm -rf hiv_founder_id_$patient
   mkdir hiv_founder_id_$patient
   rm ${patient}.list
   touch ${patient}.list
   for f in `ls -c1 *_${patient}_*.fasta`
   do 
      echo `pwd`/$f >> ${patient}.list
   done
done

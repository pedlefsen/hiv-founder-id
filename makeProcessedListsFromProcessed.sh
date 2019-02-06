#!/bin/sh
##
# read directories of the structure hiv_founder_id_nnnnnn.  Find the file called
# identify_founders.tab.  Find the 4th field of 2nd through last lines of that
# file, and put the full path of that file into a file called 
# processed_nnnnnn.list.
#
# TAH 11/15
##
for f in  `find hiv_founder_id_processed_* -name identify_founders.tab`
do
   patient=`echo $f |  cut -d_ -f5 | cut -d"/" -f1 | sort -u`
   rm processed_${patient}.list 
   touch processed_${patient}.list
   for g in `tail -n +2 $f | cut -f2`
   do 
      echo `pwd`/hiv_founder_id_processed_${patient}/$g >> processed_${patient}.list
   done
done

      

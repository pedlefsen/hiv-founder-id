 cd /home/phillipl/projects/hiv-founder-id/code/hiv-founder-id
rm -r tests/tmp/test_001
perl identify_founders.pl -CPRFHT -o tests/tmp/test_001 tests/data/ld_seqs.fasta

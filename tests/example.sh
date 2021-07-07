# run identify_founders.pl

mkdir /tmp/hf_example

cd /home/phillipl/projects/hiv-founder-id/code/hiv-founder-id

perl /home/phillipl/projects/hiv-founder-id/code/hiv-founder-id/identify_founders.pl -PRT -o /tmp/hf_example/ld_seqs_0 /home/phillipl/projects/hiv-founder-id/code/hiv-founder-id/tests/example_data/ld_seqs_0.fasta
perl /home/phillipl/projects/hiv-founder-id/code/hiv-founder-id/identify_founders.pl -PRT -o /tmp/hf_example/ld_seqs_1 /home/phillipl/projects/hiv-founder-id/code/hiv-founder-id/tests/example_data/ld_seqs_1.fasta
perl /home/phillipl/projects/hiv-founder-id/code/hiv-founder-id/identify_founders.pl -PRT -o /tmp/hf_example/ld_seqs_2 /home/phillipl/projects/hiv-founder-id/code/hiv-founder-id/tests/example_data/ld_seqs_2.fasta
perl /home/phillipl/projects/hiv-founder-id/code/hiv-founder-id/identify_founders.pl -PRT -o /tmp/hf_example/ld_seqs_3 /home/phillipl/projects/hiv-founder-id/code/hiv-founder-id/tests/example_data/ld_seqs_3.fasta
perl /home/phillipl/projects/hiv-founder-id/code/hiv-founder-id/identify_founders.pl -PRT -o /tmp/hf_example/ld_seqs_4 /home/phillipl/projects/hiv-founder-id/code/hiv-founder-id/tests/example_data/ld_seqs_4.fasta

# collate identify_founders.tab files
cp /tmp/hf_example/ld_seqs_0/identify_founders.tab /tmp/hf_example/identify_founders.tab
tail -n 1 /tmp/hf_example/ld_seqs_1/identify_founders.tab >> /tmp/hf_example/identify_founders.tab
tail -n 1 /tmp/hf_example/ld_seqs_2/identify_founders.tab >> /tmp/hf_example/identify_founders.tab
tail -n 1 /tmp/hf_example/ld_seqs_3/identify_founders.tab >> /tmp/hf_example/identify_founders.tab
tail -n 1 /tmp/hf_example/ld_seqs_4/identify_founders.tab >> /tmp/hf_example/identify_founders.tab

# run estimateInfectionTime.R

/home/phillipl/projects/hiv-founder-id/code/hiv-founder-id/estimateInfectionTime.R --model_structure=slope --identify_founders.tab=/tmp/hf_example/identify_founders.tab --estimator=pfitter

./estimateInfectionTime.R --model_structure=full --identify_founders.tab=/tmp/hf_example/identify_founders.tab --vl_file=tests/example_data/vl_diff.csv --bounds_file=tests/example_data/bounds_diff.csv --estimator=pfitter

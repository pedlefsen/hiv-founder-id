# run identify_founders.pl

cd /home/docker/hiv-founder-id

perl /home/docker/hiv-founder-id/identify_founders.pl -PRT -o /home/docker/example/ld_seqs_0 /home/docker/example/ld_seqs_0.fasta
perl /home/docker/hiv-founder-id/identify_founders.pl -PRT -o /home/docker/example/ld_seqs_1 /home/docker/example/ld_seqs_1.fasta
perl /home/docker/hiv-founder-id/identify_founders.pl -PRT -o /home/docker/example/ld_seqs_2 /home/docker/example/ld_seqs_2.fasta
perl /home/docker/hiv-founder-id/identify_founders.pl -PRT -o /home/docker/example/ld_seqs_3 /home/docker/example/ld_seqs_3.fasta
perl /home/docker/hiv-founder-id/identify_founders.pl -PRT -o /home/docker/example/ld_seqs_4 /home/docker/example/ld_seqs_4.fasta

# collate identify_founders.tab files
cp /home/docker/example/ld_seqs_0/identify_founders.tab /home/docker/example/identify_founders.tab
tail -n 1 /home/docker/example/ld_seqs_1/identify_founders.tab >> /home/docker/example/identify_founders.tab
tail -n 1 /home/docker/example/ld_seqs_2/identify_founders.tab >> /home/docker/example/identify_founders.tab
tail -n 1 /home/docker/example/ld_seqs_3/identify_founders.tab >> /home/docker/example/identify_founders.tab
tail -n 1 /home/docker/example/ld_seqs_4/identify_founders.tab >> /home/docker/example/identify_founders.tab

# run estimateInfectionTime.R

/home/docker/hiv-founder-id/estimateInfectionTime.R --model_structure=slope --identify_founders.tab=/home/docker/example/identify_founders.tab --estimator=pfitter

/home/docker/hiv-founder-id/estimateInfectionTime.R --model_structure=full --identify_founders.tab=/home/docker/example/identify_founders.tab --vl_file=/home/docker/example/vl_diff.csv --bounds_file=/home/docker/example/bounds_diff.csv --estimator=pfitter

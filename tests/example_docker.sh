cd /home/docker/hiv-founder-id

perl /home/docker/hiv-founder-id/identify_founders.pl -PRT -o /home/docker/example/example_1 /home/docker/example/example_1.fasta
perl /home/docker/hiv-founder-id/identify_founders.pl -PRT -o /home/docker/example/example_2 /home/docker/example/example_2.fasta
perl /home/docker/hiv-founder-id/identify_founders.pl -PRT -o /home/docker/example/example_3 /home/docker/example/example_3.fasta
perl /home/docker/hiv-founder-id/identify_founders.pl -PRT -o /home/docker/example/example_4 /home/docker/example/example_4.fasta
perl /home/docker/hiv-founder-id/identify_founders.pl -PRT -o /home/docker/example/example_5 /home/docker/example/example_5.fasta
perl /home/docker/hiv-founder-id/identify_founders.pl -PRT -o /home/docker/example/example_6 /home/docker/example/example_6.fasta
perl /home/docker/hiv-founder-id/identify_founders.pl -PRT -o /home/docker/example/example_7 /home/docker/example/example_7.fasta

# collate identify_founders.tab files
cp        /home/docker/example/example_1/identify_founders.tab    /home/docker/example/identify_founders.tab
tail -n 1 /home/docker/example/example_2/identify_founders.tab >> /home/docker/example/identify_founders.tab
tail -n 1 /home/docker/example/example_3/identify_founders.tab >> /home/docker/example/identify_founders.tab
tail -n 1 /home/docker/example/example_4/identify_founders.tab >> /home/docker/example/identify_founders.tab
tail -n 1 /home/docker/example/example_5/identify_founders.tab >> /home/docker/example/identify_founders.tab
tail -n 1 /home/docker/example/example_6/identify_founders.tab >> /home/docker/example/identify_founders.tab
tail -n 1 /home/docker/example/example_7/identify_founders.tab >> /home/docker/example/identify_founders.tab

# run estimateInfectionTime.R

/home/docker/hiv-founder-id/estimateInfectionTime.R --model_structure=slope --identify_founders.tab=/home/docker/example/identify_founders.tab --estimator=pfitter

/home/docker/hiv-founder-id/estimateInfectionTime.R --model_structure=slope --identify_founders.tab=/home/docker/example/identify_founders.tab --bounds_file=/home/docker/example/bounds_diff.csv --estimator=pfitter

/home/docker/hiv-founder-id/estimateInfectionTime.R --model_structure=full --identify_founders.tab=/home/docker/example/identify_founders.tab --vl_file=/home/docker/example/vl_diff.csv --bounds_file=/home/docker/example/bounds_diff.csv --estimator=pfitter

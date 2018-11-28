#!/usr/bin/Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

#make_option("--non_overlapping", 
#            action = "store_true",
#            default = FALSE,
#            help = "If your forwards and reverse reads do not overlap, include this flag. Also, if they do overlap, but you do not want them to be merged, include this flag."),

#make_option("--max_seq",
#            default = 0,
#            help = paste("Maximum number of sequences to read from the input files. By default the ",
#                         "forward and reverse read files produced by a MiSeq are in the same order, ",
#                         "so just reading out the first 'max_seq' sequences from each file should ",
#                         "yield 'max_seq' pair-end reads in the dataset.",
#                         sep = "")),

make_option("--pipeline-dir", help = "Path to the root of the pipeline folder"),

make_option('--print-spec',
            action = 'store_true',
            help = 'Prints out a list of all the tests in the spec file')

)

opt <- parse_args(OptionParser(option_list = option_list,
  description = 'Script to control the testing for the hiv-founder-pipeline',
  epilogue = 'Example Calls:
./test_pipeline.R --pipeline-dir="/home/phillipl/projects/hiv-founder-id/code/hiv-founder-id" --print-spec
')
)

print(opt)


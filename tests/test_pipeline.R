#!/usr/bin/Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(knitr))

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

make_option("--pipeline_dir", help = "Path to the root of the pipeline folder"),

make_option('--print_spec',
            action = 'store_true',
            help = 'Prints out a list of all the tests in the spec file'),

make_option('--compare_spec_to_command_scripts',
            action = 'store_true',
            help = 'Compares the command scripts in the commands folder to the test specifications in the test_spec.csv file')
)

opt <- parse_args(OptionParser(option_list = option_list,
  description = 'Script to control the testing for the hiv-founder-pipeline',
  epilogue = 'Example Calls:
./test_pipeline.R --pipeline-dir="/home/phillipl/projects/hiv-founder-id/code/hiv-founder-id" --print-spec
')
)

print(opt)

pipeline_dir <- opt$pipeline_dir
if (!dir.exists(pipeline_dir)){
  stop('ERROR: Invalid dir passed to --pipeline_dir')
}

test_utils_file_name <- paste(pipeline_dir, '/tests/dev/test_utils.R', sep = '')
if (!file.exists(test_utils_file_name)){
  stop(paste('ERROR: The test utilities script, "', test_utils_file_name, '", does not exists', sep = ''))
}
source(test_utils_file_name)

test_spec_file_name <- paste(pipeline_dir, '/tests/test_specs.csv', sep = '')
if (!file.exists(test_spec_file_name)){
  stop(paste('ERROR: The test utilities script, "', test_spec_file_name, '", does not exists', sep = ''))
}

if (opt$print_spec){
  cat('
Test Specification File:
========================')
  test_spec <- read_test_spec_file(pipeline_dir = pipeline_dir)
  kable(test_spec)
}



#!/usr/bin/Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(knitr))

option_list <- list(

make_option("--pipeline_dir", help = "Path to the root of the pipeline folder"),

make_option('--print_spec',
            action = 'store_true',
            default = FALSE,
            help = 'Prints out a list of all the tests in the spec file'),

make_option('--compare_spec_to_command_scripts',
            action = 'store_true',
            default = FALSE,
            help = 'Compares the command scripts in the commands folder to the test specifications in the test_spec.csv file'),

make_option('--run_all_test_command_scripts',
            action = 'store_true',
            default = FALSE,
            help = 'Generate the call that will produce all the pipeline runs required to perform the testing')

)

opt <- parse_args(OptionParser(option_list = option_list,
  description = 'Script to control the testing for the hiv-founder-pipeline',
  epilogue = 'Example Calls:
./test_pipeline.R --pipeline-dir="/home/phillipl/projects/hiv-founder-id/code/hiv-founder-id" --print-spec
'))

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

if (opt$run_all_test_command_scripts){
  cat(paste('
Due to a weird interaction between R, bash and the pipeline, this script cannot execute the appropriate calls. Copy and paste the following command into a terminal to perform all the runs.
WARNING: This may take a long time.

cd ', pipeline_dir, '/tests/commands; ./run_all.sh

', sep = ''))
}


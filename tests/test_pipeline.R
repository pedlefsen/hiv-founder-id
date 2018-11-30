#!/usr/bin/Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(rmarkdown))

option_list <- list(

make_option('--compact_help',
            action = 'store_true',
            default = FALSE,
            help = 'Prints out compact list of viable options'),

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
            help = 'Generate the call that will produce all the pipeline runs required to perform the testing'),

make_option('--verbose',
            action = 'store_true',
            default = FALSE,
            help = 'Generate verbose output'),

make_option('--build_command_scripts',
            action = 'store_true',
            default = FALSE,
            help = 'Repopulate the commands folder based on the contents of the test_spec.csv file'),

make_option('--check_freshness',
            action = 'store_true',
            default = FALSE,
            help = 'Computes the number of hours that has elapsed since the last time each test specified in the test_spec file was run'),

make_option('--build_test_doc',
            action = 'store_true',
            default = FALSE,
            help = 'Builds the knitr document that presents the test results')

)

opt <- parse_args(OptionParser(option_list = option_list,
  description = 'Script to control the testing for the hiv-founder-pipeline',
  epilogue = 'Example Calls:
./test_pipeline.R --pipeline-dir="/home/phillipl/projects/hiv-founder-id/code/hiv-founder-id" --print-spec
'))

if (opt$compact_help){
  for (i in option_list){print(i@long_flag)}
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

print(opt)

pipeline_dir <- opt$pipeline_dir
if (!dir.exists(pipeline_dir)){
  stop('ERROR: Invalid dir passed to --pipeline_dir')
}

test_utils_file_name <- paste(pipeline_dir, '/tests/test_utils.R', sep = '')
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
  print(kable(test_spec))
}

if (opt$run_all_test_command_scripts){
  cat(paste('
Due to a weird interaction between R, bash and the pipeline, this script cannot execute the appropriate calls. Copy and paste the following command into a terminal to perform all the runs.
WARNING: This may take a long time.

cd ', pipeline_dir, '/tests/commands; bash run_all.sh

', sep = ''))
}

if (opt$compare_spec_to_command_scripts){
  cat('
Comparison of command scripts and test specification file:
==========================================================
')
  comparison_result <- compare_command_scripts_and_spec_file(pipeline_dir, verbose = opt$verbose)
  print(kable(comparison_result))
  if (any(comparison_result$result == 'mis_match')){
    cat('
Mismatches found: 
Common cause is tilde expansion - inspect this by setting the "--verbose" option.
If tilde expansion is not the case, you should almost certainly rebuild the command scripts with the "--build_command_scripts" option')
  }
}

if (opt$build_command_scripts){
  cat('
Rebuilding all the scripts in the commands folder
=================================================
The run_all.sh script was produced.
Additionally, scripts were produced for the following tests:

')
  scripts_produced <- write_all_commands_from_spec(pipeline_dir)
  write_run_all_command_script(pipeline_dir)
  print(scripts_produced)
}

if (opt$check_freshness){
  cat('
Checking that all the commands were recently executed
')
  last_run_times <- all_last_run_times(pipeline_dir)
  print(kable(last_run_times))
}

if (opt$build_test_doc){
  cat('
Building knitr document
')
  render(paste(pipeline_dir, '/tests/dev_tests.Rmd', sep = ''), 
         output_dir = paste(pipeline_dir, '/tests', sep = ''))
}

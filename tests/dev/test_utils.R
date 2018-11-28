#' Read test specification file
#'
#' Reads the test specification file and optionally filters it according to the given test_name
#' @export

read_test_spec_file <- function(pipeline_dir, test_name = NULL){
  spec_file <- read.csv(paste(pipeline_dir, '/tests/test_specs.csv', sep = ''),
                        stringsAsFactors = FALSE)
  if (!is.null(test_name)){
    return(spec_file[spec_file$test_name == test_name, ])
  }
  return(spec_file)
}

#' Generate command
#'
#' Given the input parameters for a command, generate the command character vectors that will eventually get written into the command script.
#' @export

generate_command <- function(test_name, fasta_file, command_flags, 
                             test_description, pipeline_dir){
  command <- paste('cd ', pipeline_dir, '
rm -r ', pipeline_dir, '/tests/tmp/', test_name, '
perl ', pipeline_dir, '/identify_founders.pl ', command_flags, ' -o ', pipeline_dir, '/tests/tmp/', test_name, ' ', pipeline_dir, '/tests/data/', fasta_file, 
  sep = '')
  return(command)
}

#' Read command file
#'
#' Given the name of the test and the base folder, read the command script into
#' a character vector.  
#' @export

read_command_script <- function(pipeline_dir, test_name){
  command_file <- paste(pipeline_dir, '/tests/commands',
                        '/', test_name, '.sh', sep = '')
  command_from_file <- read.delim(command_file, header = F, stringsAsFactors = FALSE)
  command_from_file <- paste(command_from_file[,1], collapse = '\n', sep = '\n')
  return(command_from_file)
}

#' Compares the command scripts to the test spec file
#'
#' Read through the test spec file and generate the related command for each
#' test. This is compared to the command scripts in the commands folder and any
#' scripts that do not match the generated command are returned.
#'
#' @export

compare_command_scripts_and_spec_file <- function(pipeline_dir, verbose = FALSE){
  test_specs <- read_test_spec_file(pipeline_dir = pipeline_dir)
  comparison_results <- NULL

  for (indx in 1:nrow(test_specs)){
    c_test_name <- test_specs[indx, 'test_name']
    c_fasta_file <- test_specs[indx, 'fasta_file']
    c_command_flags <- test_specs[indx, 'command_flags']
    c_test_description <- test_specs[indx, 'test_description']

    generated_command <- generate_command(test_name = c_test_name,
                                          fasta_file = c_fasta_file,
                                          command_flags = c_command_flags,
                                          test_description = c_test_description,
                                          pipeline_dir = pipeline_dir)

    command_from_file <-
      read_command_script(pipeline_dir = pipeline_dir,
                          test_name = c_test_name)

    if (generated_command != command_from_file){
      c_result <- 'mis_match'
    } else {
      c_result <- 'match'
    }

    if (verbose){
      print(paste('Test name: ', c_test_name, sep = ''))
      print(generated_command)
      print(command_from_file)
    }

    comparison_results <- rbind(comparison_results,
      data.frame(test_name = c_test_name,
                 result = c_result,
                 stringsAsFactors = FALSE)
      )
  }

  return(comparison_results)
}

#' Write a command script
#'
#' Given a command (character vector of the commands, as produced by generate_command), write a bash script into the commands folder.
#' @export

write_command_script <- function(pipeline_dir, test_name, command){
  command_file <- paste(pipeline_dir, '/tests/commands',
                        '/', test_name, '.sh', sep = '')
  write(command, command_file)
}

#' Write a script that will run all the command scripts
#'
#' Produces a scripts in the commands folder that will run all the command scripts if it is called with bash.
#' @export

write_run_all_command_script <- function(pipeline_dir){
  command_file <- paste(pipeline_dir, '/tests/commands',
                        '/run_all.sh', sep = '')
  command <- 'for f in test_*.sh; do
  echo "$f"
  bash "$f" -H
done'
  write(command, command_file)
}

#' Generate command scripts from test_spec file
#'
#' Reads the test_spec file and generates a command script for each test in the commands folder
#' @export

write_all_commands_from_spec <- function(pipeline_dir){
  test_specs <- read_test_spec_file(pipeline_dir = pipeline_dir)
  test_scripts_produced <- NULL
  for (indx in 1:nrow(test_specs)){
    c_test_name <- test_specs[indx, 'test_name']
    c_fasta_file <- test_specs[indx, 'fasta_file']
    c_command_flags <- test_specs[indx, 'command_flags']
    c_test_description <- test_specs[indx, 'test_description']

    generated_command <- generate_command(test_name = c_test_name,
                                          fasta_file = c_fasta_file,
                                          command_flags = c_command_flags,
                                          test_description = c_test_description,
                                          pipeline_dir = pipeline_dir)

    write_command_script(pipeline_dir = pipeline_dir,
                         test_name = c_test_name,
                         command = generated_command)
    test_scripts_produced <- c(test_scripts_produced, c_test_name)
  }
  return(test_scripts_produced)
}

#' Retrieves the time since the last run of the command scripts
#'
#' Inspects the timestamps on the output from the command scripts and computes the number of hours since the last time the commands were run.
#' @export

last_run_time <- function(pipeline_dir, test_name){
  result_folder <- paste(pipeline_dir, '/tests/tmp',
                         '/', test_name, sep = '')
  last_run <- min(file.info(dir(result_folder, full.names = TRUE))$ctime)
  time_since_last_run <- as.numeric(difftime(Sys.time(), last_run), units = 'hours')
  return(time_since_last_run)
}

#' Computes the elapsed time since the previous time the test_spec was run
#'
#' Load the test spec file and computes the last time each specified test was run.
#' @export

all_last_run_times <- function(pipeline_dir){
  test_spec <- read_test_spec_file(pipeline_dir)
  last_run_times <- NULL
  for (c_test_name in test_spec$test_name){
    last_run_times <- rbind(last_run_times,
      data.frame(test_name = c_test_name,
                 time_since_last_run = last_run_time(pipeline_dir = pipeline_dir,
                                                     test_name = c_test_name),
                 stringsAsFactors = FALSE)
      )
  }
  return(last_run_times)
}

#' Conducts a test on the hiv-founder-pipeline
#'
#' Given a test name and specification, check that the command to run the test
#' matches the command file stored in the repo and then test that the result
#' from running that command file matches the stored reference values.
#'
#' @export

conduct_test <- function(c_test_name, fasta_file, command_flags, test_description, pipeline_dir, write_new_command = FALSE){
  c_result <- NULL
  
  ref_folder <- paste(pipeline_dir, '/tests/ref_results',
                      '/', c_test_name, sep = '')
  result_folder <- paste(pipeline_dir, '/tests/tmp',
                         '/', c_test_name, sep = '')
  command_file <- paste(pipeline_dir, '/tests/commands',
                        '/', c_test_name, '.sh', sep = '')
  
  ref_file_names <- list.files(ref_folder)
  result_file_names <- list.files(result_folder)

  time_since_last_run <- last_run_time(pipeline_dir = pipeline_dir,
                                       test_name = c_test_name)
  if (time_since_last_run > 20){
    return('Time since last run is more than 20 hours, rerun tests')
  }
  c_result <- rbind(c_result,
    data.frame(name = 'time_since_last_run',
               result = TRUE,
               obs = time_since_last_run,
               expected = '< 20',
               note = 'Tests must be run daily',
               stringsAsFactors = FALSE)
    )

  command <- generate_command(test_name = c_test_name, 
                              fasta_file = fasta_file,
                              command_flags = command_flags,
                              test_description = test_description,
                              pipeline_dir = pipeline_dir)

  if (write_new_command){
    write_command_script(pipeline_dir, c_test_name, command)
    stop('New command written, now rerun bash scripts')
  }

  command_from_file <- read_command_script(pipeline_dir = pipeline_dir,
                                           test_name = c_test_name)

  commands_match <- command == command_from_file
  if (!commands_match){
    stop('Constructed command and command from command script do not match')
  }

  c_result <- rbind(c_result,
    data.frame(name = 'output_folder_exists',
               result = dir.exists(paste(pipeline_dir, '/tests/tmp/', c_test_name, sep = '')),
               obs = NA_character_,
               expected = NA_character_,
               note = 'We expect the output folder to exists after calling the pipeline',
               stringsAsFactors = FALSE)
    )

  for (c_file in ref_file_names){
    c_result <- rbind(c_result,
      data.frame(name = paste(c_file, '_exists', sep = ''),
                 result = c_file %in% result_file_names,
                 obs = NA_character_,
                 expected = NA_character_,
                 note = 'We expect all the required files to have been created',
                 stringsAsFactors = FALSE)
      )
  }

  for (c_file in ref_file_names){
    ref_md5 <- md5sum(paste(ref_folder, '/', c_file, sep = ''))
    result_md5 <- md5sum(paste(result_folder, '/', c_file, sep = ''))

    names(ref_md5) <- NULL
    names(result_md5) <- NULL

    c_result <- rbind(c_result,
      data.frame(name = paste(c_file, '_matches_ref', sep = ''),
                 result = c_file %in% result_file_names,
                 obs = result_md5,
                 expected = ref_md5,
                 note = 'Produced files are expected to be identical to the reference files',
                 stringsAsFactors = FALSE)
      )
  }

  c_result$name <- gsub("(.{15})", "\\1 ", c_result$name)
  c_result$obs <- gsub("(.{17})", "\\1 ", c_result$obs)
  c_result$expected <- gsub("(.{17})", "\\1 ", c_result$expected)

  return(list(test_name = c_test_name,
              fasta_file = fasta_file,
              command_flags = command_flags,
              test_description = test_description,
              command = command,
              test_result = c_result)
  )
}

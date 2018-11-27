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

  last_run <- min(file.info(dir(result_folder, full.names = TRUE))$ctime)
  time_since_last_run <- as.numeric(difftime(Sys.time(), last_run), units = 'hours')
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
    write(command, command_file)
    stop('New command written, now rerun bash scripts')
  }

#  command_from_file <- read.delim(command_file, header = F, stringsAsFactors = FALSE)
#  command_from_file <- paste(command_from_file[,1], collapse = '\n', sep = '\n')
  command_from_file <- read_command_script(pipeline_dir = pipeline_dir,
                                           test_name = c_test_name)

  commands_match <- command == command_from_file
  if (!commands_match){
    print(command)
    print('-----')
    print(command_from_file)
    print('-----')
    return(paste('Commands do not match:
Command 1:
', command, '
-------------
Command 2:
', command_from_file, '
', sep = ''))
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

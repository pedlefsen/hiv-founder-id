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

  command <- paste('cd ', pipeline_dir, '
rm -r tests/tmp/', c_test_name, '
perl identify_founders.pl ', command_flags, ' -o tests/tmp/', c_test_name, ' tests/data/', fasta_file, 
  sep = '')

  if (write_new_command){
    write(command, command_file)
    stop('New command written, now rerun bash scripts')
  }

  file_command <- read.delim(command_file, header = F, stringsAsFactors = FALSE)
  file_command <- paste(file_command[,1], collapse = '\n', sep = '\n')

  commands_match <- command == file_command
  if (!commands_match){
    print(command)
    print('-----')
    print(file_command)
    print('-----')
    return(paste('Commands do not match:
Command 1:
', command, '
-------------
Command 2:
', file_command, '
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

  return(list(test_name = c_test_name,
              fasta_file = fasta_file,
              command_flags = command_flags,
              test_description = test_description,
              command = command,
              test_result = c_result)
  )
}

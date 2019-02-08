#!/usr/bin/Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option("--input_file", 
              help = "Path to input file and file name"),
  make_option("--output_file", 
              help = "Path and name of output file"),
  make_option("--reduplicate",
              action = "store_true",
              default = FALSE,
              help = "Instead of removing duplicates, consult an existing duplicates_file and re-introduce previously removed duplicates."),
  make_option("--duplicates_file",
              default = NULL,
              help = "The name of the file that will track the names of the duplicates that were removed. The format is a csv file with the first column the name of the sequences that was left in the file and the second column containing a the name of a SINGLE sequence that was removed. If you are re-inserting duplicates, this must already exist and will be read as input. If left blank, then the .fasta extension of the output_file will be replaced with .csv. If --reduplicate is specified, then the .fasta extension of the input_file will be replaced with .csv."),
  make_option("--verbose",
              default = FALSE,
              action = "store_true",
              help = "Print verbose output to STDOUT."),
  make_option("--but_first_remove",
              default = NULL,
              help = "A tab seperated file with a header row and with the first column specifying the names of the sequences that must be removed before the reduplication operation.")
)

if (FALSE){
  # some debugging code that does not get run

  opt <- list(input_file = "tests/data/ld_seqs.fasta", 
              output_file = '/tmp/out_new.fasta',
              duplicates_file = '/tmp/out_new.csv',
              verbose = TRUE)

  opt <- list(input_file = "/tmp/outa.fasta", 
              output_file = '/tmp/reduped.fasta',
              duplicates_file = '/tmp/out.csv',
              reduplicate = TRUE,
              verbose = TRUE)
  
  opt <- list(input_file = "/tmp/out_miss.fasta", 
              output_file = '/tmp/reduped_miss.fasta',
              duplicates_file = '/tmp/out.csv',
              reduplicate = TRUE,
              verbose = TRUE)

}

opt <- parse_args(OptionParser(option_list = option_list,
  description = "Removes duplicates and appends _x to the sequence names where x is the number of sequences that were identical to the one that was left in the dataset. If --reduplicate is specified, then duplicates will be re-inserted into the file.",
  epilogue = "Example Call:
ElimDupes.R --input_file input.fasta --output_file output.fasta --duplicates_file counts_output.csv

File extensions must be specified in lower case."))

suppressPackageStartupMessages(library("hypermutR"))

if(!grepl('.fasta$', opt$input_file)){stop('input_file must end in .fasta (lowercase)')}
if(!grepl('.fasta$', opt$output_file)){stop('ouput_file must end in .fasta (lowercase)')}

if (opt$verbose){
  print("The parsed options:")
  print(opt)
}

###################################################################
# D:  DEDUPLICATION
###################################################################

deduplicate_operation <- function(opt){
  dat <- read.fasta(opt$input_file)
  dat <- sapply(dat, function(x){paste(x, collapse = '')})
  
  result <- deduplicate_seqs(dat)
  uniq_seqs <- sapply(result, function(x){x$the_seq})
  names(uniq_seqs) <- sapply(result, function(x){paste(x$dup_names[1], length(x$dup_names), sep = '_')})
  
  dups_counts <- table(sapply(result, function(x){length(x$dup_names)}))
  
  dups_tab <- data.frame(kept_seq = character(0),
                         dup_seq = character(0),
                         stringsAsFactors = FALSE)
  for (i in 1:length(result)){
    dups_tab <- rbind(dups_tab,
      data.frame(kept_seq = result[[i]]$dup_names[1],
                 dup_seq = result[[i]]$dup_names,
                 stringsAsFactors = FALSE))
  }
  stopifnot(nrow(dups_tab) == length(dat))
  return(list(dat = dat,
              uniq_seqs = uniq_seqs,
              dups_counts = dups_counts,
              dups_tab = dups_tab))
}

deduplicate_verbose_output <- function(dat, uniq_seqs, dups_counts){
  print(paste("Number of input sequences: ",  length(dat),       sep = ''))
  print(paste("Number of output sequences: ", length(uniq_seqs), sep = ''))

  for (i in 1:length(dups_counts)){
    num_dups <- as.numeric(names(dups_counts[i]))-1
    num_uniq_seqs <- dups_counts[i]
    first_bit <- ifelse(num_uniq_seqs == 1,
                   "There is only 1 unique sequence with",
                   paste("There are ", num_uniq_seqs, " unique sequences each with", sep = ''))
    last_bit <- ifelse(num_dups == 0,
                  "no duplicates",
                  ifelse(num_dups == 1,
                    "only one duplicate",
                    paste(num_dups, " duplicates.", sep = '')))
    print(paste(first_bit, last_bit, sep = ' '))
  }
  print(paste("Writing output deduplicated sequence file to ", opt$output_file, sep = ''))
  if (is.null(opt$duplicates_file)){
    print(paste("Writing the duplicates table to ", gsub(".fasta", ".csv", opt$output_file), sep = ''))
  } else {
    print(paste("Writing the duplicates table to ", opt$duplicates_file, sep = ''))
  }
}

write_deduplicate_outputs <- function(opt, uniq_seqs, dups_tab){
  write.fasta(sequences = as.list(uniq_seqs),
              names = names(uniq_seqs),
              file.out = opt$output_file)
  if (is.null(opt$duplicates_file)){
    opt$duplicates_file <- gsub(".fasta", ".csv", opt$output_file)
  }
  write.csv(dups_tab, opt$duplicates_file, row.names = FALSE)
}

deduplicate_wrapper <- function(opt){
  result <- deduplicate_operation(opt = opt)
  if (opt$verbose){
    deduplicate_verbose_output(dat = result$dat,
                               uniq_seqs = result$uniq_seqs,
                               dups_counts = result$dups_counts)
  }
  print(paste("Number of duplicates removed: ", 
              length(result$dat) - length(result$uniq_seqs), 
              sep = '')
  )
  write_deduplicate_outputs(opt = opt,
                            uniq_seqs = result$uniq_seqs,
                            dups_tab = result$dups_tab)
}

###################################################################
# R:  REDUPLICATION
###################################################################

reduplicate_operation <- function(opt){
  dat <- read.fasta(opt$input_file)
  dat <- sapply(dat, function(x){paste(x, collapse = '')})
  if (!is.null(opt$but_first_remove)){
    to_remove <- read.delim(opt$but_first_remove)
    dat <- dat[!(names(dat) %in% to_remove[,1])]
  }

  dups_tab <- read.csv(opt$duplicates_file, stringsAsFactors = FALSE)

  type1 <- sum(names(dat) %in% dups_tab$kept_seq)
  type2 <- sum(gsub("_[0-9]+$", "", names(dat)) %in% dups_tab$kept_seq)
  if (type2 > type1){
    #Sequence file in elim dupes and dups_table is not
    names(dat) <- gsub("_[0-9]+$", "", names(dat))
    if (opt$verbose){
      print("The Sequence file in ElimDupes format while the duplicates table is not. Stripping _[0-9]+$ from the names of the sequences in the sequence file.")
    }
  }

  not_found <- NULL
  found <- NULL

  reduped <- character(0)
  for (i in 1:nrow(dups_tab)){
    match_indx <- which(names(dat) == dups_tab[i,1])
    stopifnot(length(match_indx) <= 1)
    if (length(match_indx) == 0){
      not_found <- c(not_found, dups_tab[i,1])
    } else {
      reduped <- c(reduped, dat[match_indx])
      names(reduped)[length(reduped)] <- dups_tab[i,2]
      found <- c(found, paste(dups_tab[i,1], " for ", dups_tab[i,2], sep = ''))
    }
  }

  orphaned_inputs <- names(dat)[!(names(dat) %in% unique(dups_tab$kept_seq))]
  if (length(orphaned_inputs) > 0){
    x <- length(unique(orphaned_inputs))
    z <- length(input_dat)
    print(paste(x, " out of ", z, 
                " of the sequences in the input file were NOT found in the duplicates table.", 
                sep = ''))
    print("(This is a serious problem!)")
    print(orphaned_inputs)
  }

  return(
    list(input_dat = dat,
         reduped_seqs = reduped,
         not_found = not_found,
         found = found,
         orphaned_inputs = orphaned_inputs)
  )
}

reduplicate_write_output <- function(opt, reduped_seqs){
  write.fasta(sequences = as.list(reduped_seqs),
              names = names(reduped_seqs),
              file.out = opt$output_file)
}

reduplicate_verbose_output <- function(input_dat, not_found, found, orphaned_inputs){
  print("The sequences in the input file were matched against the duplicates table.")
  x <- length(unique(orphaned_inputs))
  z <- length(input_dat)
  print(paste(x, " out of ", z, 
              " of the sequences in the input file were NOT found in the duplicates table.", 
              sep = ''))
  x <- length(unique(not_found))
  z <- x + length(unique(found))
  print(paste(x, " out of ", z, 
              " of the sequences listed in the duplicates table were NOT found in the input sequences.", 
              sep = ''))

  if (length(not_found) > 0){
    print("Sequences listed in duplicates table not in the input data:")
    print("These were probably legitimately removed due to hypermutation etc.")
    not_found_tab <- table(not_found)
    for (i in 1:length(not_found_tab)){
      print(paste(names(not_found_tab)[i], " which has ", as.numeric(not_found_tab)[i], " duplicate(s).", sep = ''))
    }
  }
}

reduplicate_wrapper <- function(opt){
  result <- reduplicate_operation(opt = opt)
  if (opt$verbose){
    reduplicate_verbose_output(input_dat = result$input_dat,
                               not_found = result$not_found,
                               found = result$found,
                               orphaned_inputs = result$orphaned_inputs)
  }
  reduplicate_write_output(opt = opt,
                           reduped_seqs = result$reduped_seqs)
}

###############################################################
# Main "Loop"
###############################################################

if (opt$reduplicate){
  reduplicate_wrapper(opt = opt)
} else {
  deduplicate_wrapper(opt = opt)
}


# Point this script to a folder containing the output of Paul's hiv_founder_pipeline that removed hypermutation and recombination.

# It loops over the hiv_founder_id_nnnnn folders,
# inspects the identify_founders.tab
# extracting the last file name
# reads that file and concatenates the sequences into the all_sequences character vector.

library(seqinr)
library(stringr)
library(dplyr)

#SEQUENCES.DIR <- "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/";
#RESULTS.DIRNAME <- "raw_edited_20160216";

SEQUENCES.DIR <- "/fast/bakeoff_merged_analysis_sequences_results/results/";
#SEQUENCES.DIR <- "/fast/bakeoff_merged_analysis_sequences_results_2019/results/";
RESULTS.DIRNAME <- "raw_fixed";

THE.SEQUENCES.DIR <- SEQUENCES.DIR; # to avoid "promise already under evaluation" errors

prepareGeneSpecificDataSubset <- function (
  SEQUENCES.DIR = THE.SEQUENCES.DIR,
  results.dirname = RESULTS.DIRNAME,
  use.processed.lists = FALSE,
  force.recomputation = TRUE,
  genes = c( "GAG", "ENV" ),
  times = c( "1m", "6m", "1m6m" )
)
{

    ## For each results dir defined by time, call prepareGeneSpecificDataSubset.in.dir.
    for( the.time in times ) {
        prepareGeneSpecificDataSubset.in.dir(
            the.time,
            SEQUENCES.DIR = SEQUENCES.DIR,
            results.dirname = results.dirname,
            use.processed.lists = use.processed.lists,
            force.recomputation = force.recomputation,
            genes = genes
        );
    } # End foreach the.time
} # prepareGeneSpecificDataSubset (..)

prepareGeneSpecificDataSubset.in.dir <- function (
  the.time,
  SEQUENCES.DIR = THE.SEQUENCES.DIR,
  results.dirname = RESULTS.DIRNAME,
  use.processed.lists = FALSE,
  force.recomputation = TRUE,
  genes = c( "GAG", "ENV" )
)
{
    the.dir <-
        paste( SEQUENCES.DIR, results.dirname, "/nflg/", the.time, "/", sep = "" );
    output.dir <-
        paste( SEQUENCES.DIR, results.dirname, "/prepareGeneSpecificDataSubset/", sep = "" );
    dir.create( output.dir, showWarnings = FALSE );
    
    all.sequences.df.file <- paste(output.dir, the.time, "_all_sequences.csv", sep = '');
    all.sequences.fasta.file <- paste(output.dir, the.time, "_all_sequences.fasta", sep = '');

    if( force.recomputation || !file.exists( all.sequences.df.file ) || !file.exists( all.sequences.fasta.file ) ) {
        all.sequences.df <- data.frame(in_file = character(0),
                                  seq_name = character(0),
                                  n_non_gap = numeric(0),
                                  stringsAsFactors = FALSE);
        
        one.big.alignment.of.all.sequences.for.all.ppts <- NULL;
        
        # List files, one per participant, each contains the source fasta files for that participant, one per line.
        if( use.processed.lists ) {
            list.files <- dir( the.dir, pattern = "^processed_[0-9]+\\.list$" );
        } else {
            list.files <- dir( the.dir, pattern = "^[0-9]+\\.list$" );
        }
        names( list.files ) <- gsub( "^(processed_)?([0-9]+)\\.list$", "\\2", list.files );

        for( .ppt in names( list.files ) ) {
          #print("=========================")
          #print(.ppt)
          source.files.for.ppt <-
              read.delim( paste( the.dir, list.files[ .ppt ], sep = "" ), header = FALSE, sep = "\t", stringsAsFactors = FALSE )[, 1];
        
          for( source.file in source.files.for.ppt ) {
            #print(source.file)
            print( paste(source.file, file.exists( source.file ), sep = ' -> ') )
        
            source.file.contents <-
                read.fasta( source.file, seqtype = 'DNA', as.string = TRUE );
            seq.names <- sapply(source.file.contents, function(source.file.contents){attr(source.file.contents, "name")});
            n.non.gaps <- sapply(source.file.contents, function(source.file.contents){nchar(source.file.contents) - str_count(source.file.contents, "-")});
        
            all.sequences.df <-
                rbind(all.sequences.df,
                      data.frame(
                                 in_file = source.file,
                                 seq_name = seq.names,
                                 n_non_gap = n.non.gaps,
                                 stringsAsFactors = FALSE )
                      );
            one.big.alignment.of.all.sequences.for.all.ppts <-
                c( one.big.alignment.of.all.sequences.for.all.ppts,
                  source.file.contents );
          } # End foreach source.file for .ppt
        } # End foreach .ppt
        row.names( all.sequences.df ) <- NULL;
        write.csv( all.sequences.df, all.sequences.df.file, row.names = FALSE );

        all.sequences.alignment <-
            as.character( one.big.alignment.of.all.sequences.for.all.ppts );
        names( all.sequences.alignment ) <- sapply( one.big.alignment.of.all.sequences.for.all.ppts, function(x){attr(x, "name")} );
    
        write.fasta( as.list( all.sequences.alignment ), names( all.sequences.alignment ), all.sequences.fasta.file )
    } # End if we need to recompute the alignment and/or data frame of all source seqs.
    
    # load in the specific files and split them up again.
    
    which_region <- 'ENV'
    all.sequences.df <- read.csv( all.sequences.df.file, stringsAsFactors = FALSE );
    ## ERE I AM, DEBUGGING.
    region_seq_file <- read.fasta(paste(output.dir, 
                                        'from_lanl', 
                                        paste(which_region, '.NA.FASTA', sep = ''),
                                        sep = '/'),
                                  seqtype = 'DNA',
                                  as.string = TRUE)
    
    region_seq_file_vec <- as.character(region_seq_file)
    names(region_seq_file_vec) <- sapply(region_seq_file, function(x){attr(x, 'name')})
    
    region_resplit <- list()
    
    all.sequences.df$found <- 'no'
    for (i in 1:nrow(all.sequences.df)){
      
      if (all.sequences.df$seq_name[i] %in% names(region_seq_file_vec)){
        all.sequences.df$found[i] <- 'yes'
        the_seq <- region_seq_file_vec[all.sequences.df$seq_name[i] == names(region_seq_file_vec)]
        prop_gaps <- str_count(the_seq, "-") / nchar(the_seq)
        all.sequences.df$prop_gaps[i] <- prop_gaps
    
        # only save sequences with less than 90% gaps
        if (prop_gaps < 0.9){
          if (is.null(region_resplit[[all.sequences.df$in_file[i]]])){
            region_resplit[[all.sequences.df$in_file[i]]] <- the_seq
          } else {
            region_resplit[[all.sequences.df$in_file[i]]] <- c(region_resplit[[all.sequences.df$in_file[i] ]],
                                                          the_seq)
          }
        }
      }
    }
    
    gap_stats <- all.sequences.df %>% 
      group_by(in_file) %>%
      summarize(n = n(),
                avg_gaps = mean(prop_gaps),
                max_gaps = max(prop_gaps),
                min_gaps = min(prop_gaps),
                n_seqs_included = sum(prop_gaps < 0.9))
    
    dir.create( paste(output.dir, which_region, sep = '/') )
    for (source.file in names(region_resplit)){
      write.fasta(as.list(region_resplit[[source.file]]),
                  names(region_resplit[[source.file]]),
                  paste(output.dir, which_region, 
                        gsub(".fasta", 
                             paste("_", which_region, ".fasta", sep = ''), 
                             source.file), 
                        sep = '/'))
    }
    
    write.csv(gap_stats, 
              paste(output.dir, 
                    which_region,
                    paste('gap_stats_', which_region, '.csv', sep = ''),
                    sep = '/'),
              row.names = FALSE)

} # prepareGeneSpecificDataSubset.in.dir (..)

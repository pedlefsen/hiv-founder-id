# Point this script to a folder containing the output of Paul's hiv_founder_pipeline that removed hypermutation and recombination.

# It loops over the hiv_founder_id_nnnnn folders,
# inspects the identify_founders.tab
# extracting the last file name
# reads that file and concatenates the sequences into the all_sequences character vector.

library(seqinr)
library(stringr)
library(dplyr)

SEQUENCES.DIR <- "/fast/bakeoff_merged_analysis_sequences_unfiltered/results/";
#SEQUENCES.DIR <- "/fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/";
#SEQUENCES.DIR <- "/fast/bakeoff_merged_analysis_sequences_filtered2019/results/";
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

    ## Run the file at genecutter.
    gene.cutter.command <-
        paste( "perl -w ./runGeneCutterOnline.pl -pDV", all.sequences.fasta.file, output.dir );
    cat( paste( "Running", gene.cutter.command ), fill = TRUE );
    system( gene.cutter.command );

#    all.sequences.fasta.file.short <-
#        gsub( "^.*?\\/?([^\\/]+?)$", "\\1", all.sequences.fasta.file, perl = TRUE );
#    all.sequences.fasta.file.short.nosuffix <-
#        gsub( "^([^\\.]+)(\\..+)?$", "\\1", all.sequences.fasta.file.short, perl = TRUE );
#    all.sequences.fasta.file.short.nosuffix <-
#        paste( the.time, "_all_sequences", sep = "" );

    ## ERE I AM. TESTING THE FOLLOWING FOR EXPANDING THE FILES:
    
    ## First expand the output zip file.
    system( paste( "cd ", output.dir, "; mkdir ", the.time, "; cd ", the.time, "; unzip ../", the.time, "_all_sequences_allnucs.zip" ), sep = "" );
    
    ## Load in the specific files and split them up again.

    all.sequences.df <- read.csv( all.sequences.df.file, stringsAsFactors = FALSE );
    for( the.gene in genes ) {
        print( the.gene );
        
        the.gene.fasta <- read.fasta(paste(output.dir, 
                                        the.time, 
                                        paste(the.gene, '.NA.FASTA', sep = ''),
                                        sep = '/'),
                                  seqtype = 'DNA',
                                  as.string = TRUE);
    
        the.gene.seqs <- as.character(the.gene.fasta)
        names(the.gene.seqs) <- sapply(the.gene.fasta, function(x){attr(x, 'name')})
        
        the.gene.seqs.by.source.file <- list();
        all.sequences.df$found <- 'no'
        the.gene.seqs.names <- names(the.gene.seqs);
        for (i in 1:nrow(all.sequences.df)){
          
          if (all.sequences.df$seq_name[i] %in% the.gene.seqs.names){
            all.sequences.df$found[i] <- 'yes';
            the.seq <- the.gene.seqs[all.sequences.df$seq_name[i] == the.gene.seqs.names];
            gap.frequency <- str_count(the.seq, "-") / nchar(the.seq);
            all.sequences.df$gap.frequency[i] <- gap.frequency;
        
            # only save sequences with less than 90% gaps
            if (gap.frequency < 0.9){
              if (is.null(the.gene.seqs.by.source.file[[ all.sequences.df$in_file[i] ]])){
                the.gene.seqs.by.source.file[[ all.sequences.df$in_file[i] ]] <- the.seq;
              } else {
                  the.gene.seqs.by.source.file[[ all.sequences.df$in_file[i] ]] <-
                      c( the.gene.seqs.by.source.file[[ all.sequences.df$in_file[i] ]], the.seq );
              }
            }
          }
        } # End sorting the newly cut gene-specific sequences by their source file
        
        gap.stats <- all.sequences.df %>% 
          group_by(in_file) %>%
          summarize(n = n(),
                    avg_gaps = mean(gap.frequency),
                    max_gaps = max(gap.frequency),
                    min_gaps = min(gap.frequency),
                    n_seqs_included = sum(gap.frequency < 0.9))
        
        dir.create( paste(output.dir, the.gene, sep = '/') )
        for (source.file in names(the.gene.seqs.by.source.file)){
          write.fasta(as.list(the.gene.seqs.by.source.file[[source.file]]),
                      names(the.gene.seqs.by.source.file[[source.file]]),
                      paste(output.dir, the.gene, 
                            gsub(".fasta", 
                                 paste("_", the.gene, ".fasta", sep = ''), 
                                 source.file), 
                            sep = '/'))
        }
        
        write.csv(gap.stats, 
                  paste(output.dir, 
                        the.gene,
                        paste('gap_stats_', the.gene, '.csv', sep = ''),
                        sep = '/'),
                  row.names = FALSE)
    } # End foreach the.gene
    

} # prepareGeneSpecificDataSubset.in.dir (..)


library(seqinr)
library(stringr)
library(dplyr)

EXCLUSION.PATTERN = "env";

#SEQUENCES.DIR <- "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/";
#RESULTS.DIRNAME <- "raw_edited_20160216";

#SEQUENCES.DIR <- "/fast/bakeoff_merged_analysis_sequences_results/results/";
SEQUENCES.DIR <- "/fast/bakeoff_merged_analysis_sequences_results_2019/results/";
RESULTS.DIRNAME <- "raw_fixed";

THE.SEQUENCES.DIR <- SEQUENCES.DIR; # to avoid "promise already under evaluation" errors

cleanLists <- function (
  exclusion.pattern = EXCLUSION.PATTERN,
  SEQUENCES.DIR = THE.SEQUENCES.DIR,
  results.dirname = RESULTS.DIRNAME,
  use.processed.lists = TRUE,
  regions = c( "nflg", "v3", "rv217_v3" ),
  times = c( "1m", "6m", "1m6m" )
)
{

    ## For each results dir defined by region and time, call cleanLists.in.dir.
    for( the.region in regions ) {
        for( the.time in times ) {
            cleanLists.in.dir(
                the.region,
                the.time,
                exclusion.pattern = exclusion.pattern,
                SEQUENCES.DIR = SEQUENCES.DIR,
                results.dirname = results.dirname,
                use.processed.lists = use.processed.lists
            );
        } # End foreach the.time
    } # End foreach the.region
} # cleanLists (..)

cleanLists.in.dir <- function (
  the.region,
  the.time,
  exclusion.pattern = EXCLUSION.PATTERN,
  SEQUENCES.DIR = THE.SEQUENCES.DIR,
  results.dirname = RESULTS.DIRNAME,
  use.processed.lists = TRUE
)
{
    the.dir <-
        paste( SEQUENCES.DIR, results.dirname, "/", the.region, "/", the.time, "/", sep = "" );
    if( !file.exists( the.dir) ) {
        return;
    }
    
    print("-------------------------")
    print( the.region )
    print( the.time )
    # List files, one per participant, each contains the source fasta files for that participant, one per line.
    if( use.processed.lists ) {
        list.files <- dir( the.dir, pattern = "^processed_[0-9]+\\.list$" );
    } else {
        list.files <- dir( the.dir, pattern = "^[0-9]+\\.list$" );
        }
    names( list.files ) <- gsub( "^(processed_)?([0-9]+)\\.list$", "\\2", list.files );

    for( .ppt in names( list.files ) ) {
        print(.ppt)

        source.files.for.ppt <-
            read.delim( paste( the.dir, list.files[ .ppt ], sep = "" ), header = FALSE, sep = "\t", stringsAsFactors = FALSE )[, 1];

        filtered.source.files <- c();
        for( source.file in source.files.for.ppt ) {
            #print(source.file)
            if( length( grep( exclusion.pattern, source.file ) ) > 0 ) {
                cat( paste( "Excluding", source.file, "because it matches the exclusion pattern." ), fill = TRUE );
                next;
            }
            if( !file.exists( source.file ) ) {
                cat( paste( "Excluding", source.file, "because FILE NOT FOUND." ), fill = TRUE );
                next;
            }
            filtered.source.files <- c( filtered.source.files, source.file );
        } # End foreach source.file

        if( length( filtered.source.files ) < length( source.files.for.ppt ) ) {
            stopifnot( length( filtered.source.files ) > 0 );
            names( filtered.source.files ) <- NULL;
            write.table( filtered.source.files, paste( the.dir, list.files[ .ppt ], sep = "" ), quote = FALSE, col.names = FALSE, row.names = FALSE );
        }
    } # End foreach .ppt
} # cleanLists.in.dir (..)

## Here is where the action is.
cleanLists()

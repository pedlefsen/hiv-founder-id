library( "ROCR" ) # for "prediction" and "performance"

GOLD.STANDARD.DIR <- "/fh/fast/edlefsen_p/bakeoff/gold_standard";
#RESULTS.DIR <- "/fh/fast/edlefsen_p/bakeoff_merged_analysis_sequences/raw_fixed";
RESULTS.DIR <- "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/";
THE.TIMES <- c( "1m", "6m", "1m6m" );

# Here "the.study" must be "rv217" or "caprisa002" or (special) "rv217_v3"
evaluateIsMultiple <- function ( the.study, output.dir = NULL, output.file = NULL, output.file.append = FALSE ) {

    ## Read in the gold standards.
    # From this file we read in indicators of whether to use the multiple- or single-founder true profile.
    if( ( the.study == "rv217" ) || ( the.study == "rv217_v3" ) ) {
      gold.standards.in <- read.csv( paste( GOLD.STANDARD.DIR, "/rv217/RV217_gold_standards.csv", sep = "" ) );
    } else if( the.study == "caprisa002" ) {
      # From this file we read in indicators of whether to use the multiple- or single-founder true profile.
      gold.standards.in <- read.csv( paste( GOLD.STANDARD.DIR, "/caprisa_002/caprisa_002_gold_standards.csv", sep = "" ) );
    } else {
        stop( paste( "unrecognized value of the.study:", the.study ) );
    }
    gold.is.multiple <- gold.standards.in[ , "gold.is.multiple" ];
    names( gold.is.multiple ) <- gold.standards.in[ , "ptid" ];

    if( !is.null( output.file ) ) {
      if( length( grep( "^(.*?)\\/[^\\/]+$", output.file ) ) == 0 ) {
          output.file.path <- NULL;
          output.file.path.is.absolute <- NA;
      } else {
          output.file.path <-
              gsub( "^(.*?)\\/[^\\/]+$", "\\1", output.file );
          output.file.path.is.absolute <- ( substring( output.file.path, 1, 1 ) == "/" );
      }
      output.file.short <-
          gsub( "^.*?\\/?([^\\/]+?)$", "\\1", output.file, perl = TRUE );

      if( !is.null( output.file.path ) && output.file.path.is.absolute ) {
          output.dir <- output.file.path;
      } else if( is.null( output.dir ) ) {
        if( is.null( output.file.path ) ) {
            output.dir <- ".";
        } else {
            output.dir <- output.file.path;
        }
      } else {
          output.dir <- paste( output.dir, output.file.path, sep = "/" );
      }
      output.file <- output.file.short;
    } else { # is.null( output.file )
      output.file <- paste( the.study, "_evaluateIsMultiple.tab", sep = "" );
    }
    if( is.null( output.dir ) || ( nchar( output.dir ) == 0 ) ) {
      output.dir <- ".";
    }
    ## Remove "/" from end of output.dir
    output.dir <-
      gsub( "^(.*?)\\/+$", "\\1", output.dir );

    evaluateIsMultiple.OneFile <- function ( identify.founders.tab.file ) {
        cat( identify.founders.tab.file, fill = T );
        
        if( length( grep( "^(.*?)\\/[^\\/]+$", identify.founders.tab.file ) ) == 0 ) {
            identify.founders.tab.file.path <- ".";
        } else {
            identify.founders.tab.file.path <-
                gsub( "^(.*?)\\/[^\\/]+$", "\\1", identify.founders.tab.file );
        }
        identify.founders.tab.file.short <-
            gsub( "^.*?\\/?([^\\/]+?)$", "\\1", identify.founders.tab.file, perl = TRUE );
        identify.founders.tab.file.short.nosuffix <-
            gsub( "^([^\\.]+)(\\..+)?$", "\\1", identify.founders.tab.file.short, perl = TRUE );
        identify.founders.tab.file.suffix <-
            gsub( "^([^\\.]+)(\\..+)?$", "\\2", identify.founders.tab.file.short, perl = TRUE );
    
        identify.founders.in <- NULL;
        try.count <- 0;
        TOO.MANY.TRIES <- 100;
        while( ( try.count < TOO.MANY.TRIES ) && ( is.null( identify.founders.in ) || is.null( nrow( identify.founders.in ) ) || ( nrow( identify.founders.in ) == 0 ) ) ) {
           if( try.count > 0 ) {
               warning( paste( "Trouble reading results in file", identify.founders.tab.file, ".. trying again in a second." ) );
               Sys.sleep( 1 );
           }
#           identify.founders.in <- read.delim( identify.founders.tab.file, sep = "\t", header = TRUE, fill = FALSE );

           ### FILL = TRUE is needed for now because of a bug in the multi-region lines in the identify_founders.tab file, which is too short.  This is ok because we are ignoring those results for now anyway.
           identify.founders.in <- read.delim( identify.founders.tab.file, sep = "\t", header = TRUE, fill = TRUE );
        }
        if( try.count >= TOO.MANY.TRIES ) {
            stop( paste( "COULD NOT READ results in file", identify.founders.tab.file, ".. GIVING UP." ) );
        }
    
#         ## For now we only look at the caprisa002 results in "v3"
#         identify.founders.study <-
#             identify.founders.in[ grep( paste( ".*", the.study, "_([^_]+)_.+", sep = "" ), as.character( identify.founders.in[ , "infile" ] ) ), , drop = FALSE ];
        # Fix accidental temporary bug in which second round of results isn't called "multifounder" (in R it gets called .1 instead).
        identify.founders.study <- identify.founders.in;
        
        colnames( identify.founders.study ) <- 
            gsub( "^(.+)\\.1$", "multifounder.\\1", colnames( identify.founders.study ) )
        is.one.founder.methods <-
            grep( "is\\.|fits", colnames( identify.founders.study ), value = T );

        ## Fix other accidental temporary bug in which the multiregion results are all messed up.  Just remove them.
        identify.founders.study <-
          identify.founders.study[ !is.na( identify.founders.study[ , "file" ] ), , drop = FALSE ];

        if( the.study == "rv217_v3" ) {
            identify.founders.ptids <-
                gsub( paste( ".*", "rv217", "_([^_]+)_.+", sep = "" ), "\\1", as.character( identify.founders.study[ , "infile" ] ) );
        } else {
            identify.founders.ptids <-
                gsub( paste( ".*", the.study, "_([^_]+)_.+", sep = "" ), "\\1", as.character( identify.founders.study[ , "infile" ] ) );
        }
        
        ## _All_ of these are "is.one.founder", so note that in the comparison they need to be reversed.
        estimates.is.one.founder <-
            identify.founders.study[ , is.one.founder.methods, drop = FALSE ];
        ## Sometimes there are multiple entries for one ptid/sample, eg for the NFLGs there are often right half and left half (RH and LH) and sometimes additionally NFLG results.  If so, for something to be called single founder, all of the estimates (across regions, for a given test) must agree that it is one founder.

        ## Note that this also transposes it, so the ptids are in columns and the tests are in rows.
        estimates.is.one.founder.one.per.person <- 
          sapply( unique( identify.founders.ptids ), function ( .ptid ) {
            .ptid.subtable <- estimates.is.one.founder[ identify.founders.ptids == .ptid, , drop = FALSE ];
            # Use the AND condition, meaning it's only called single-founder if all rows call it single-founder.
            if( nrow( .ptid.subtable ) == 1 ) {
              return( .ptid.subtable );
            }
            apply( .ptid.subtable, 2, function( .column ) { as.numeric( all( as.logical( .column ) ) ) } );
          } );

        gold.is.one.founder.per.person <- 1-gold.is.multiple[ colnames( estimates.is.one.founder.one.per.person ) ];
        
        # Cols are gold.is.multiple
        # gold.is.multiple.tables <- 
        #     sapply( 1:length( is.one.founder.methods ), function( .i ) { table( 1 - estimates.is.one.founder[ , .i ], gold.is.multiple[ identify.founders.ptids ] ) } );
        # names( gold.is.multiple.tables ) <- is.one.founder.methods;

        mode( estimates.is.one.founder.one.per.person ) <- "numeric";
        sum.correct.among.one.founder.people <-
          apply( estimates.is.one.founder.one.per.person[ , as.logical( gold.is.one.founder.per.person ) ], 1, sum );
        sum.incorrect.among.one.founder.people <-
          apply( estimates.is.one.founder.one.per.person[ , as.logical( gold.is.one.founder.per.person ) ], 1, function( .row ) { sum( 1-.row ) } );
        sum.correct.among.multiple.founder.people <-
          apply( estimates.is.one.founder.one.per.person[ , !as.logical( gold.is.one.founder.per.person ) ], 1, function( .row ) { sum( 1-.row ) } );
        sum.incorrect.among.multiple.founder.people <-
          apply( estimates.is.one.founder.one.per.person[ , !as.logical( gold.is.one.founder.per.person ) ], 1, sum );
        
        gold.is.multiple.aucs <- 
            sapply( 1:length( is.one.founder.methods ), function( .i ) {
                performance( prediction( unlist( estimates.is.one.founder.one.per.person[ .i, ] ), gold.is.one.founder.per.person ), measure = "auc" )@y.values[[ 1 ]];
            } );
        names( gold.is.multiple.aucs ) <- is.one.founder.methods;

        return( gold.is.multiple.aucs );
    } # evaluateIsMultiple.OneFile ( identify.founders.tab.file )

    if( the.study == "rv217" ) {
        the.region <- "nflg";
    } else if( the.study == "rv217_v3" ) {
        the.region <- "rv217_v3";
    } else if( the.study == "caprisa002" ) {
        the.region <- "v3";
    } else {
        stop( paste( "again?!? unrecognized value of the.study:", the.study ) );
    }
        
    results.by.time <- lapply( THE.TIMES, function( the.time ) {
        cat( the.time, fill = TRUE );
        identify.founders.tab.file <- paste( RESULTS.DIR, the.region, the.time, "identify_founders.tab", sep = "/" );
        stopifnot( file.exists( identify.founders.tab.file ) );
        return( evaluateIsMultiple.OneFile( identify.founders.tab.file ) );
    } );
    names( results.by.time ) <- THE.TIMES;
    
    results.matrix <-
        do.call( cbind, results.by.time );
    results.matrix <-
        do.call( cbind, results.by.time );
    # No need for so much precision.  2 digits should suffice.
    results.matrix <- apply( results.matrix, 1:2, function ( .float ) { sprintf( "%0.2f", .float ) } ); 
    output.table.path <-
        paste( output.dir, "/", output.file, sep = "" );

    write.table( results.matrix, file = output.table.path, append = ( file.exists( output.table.path ) && output.file.append ), row.names = TRUE, col.names = ( !output.file.append || !file.exists( output.table.path ) ), sep = "\t", quote = FALSE );

    # Return the file name.
    return( output.table.path );
} # evaluateIsMultiple ( the.study, ... )

## Here is where the action is.
the.study <- Sys.getenv( "evaluateIsMultiple_study" );
output.table.file <- Sys.getenv( "evaluateIsMultiple_outputFilename" );
if( nchar( output.table.file ) == 0 ) {
    output.table.file <- NULL;
}
output.dir <- Sys.getenv( "evaluateIsMultiple_outputDir" );
if( nchar( output.dir ) == 0 ) {
    output.dir <- NULL;
}
append.to.output.file <- Sys.getenv( "evaluateIsMultiple_append" );
if( ( nchar( append.to.output.file ) == 0 ) || ( append.to.output.file == "0" ) || ( toupper( append.to.output.file ) == "F" ) || ( toupper( append.to.output.file ) == "FALSE" ) ) {
    append.to.output.file <- FALSE;
} else {
    append.to.output.file <- TRUE;
}

print( evaluateIsMultiple( the.study, output.dir = output.dir, output.file = output.table.file, output.file.append = append.to.output.file ) );

library( "ROCR" ) # for "prediction" and "performance"

source( "readIdentifyFounders_safetosource.R" );
source( "summarizeCovariatesOnePerParticipant_safetosource.R" );

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
    
             ## identify-founders results
             identify.founders.study <- readIdentifyFounders( identify.founders.tab.file );
             single.colnames <- grep( "\\.is\\.|fits", colnames( identify.founders.study ), perl = TRUE, value = TRUE );

        ## _All_ of these are "is.one.founder", so note that in the comparison they need to be reversed.
        estimates.is.one.founder <-
            identify.founders.study[ , single.colnames, drop = FALSE ];
        identify.founders.ptids <- rownames( identify.founders.study );
        ## Sometimes there are multiple entries for one ptid/sample, eg for the NFLGs there are often right half and left half (RH and LH) and sometimes additionally NFLG results.  If so, for something to be called single founder, all of the estimates (across regions, for a given test) must agree that it is one founder.

        estimates.is.one.founder.one.per.person <- 
          t( sapply( unique( identify.founders.ptids ), function ( .ptid ) {
            .ptid.subtable <- estimates.is.one.founder[ identify.founders.ptids == .ptid, , drop = FALSE ];
            # Use the AND condition, meaning it's only called single-founder if all rows call it single-founder.
            if( nrow( .ptid.subtable ) == 1 ) {
              return( .ptid.subtable );
            }
            apply( .ptid.subtable, 2, function( .column ) { as.numeric( all( as.logical( .column ) ) ) } );
          } ) );
        colnames( estimates.is.one.founder.one.per.person ) <- colnames( estimates.is.one.founder );

        gold.is.one.founder.per.person <-
            1-gold.is.multiple[ rownames( estimates.is.one.founder.one.per.person ) ];

        ## ERE I AM.  I'm going to add some estimates made by prediction using leave-one-out cross-validatino.
        results.covars.one.per.ppt.with.extra.cols <-
            summarizeCovariatesOnePerParticipant( identify.founders.study );
        
        .keep.cols <-
            grep( "num.*\\.seqs|totalbases", colnames( results.covars.one.per.ppt.with.extra.cols ), value = TRUE, perl = TRUE, invert = TRUE );
        ### TODO: Something else.  Just trying to get down to a reasonable set; basically there are very highly clustered covariates here and it screws up the inference.
        ## Also remove all of the mut.rate.coef except for multifounder.Synonymous.PFitter.mut.rate.coef.
        
        #.keep.cols <- c( "multifounder.Synonymous.PFitter.mut.rate.coef",
        #                 grep( "mut\\.rate\\.coef", .keep.cols, invert = TRUE, value = TRUE ) )
        #.keep.cols <- c( "multifounder.Synonymous.PFitter.mut.rate.coef", "inf.to.priv.ratio", "priv.sites", "inf.sites.clusters", "InSites.founders", "multifounder.Synonymous.PFitter.is.poisson" );
        if( the.region == "v3" ) {
            #helpful.additional.cols <- c( "priv.sites" )
            helpful.additional.cols <- c();
        } else {
            helpful.additional.cols <- c();
        }
        mut.rate.cols <- grep( "mut\\.rate\\.coef", .keep.cols, value = TRUE );
        .keep.cols <- c( helpful.additional.cols, mut.rate.cols );
        
        results.covars.one.per.ppt <-
            results.covars.one.per.ppt.with.extra.cols[ , .keep.cols, drop = FALSE ];
        results.covars.one.per.ppt.df <-
            data.frame( results.covars.one.per.ppt );
        
        regression.df <-
            cbind( data.frame( is.one.founder = gold.is.one.founder.per.person[ rownames( results.covars.one.per.ppt.df ) ] ), results.covars.one.per.ppt.df );
        
        # logistic.fit.formula <- as.formula( paste( "is.one.founder ~ ", paste( c( helpful.additional.cols, mut.rate.cols[ 1 ] ), collapse = "+" ) ) );
        # summary( logistic.fit <- glm( logistic.fit.formula, family = "binomial", data = regression.df ) );

        ## new proof of concept:
        helpful.additional.parameters.validation.results.one.per.ppt <- matrix( NA, nrow = nrow( results.covars.one.per.ppt.df ), ncol = length( mut.rate.cols ) );
        for( .row.i in 1:nrow( regression.df ) ) {
            for( .col.i in 1:length( mut.rate.cols ) ) {
                .mut.rate.coef.colname <- mut.rate.cols[ .col.i ];
                ## Ok build a regression model with no intercept, including only the helpful.additional.cols
                if( length( helpful.additional.cols ) == 0 ) {
                    .formula <- as.formula( paste( "is.one.founder ~", .mut.rate.coef.colname ) );
                } else {
                    .formula <- as.formula( paste( "is.one.founder ~", paste( c( helpful.additional.cols, .mut.rate.coef.colname ), collapse = "+" ) ) );
                }
                .pred.value <-
                    predict( glm( .formula, family = "binomial", data = regression.df[ -.row.i, ] ), regression.df[ .row.i, , drop = FALSE ] );
                helpful.additional.parameters.validation.results.one.per.ppt[ .row.i, .col.i ] <- 
                    .pred.value;
            } # End foreach .col.i
        } # End foreach .row.i
        colnames( helpful.additional.parameters.validation.results.one.per.ppt ) <-
            paste( "glm.validation", mut.rate.cols, sep = "." );
        rownames( helpful.additional.parameters.validation.results.one.per.ppt ) <-
            rownames( regression.df );
    
        estimates.is.one.founder.one.per.person <-
            cbind( estimates.is.one.founder.one.per.person,
                  helpful.additional.parameters.validation.results.one.per.ppt );
        
        mode( estimates.is.one.founder.one.per.person ) <- "numeric";
        
        # sum.correct.among.one.founder.people <-
        #   apply( estimates.is.one.founder.one.per.person[ as.logical( gold.is.one.founder.per.person ),  ], 2, sum );
        # sum.incorrect.among.one.founder.people <-
        #   apply( estimates.is.one.founder.one.per.person[ as.logical( gold.is.one.founder.per.person ),  ], 2, function( .row ) { sum( 1-.row ) } );
        # sum.correct.among.multiple.founder.people <-
        #   apply( estimates.is.one.founder.one.per.person[ !as.logical( gold.is.one.founder.per.person ),  ], 2, function( .row ) { sum( 1-.row ) } );
        # sum.incorrect.among.multiple.founder.people <-
        #   apply( estimates.is.one.founder.one.per.person[ !as.logical( gold.is.one.founder.per.person ),  ], 2, sum );
        
        isMultiple.aucs <- 
            sapply( 1:ncol( estimates.is.one.founder.one.per.person ), function( .col.i ) {
                #print( .col.i );
                performance( prediction( as.numeric( estimates.is.one.founder.one.per.person[ , .col.i ] ), gold.is.one.founder.per.person ), measure = "auc" )@y.values[[ 1 ]];
            } );
        names( isMultiple.aucs ) <- colnames( estimates.is.one.founder.one.per.person );

        return( isMultiple.aucs );
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

## First do all the stuff in README.postprocessing.txt.
## [TODO: FUTURE: Also run evaluateTimings.R because its outputs are used as inputs here]
library( "ROCR" ) # for "prediction" and "performance"
library( "parallel" ) # for mcapply

source( "readIdentifyFounders_safetosource.R" );
source( "getArtificialBounds_safetosource.R" );
source( "summarizeCovariatesOnePerParticipant_safetosource.R" );

GOLD.STANDARD.DIR <- "/fh/fast/edlefsen_p/bakeoff/gold_standard";
RESULTS.DIR <- "/fh/fast/edlefsen_p/bakeoff_analysis_results/";

regions <- c( "nflg", "v3", "rv217_v3" );
times <- c( "1m", "6m", "1m6m" );

## MARK

#' Evaluate isMultiple estimates and produce results tables.
#'
#' This function runs the BakeOff results analysis for the isMultiple results.
#'
#' @param use.glm.validate evaluate predicted values from leave-one-out cross-validation.
#' @param include.bounds.in.glm include the corresponding prior bounds in the regression equations used for cross-validation.
#' @param include.bounds.in.lasso include the corresponding prior bounds in the regression equations used for cross-validation.
#' @param include.helpful.additional.cols.in.glm include helpful.additional.cols in the glm (but don't force them to stay in the lasso)? By default this is the opposite of include.bounds.in.glm).
#' @param helpful.additional.cols [TODO: extra cols to be included in the glm. By default, "pred.days", which is a special code meaning that the results should be computed for _each_ of the sets of preditions in the set: { same-time-and-region, same-time, same-region, all-times-and-regions }.]
#' @param results.dirname the subdirectory of "/fh/fast/edlefsen_p/bakeoff/analysis_sequences" and also of "/fh/fast/edlefsen_p/bakeoff_analysis_results"
#' @param force.recomputation if FALSE (default) and if there is a saved version called isMultiple.results.by.region.and.time.Rda (under bakeoff_analysis_results/results.dirname), then that file will be loaded; otherwise the results will be recomputed and saved in that location.
#' @param partition.bootstrap.seed the random seed to use when bootstrapping samples by selecting one partition number per ptid, repeatedly; we do it this way because there are an unequal number of partitions, depending on sampling depth.
#' @param partition.bootstrap.samples the number of bootstrap replicates to conduct; the idea is to get an estimate of the variation in estimates and results (errors) across these samples.
#' @param partition.bootstrap.num.cores the number of cores to run the boostrap replicates on (defaults to all of the cores returned by parallel::detectCores()).
#' @return NULL
#' @export

evaluateIsMultiple <- function (
                             use.bounds = TRUE,
                             use.glm.validate = TRUE,
                             use.lasso.validate = TRUE,
                             include.bounds.in.glm = TRUE,
                             include.bounds.in.lasso = TRUE,
                             include.helpful.additional.cols.in.glm = !include.bounds.in.glm,
                             helpful.additional.cols = c(), # c( "pred.days" )
                             results.dirname = "raw_edited_20160216",
                             force.recomputation = FALSE,
                             partition.bootstrap.seed = 98103,
                             partition.bootstrap.samples = 100,
                             partition.bootstrap.num.cores = detectCores(),
                             regions = c( "nflg", "v3", "rv217_v3" ),
                             times = c( "1m", "6m", "1m6m" )
                            )
{
    isMultiple.results.by.region.and.time.Rda.filename <-
        paste( RESULTS.DIR, results.dirname, "/isMultiple.results.by.region.and.time.Rda", sep = "" );

    ## Read in the gold standards.
    # From this file we read in indicators of whether to use the multiple- or single-founder true profile.
    rv217.gold.standards.in <- read.csv( paste( GOLD.STANDARD.DIR, "/rv217/RV217_gold_standards.csv", sep = "" ) );
    rv217.gold.is.multiple <- rv217.gold.standards.in[ , "gold.is.multiple" ];
    names( rv217.gold.is.multiple ) <- rv217.gold.standards.in[ , "ptid" ];
    
    caprisa002.gold.standards.in <- read.csv( paste( GOLD.STANDARD.DIR, "/caprisa_002/caprisa_002_gold_standards.csv", sep = "" ) );
    caprisa002.gold.is.multiple <- caprisa002.gold.standards.in[ , "gold.is.multiple" ];
    names( caprisa002.gold.is.multiple ) <- caprisa002.gold.standards.in[ , "ptid" ];

    ## Sometimes there are multiple entries for one ptid/sample, eg
    ## for the NFLGs there are often right half and left half (RH and
    ## LH) and sometimes additionally NFLG results.  If so, for
    ## something to be called single founder, all of the estimates
    ## (across regions, for a given test) must agree that it is one
    ## founder.
    compute.results.per.person <- function ( results ) {
        identify.founders.ptids <- rownames( results );
            
        results.per.person <- 
            t( sapply( unique( identify.founders.ptids ), function ( .ptid ) {
            .ptid.subtable <- results[ identify.founders.ptids == .ptid, , drop = FALSE ];
            # Use the AND condition, meaning it's only called single-founder if all rows call it single-founder.
            if( nrow( .ptid.subtable ) == 1 ) {
              return( .ptid.subtable );
            }
            apply( .ptid.subtable, 2, function( .column ) { as.numeric( all( as.logical( .column ) ) ) } );
        } ) );
            
        colnames( results.per.person ) <- colnames( results );
        return( results.per.person );
    } # compute.results.per.person (..)

    evaluate.results.per.ppt <-
        function ( estimates.is.one.founder.per.person, gold.is.one.founder.per.person, results.covars.per.person.with.extra.cols, the.time, the.artificial.bounds = NA, ppt.suffix.pattern = "\\.([^\\.]+)?\\.[16]m(6m)?$" ) {
       ## Special: the ppt names might have suffices in estimates.is.one.founder.per.person; if so, strip off the suffix for purposes of matching ppts to the covars, etc.
       ppt.names <- rownames( estimates.is.one.founder.per.person );
       if( !is.na( ppt.suffix.pattern ) ) {
           .ppt.names <- ppt.names;
           ppt.names <- gsub( ppt.suffix.pattern, "", .ppt.names );
           names( ppt.names ) <- .ppt.names;
       }
       ## TODO: Use this for something.

       if( use.glm.validate || use.lasso.validate ) {
           results.covars.per.person.with.extra.cols <-
               cbind( estimates.is.one.founder.per.person, results.covars.per.person.with.extra.cols );
        
           .keep.cols <-
               grep( "num.*\\.seqs|totalbases", colnames( results.covars.per.person.with.extra.cols ), value = TRUE, perl = TRUE, invert = TRUE );
           single.cols <- grep( "\\.is\\.|fits", colnames( identify.founders.study ), perl = TRUE, value = TRUE );
           mut.rate.coef.cols <- grep( "mut\\.rate\\.coef", .keep.cols, value = TRUE );
           all.additional.cols <- setdiff( .keep.cols, c( single.cols, mut.rate.coef.cols ) );
          
        if( use.lasso.validate ) {
            keep.cols <- c( all.additional.cols, single.cols, mut.rate.coef.cols );
            estimate.cols <- setdiff( keep.cols, all.additional.cols );
        } else {
            keep.cols <- c( helpful.additional.cols, single.cols, mut.rate.coef.cols );
            estimate.cols <- setdiff( keep.cols, helpful.additional.cols );
        }
          
        results.covars.per.person <-
            results.covars.per.person.with.extra.cols[ , keep.cols, drop = FALSE ];
        results.covars.per.person.df <-
            data.frame( results.covars.per.person );

        regression.df <- cbind( data.frame( is.one.founder = gold.is.one.founder.per.person[ rownames( results.covars.per.person.df ) ] ), results.covars.per.person.df, lapply( the.artificial.bounds, function( .mat ) { .mat[ rownames( results.covars.per.person.df ), , drop = FALSE ] } ) );
        
        ## Ok build a regression model with no intercept, including only the helpful.additional.cols, and also the lower and upper bounds associated with either 5 weeks or 30 weeks, depending on the.time (if there's a 1m sample, uses "5weeks").
        if( the.time == "6m" ) {
            .lower.bound.colname <- "uniform_30weeks.lower";
            .upper.bound.colname <- "uniform_30weeks.upper";
        } else {
            .lower.bound.colname <- "uniform_5weeks.lower";
            .upper.bound.colname <- "uniform_5weeks.upper";
        }
        
        if( use.glm.validate ) {
            glm.validation.estimates.is.one.founder.per.person <- matrix( NA, nrow = nrow( results.covars.per.person.df ), ncol = length( estimate.cols ) );
        }
        if( use.lasso.validate ) {
            lasso.validation.estimates.is.one.founder.per.person <- matrix( NA, nrow = nrow( results.covars.per.person.df ), ncol = length( estimate.cols ) );
        }
        for( .row.i in 1:nrow( regression.df ) ) {
            regression.df.without.row.i <-
                regression.df[ -.row.i, , drop = FALSE ];
            for( .col.i in 1:length( estimate.cols ) ) {
                .estimate.colname <- estimate.cols[ .col.i ];
                if( use.glm.validate ) {
                  # covariates for glm
                  .covariates.glm <- c();
                  if( include.helpful.additional.cols.in.glm ) {
                      .covariates.glm <-
                          c( .covariates.glm, helpful.additional.cols );
                  }
                  if( include.bounds.in.glm ) {
                      ## This was causing a problem because of perfect collinearity, so I've taken the upper bound away.
                      .covariates.glm <-
                          c( .covariates.glm, .lower.bound.colname );
                      #     c( .covariates.glm, .lower.bound.colname, .upper.bound.colname );
                  }
                  # glm:
                  .covars.to.exclude <- apply( regression.df.without.row.i, 2, function ( .col ) {
                      return( ( var( .col ) == 0 ) || ( sum( !is.na( .col ) ) <= 1 ) );
                  } );
                  .retained.covars <- setdiff( colnames( regression.df.without.row.i ), names( which( .covars.to.exclude ) ) );
                  if( .estimate.colname %in% .retained.covars ) {
                    .df <- regression.df.without.row.i[ , .retained.covars, drop = FALSE ];
                    .formula <- as.formula( paste( "is.one.founder ~ ", paste( intersect( .retained.covars, c( .covariates.glm, .estimate.colname ) ), collapse = "+" ) ) );
                    .pred.value.glm <- predict( glm( .formula, family = "binomial", data = regression.df.without.row.i ), regression.df[ .row.i, , drop = FALSE ] );
                    glm.validation.estimates.is.one.founder.per.person[ .row.i, .col.i ] <- 
                        .pred.value.glm;
                  } # If the estimate can be used
    
                } # End if use.glm.validate
    
                if( use.lasso.validate ) {
                  # covariates for lasso
                  .covariates.lasso <- c( all.additional.cols );
                  if( include.bounds.in.lasso ) {
                      ## This was causing a problem because of perfect collinearity, so I've taken the upper bound away.
                      .covariates.lasso <-
                          c( .covariates.lasso, .lower.bound.colname );
                      #     c( .covariates.lasso, .lower.bound.colname, .upper.bound.colname );
                  }
                  # lasso:
                  .mat1 <- as.matrix( regression.df.without.row.i[ , c( .covariates.lasso, .estimate.colname ) ] );
                  .out <- regression.df.without.row.i[[ "is.one.founder" ]];
    
                  .covars.to.exclude <- apply( .mat1, 2, function ( .col ) {
                      return( ( var( .col ) == 0 ) || ( sum( !is.na( .col ) ) <= 1 ) );
                  } );
                  .retained.covars <- setdiff( colnames( .mat1 ), names( which( .covars.to.exclude ) ) );
                  if( .estimate.colname %in% .retained.covars ) {
                    .mat1 <- .mat1[ , setdiff( colnames( .mat1 ), names( which( .covars.to.exclude ) ) ), drop = FALSE ];
                    # penalty.factor = 0 to force the .estimate.colname variable.
    
                    tryCatch( {
                    cv.glmnet.fit <- cv.glmnet( .mat1, .out, family = "binomial",
                                               penalty.factor = as.numeric( colnames( .mat1 ) != .estimate.colname ) );
                    .pred.value.lasso <- predict( cv.glmnet.fit, newx = as.matrix( regression.df[ .row.i, colnames( .mat1 ), drop = FALSE ] ), s = "lambda.min" );
                            lasso.validation.estimates.is.one.founder.per.person[ .row.i, .col.i ] <<- 
                        .pred.value.lasso;
                     },
                     error = function( e ) {
                         warning( paste( "lasso failed with error", e, "\nReverting to simple regression vs", .estimate.colname ) );
                         .formula <- as.formula( paste( "is.one.founder ~", .estimate.colname ) );
                         .pred.value.lasso <- predict( glm( .formula, family = "binomial", data = regression.df.without.row.i ), regression.df[ .row.i, , drop = FALSE ] );
                            lasso.validation.estimates.is.one.founder.per.person[ .row.i, .col.i ] <<- 
                        .pred.value.lasso;
                     },
                        finally = {
                        }
                        );
                  } # End if the estimate variable is usable
                } # End if use.lasso.validate
    
            } # End foreach .col.i
        } # End foreach .row.i
        if( use.glm.validate ) {
            colnames( glm.validation.estimates.is.one.founder.per.person ) <-
                paste( "glm.validation.results", estimate.cols, sep = "." );
            rownames( glm.validation.estimates.is.one.founder.per.person ) <-
                rownames( regression.df );
            estimates.is.one.founder.per.person <-
                cbind( estimates.is.one.founder.per.person,
                      glm.validation.estimates.is.one.founder.per.person );
        }
        if( use.lasso.validate ) {
            colnames( lasso.validation.estimates.is.one.founder.per.person ) <-
                paste( "lasso.validation.results", estimate.cols, sep = "." );
            rownames( lasso.validation.estimates.is.one.founder.per.person ) <-
                rownames( regression.df );
            estimates.is.one.founder.per.person <-
                cbind( estimates.is.one.founder.per.person,
                      lasso.validation.estimates.is.one.founder.per.person );
        }
      } # End if use.glm.validate || use.lasso.validate
      
       mode( estimates.is.one.founder.per.person ) <- "numeric";
       
      ## For fairness in evaluating when some methods
      ## completely fail to give a result, (so it's NA
      ## presently), we change all of these estimates
      ## from NA to 0.  When we put bounds on the
      ## results, below, they will be changed from 0 to
      ## a boundary endpoint if 0 is outside of the
      ## bounds.
      estimates.is.one.founder.per.person.zeroNAs <- apply( estimates.is.one.founder.per.person, 1:2, function( .value ) {
          if( is.na( .value ) ) {
              0
          } else {
              .value
          }
      } );
       
        ## EVALUATION
        
        # sum.correct.among.one.founder.people <-
        #   apply( estimates.is.one.founder.per.person[ as.logical( gold.is.one.founder.per.person ),  ], 2, sum );
        # sum.incorrect.among.one.founder.people <-
        #   apply( estimates.is.one.founder.per.person[ as.logical( gold.is.one.founder.per.person ),  ], 2, function( .row ) { sum( 1-.row ) } );
        # sum.correct.among.multiple.founder.people <-
        #   apply( estimates.is.one.founder.per.person[ !as.logical( gold.is.one.founder.per.person ),  ], 2, function( .row ) { sum( 1-.row ) } );
        # sum.incorrect.among.multiple.founder.people <-
        #   apply( estimates.is.one.founder.per.person[ !as.logical( gold.is.one.founder.per.person ),  ], 2, sum );
        
        isMultiple.aucs <- 
            sapply( 1:ncol( estimates.is.one.founder.per.person ), function( .col.i ) {
                #print( .col.i );
                performance( prediction( as.numeric( estimates.is.one.founder.per.person[ , .col.i ] ), gold.is.one.founder.per.person ), measure = "auc" )@y.values[[ 1 ]];
            } );
        names( isMultiple.aucs ) <- colnames( estimates.is.one.founder.per.person );

        return( isMultiple.aucs );
    } # evaluate.results.per.ppt (..)
    
    evaluate.is.multiple.for.region.and.time <- function ( the.region, the.time ) {
        identify.founders.tab.file <- paste( RESULTS.DIR, results.dirname, the.region, the.time, "identify_founders.tab", sep = "/" );
        stopifnot( file.exists( identify.founders.tab.file ) );

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

        if( the.region == "v3" ) {
            gold.is.one.founder.per.person <-
                1-caprisa002.gold.is.multiple[ rownames( estimates.is.one.founder.per.person ) ];
        } else {
            stopifnot( ( the.region == "nflg" ) || ( length( grep( "rv217", the.region ) ) > 0 ) );
            gold.is.one.founder.per.person <-
                1-rv217.gold.is.multiple[ rownames( estimates.is.one.founder.per.person ) ];
        }
        
        if( is.na( partition.size ) ) {
            ## Now the issue is that there are multiple input files per ppt, eg for the NFLGs ther are often "LH" and "RH" files.  What to do?  The number of sequences varies.  Do a weighted average.
            ## Sometimes there are multiple entries for one ptid/sample, eg for the NFLGs there are often right half and left half (RH and LH) and sometimes additionally NFLG results.  If so, for something to be called single founder, all of the estimates (across regions, for a given test) must agree that it is one founder.
            estimates.is.one.founder.per.person <-
                compute.estimates.is.one.founder.per.person( estimates.is.one.founder );
            
            ## ERE I AM.  I'm going to add some estimates made by prediction using leave-one-out cross-validation.
            results.covars.per.person.with.extra.cols <-
                summarizeCovariatesOnePerParticipant( identify.founders.study );
            
          if( use.bounds ) {
              the.artificial.bounds <- getArtificialBounds( the.region, the.time, results.dirname );
              
              return( list( estimates.is.one.founder.per.person = estimates.is.one.founder.per.person, gold.is.one.founder.per.person = gold.is.one.founder.per.person, results.covars.per.person.with.extra.cols = results.covars.per.person.with.extra.cols, bounds = the.artificial.bounds, evaluated.results = evaluate.results.per.ppt( estimates.is.one.founder.per.person, gold.is.one.founder.per.person, results.covars.per.person.with.extra.cols, the.time, the.artificial.bounds ) ) );
          } else {
              return( list( estimates.is.one.founder.per.person = estimates.is.one.founder.per.person, gold.is.one.founder.per.person = gold.is.one.founder.per.person, results.covars.per.person.with.extra.cols = results.covars.per.person.with.extra.cols, evaluated.results = evaluate.results.per.ppt( estimates.is.one.founder.per.person, gold.is.one.founder.per.person, results.covars.per.person.with.extra.cols, the.time ) ) );
          }
        } else { # else !is.na( partition.size )
            ## ERE I AM.
            stop( "NOT IMPLEMENTED: EVALUATING PARTITIONS FOR ISMULTIPLE" );
        }
    } # evaluate.is.multiple.for.region.and.time ( the.region, the.time )

# Here "the.study" must be "rv217" or "caprisa002" or (special) "rv217_v3"
evaluateIsMultiple <- function ( the.study, output.dir = NULL, output.file = NULL, output.file.append = FALSE ) {


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
        return( evaluate.is.multiple.for.region.and.time( the.region, the.time ) );
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

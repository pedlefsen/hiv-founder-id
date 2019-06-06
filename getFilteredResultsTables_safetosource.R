#RESULTS.DIR <- "/fh/fast/edlefsen_p/bakeoff_analysis_results/";
#RESULTS.DIRNAME <- "raw_edited_20160216";

RESULTS.DIR <- "/fast/bakeoff_merged_analysis_sequences_results/results/";
#RESULTS.DIR <- "/fast/bakeoff_merged_analysis_sequences_results_2019/results/";
RESULTS.DIRNAME <- "raw_fixed";

THE.RESULTS.DIR <- RESULTS.DIR; # to avoid "promise already under evaluation" errors

repeatedRowsToColumns <- function ( the.matrix, pattern = "(?:glm|lasso|step).*.validation.results." ) {
    rownames.sans.patterns <- gsub( paste( "^(.*?)", pattern, "(.*)$", sep = "" ), "\\1\\2", rownames( the.matrix ) )
    pattern.matches.by.row <- gsub( paste( "^.*?(", pattern, ").*$", sep = "" ), "\\1", rownames( the.matrix ) );
    pattern.matches.by.row[ rownames.sans.patterns == pattern.matches.by.row ] <- "";

    result.submatrices <- lapply( unique( pattern.matches.by.row ), function ( .pattern ) {
        .rm <- the.matrix[ pattern.matches.by.row == .pattern, , drop = FALSE ];
        colnames( .rm ) <- paste( .pattern, colnames( the.matrix ), sep = "" );
        rownames( .rm ) <-
            gsub( paste( "^(.*?)", pattern, "(.*)$", sep = "" ), "\\1\\2", rownames( .rm ) );
        return( .rm );
    } );
    names( result.submatrices ) <- unique( pattern.matches.by.row );
    
    result.matrix <- matrix( NA, nrow = length( unique( rownames.sans.patterns ) ), ncol = length( result.submatrices ) * ncol( the.matrix ) );
    rownames( result.matrix ) <- unique( rownames.sans.patterns );
    colnames( result.matrix ) <- sapply( colnames( the.matrix ), function( .suffix ) { paste( names( result.submatrices ), .suffix, sep = "" ); } );
    .result.ignored <- lapply( result.submatrices, function( .result.submatrix ) {
        result.matrix[ rownames( .result.submatrix ), colnames( .result.submatrix ) ] <<-
            as.matrix( .result.submatrix );
        lapply( colnames( .result.submatrix ), function( .column ) {
          return( NULL );
        } )
        return( NULL );
    } );
    return( result.matrix );
} # repeatedRowsToColumns ( .. )

### Read results tables.
## setting the.time to "1m.6m" will return pooled results over those times
## rowname.pattern.map is a list of elements to be used in gsubs on the rownames of the results, iteratively in order.
getFilteredResultsTables <- function (
    out.tab.file.suffix, the.region, the.time, the.bounds.type = "unbounded", to.region = NULL, RESULTS.DIR = THE.RESULTS.DIR, results.dirname = RESULTS.DIRNAME, zeroNAs = TRUE, sort.column = "rmse", column.pattern = NA, rowname.pattern.map = list( "\\.(days|time)\\.est" = "", "\\.mut\\.rate\\.coef" = "", "multifounder\\." = "(w/in clusts) ", "Synonymous\\." = "(syn) ", "is\\.poisson" = "fits", "is\\.starlike" = "star-like", "is.one.founder" = "single-founder", "\\." = " " )
) {
    ## HACK: the isMultiple results don't have a zeroNAs option.
    if( length( grep( "sMultiple", out.tab.file.suffix ) ) > 0 ) {
        zeroNAs <- FALSE;
    }
    
    if( is.null( to.region ) || is.na( to.region ) ) {
        ## if to.region is not defined then it means we should use the single-region results.
        infile.results <- paste( RESULTS.DIR, results.dirname, "/", the.region, "/", the.time, "/", the.bounds.type, out.tab.file.suffix, sep = "" );
    } else { # is.null( to.region ) .. else ..
        ## if to.region is defined then it means we should use the pooled-over-multiple-regions results.
        from.region <- the.region;
        infile.results <- paste( RESULTS.DIR, results.dirname, "/", from.region, "_and_", to.region, "_", the.time, "_", the.bounds.type, out.tab.file.suffix, sep = "" )
    } # End if is.null( to.region ) .. else ..
    
    results.in <- read.table( infile.results, sep = "\t" );
    if( zeroNAs ) {
        results <- results.in[ , grep( "zeroNAs$", colnames( results.in ) ), drop = FALSE ];
        stopifnot( ncol( results ) > 1 );
        colnames( results ) <-
            gsub( "\\.zeroNAs$", "", colnames( results ) );
    } else {
        results <- results.in[ , grep( "zeroNAs$", colnames( results.in ), invert = TRUE ), drop = FALSE ];
    }

    # Filter out anything that's impossible -- so that's all the deterministic bounds as well as any non-matching time.  Also we use 30week snever 20weeks for the 6 month time point.
    if( the.time == "1m" ) {
        the.times.it.aint <- c( "20weeks", "30weeks", "1m5weeks_6m30weeks", "sixmonths", "1monemonth_6msixmonths" );
    } else if( the.time == "1m6m" ) {
        the.times.it.aint <- c( "20weeks", "1m5weeks_6m30weeks", "1monemonth_6msixmonths" );
    } else if( the.time == "6m" ) {
        the.times.it.aint <- c( "5weeks", "20weeks", "1m5weeks_6m30weeks", "onemonth", "1monemonth_6msixmonths" );
    } else if( the.time == "1m.6m" ) {
        the.times.it.aint <- c( "\\.5weeks", "\\.30weeks", "20weeks", "\\.onemonth", "\\.sixmonths" );
    }
    ## Also exclude deterministic bounds
    ## Also exclude "lower" and "upper" versions of results, which appear to be redundant.
    ## Also exclude "DS.Starphy" versions of results, which are redundant w PFitter (here, because it's just the est / mut rate coef)
    ## Also exclude all "Starphy" versions of results, which are close enough to redundant w PFitter that it's not worth cluttering the output.
    results.filtered <- results[ grep( paste( c( "(DS)?Star[Pp]hy(Test)?", "DS\\.Star[Pp]hy", "lower", "upper", "deterministic", the.times.it.aint ), collapse = "|" ), rownames( results ), invert = TRUE ), , drop = FALSE ];

    ## Maybe also exclude some columns.
    if( !is.null( column.pattern ) && !is.na( column.pattern ) && ( column.pattern != "" ) ) {
        results.filtered <- results.filtered[ , grep( column.pattern, colnames( results.filtered ) ), drop = FALSE ];
    }
    
    ## Maybe also rename some rows.
    if( !is.null( rowname.pattern.map ) && !is.na( rowname.pattern.map ) ) {
        for( .pattern in names( rowname.pattern.map ) ) {
            print( paste( "renaming rows according to pattern '", .pattern, "' => '", unlist( rowname.pattern.map[ .pattern ] ), "'", sep = "" ) );
            .rowname.matches <- grep( .pattern, rownames( results.filtered ), value = TRUE );
            if( length( .rowname.matches ) == 0 ) {
                next;
            }
            # Now fix 'em.
            rownames( results.filtered ) <- gsub( .pattern, unlist( rowname.pattern.map[ .pattern ] ), rownames( results.filtered ) );
        } # End foreach .pattern
    }
    
    if( !is.null( sort.column ) && ( length( sort.column ) > 0 ) && !is.na( sort.column ) ) {
        stopifnot( length( sort.column ) == 1 );
        stopifnot( sort.column %in% names( results.filtered ) );
        results.filtered.sorted <-
            results.filtered[ order( results.filtered[[ sort.column ]] ), , drop = FALSE ];
        return( repeatedRowsToColumns( results.filtered.sorted ) );
    }

    return( repeatedRowsToColumns( results.filtered ) );
} # getFilteredResultsTables (..)

### Get uses of parameters aggregated over lasso runs. see also 
## out.file.prefix should be "isMultiple" or "Timings"/
getFilteredLassoUsageTables <- function (
    out.file.prefix, the.region, the.time, the.bounds.type = "unbounded", to.region = NULL, RESULTS.DIR = THE.RESULTS.DIR, results.dirname = RESULTS.DIRNAME, column.pattern = NA, rowname.pattern.map = list( "\\.(days|time)\\.est" = "", "\\.mut\\.rate\\.coef" = "", "multifounder\\." = "(w/in clusts) ", "Synonymous\\." = "(syn) ", "is\\.poisson" = "fits", "is\\.starlike" = "star-like", "is.one.founder" = "single-founder", "\\." = " " ), colname.pattern.map = list( "inf\\.sites" = "InSites", "multifounder\\." = "(w/in clusts) ", "Synonymous\\." = "(syn) ", "is\\.poisson" = "fits", "is\\.starlike" = "star-like", "is.one.founder" = "single-founder", "\\.hd" = " HD", "\\." = " " )
) {
    results.by.region.and.time.Rda.filename <-
        paste( RESULTS.DIR, results.dirname, "/", out.file.prefix, ".results.by.region.and.time.Rda", sep = "" );
    load( file = results.by.region.and.time.Rda.filename );
    
    if( is.null( to.region ) || is.na( to.region ) ) {
        ## if to.region is not defined then it means we should use the single-region results.
        .lasso.coefs <- results.by.region.and.time[[ the.region ]][[ the.time ]][[ "evaluated.results" ]][[ the.bounds.type ]][[ "lasso.coefs" ]];
    } else { # is.null( to.region ) .. else ..
        ## if to.region is defined then it means we should use the pooled-over-multiple-regions results.
        from.region <- the.region;
        .lasso.coefs <- results.by.region.and.time[[ "results.across.regions.by.time" ]][[ the.region ]][[ to.region ]][[ the.time ]][[ "evaluated.results" ]][[ the.bounds.type ]][[ "lasso.coefs" ]];
    } # End if is.null( to.region ) .. else ..

    ## Special case if .lasso.coefs has options for eg "withbounds" then use that.
    if( "lasso.withbounds" %in% names( .lasso.coefs ) ) {
        .lasso.coefs <- .lasso.coefs[[ "lasso.withbounds" ]];
    }

    
## first index is ppt, second index is special estimator.
.mat.per.ppt <- 
    lapply( .lasso.coefs, function ( .lasso.coefs.for.ppt ) {
        #print( "hi" );
  .hi <- lapply( .lasso.coefs.for.ppt, function( .sparse.matrix ) {
    if( length( dim( .sparse.matrix ) ) == 0 ) {
      return( NA );
    }
        #print( "there" );
      .foo <- as.numeric( .sparse.matrix );
      names( .foo ) <- rownames( .sparse.matrix );
      return( .foo );
  } );
  all.evaluators <- names( .hi ) <- names( .lasso.coefs.for.ppt );
  ## This is all but the last covar listed, which is in the max position if it's there.
  .special.covar.position <- max( sapply( .hi, length ) );
  all.covars <- unique( unlist( lapply( .hi, function( .lst ) { if( length( .lst ) < .special.covar.position ) { names( .lst ) } else { names( .lst[ -.special.covar.position ] ) } } ) ) );
  .mat <- matrix( NA, nrow = length( all.covars ), ncol = length( all.evaluators ) );
  rownames( .mat ) <- all.covars;
  colnames( .mat ) <- all.evaluators;
  for( .evaluator in all.evaluators ) {
      #print( .evaluator );
      .n <- names( .hi[[ .evaluator ]] );
      #print( .n );
      #print( .hi[[ .evaluator ]] );
      ..hi <- unlist( .hi[[ .evaluator ]] );
      if( length( ..hi ) >= .special.covar.position ) {
          ..hi <- ..hi[ -.special.covar.position ];
      }
      #print(names( ..hi ) );
      .mat[ names( ..hi ), .evaluator ] <- ..hi;
  }
  .mat[ is.na( .mat ) ] <- 0;
  return( .mat );
} );
# This returns the average, over all models evaluated when a particular ppt is excluded, of the uses of each covariate by the lasso-selected model.
get.uses.by.ppt.for.evaluator <- function ( the.evaluator ) {
    .logical.mat <- do.call( rbind, lapply( .mat.per.ppt, function ( .mat.for.ppt ) { .rv <- as.logical( .mat.for.ppt[ , the.evaluator ] ); .rv[ is.na( .rv ) ] <- 0; return( .rv ); } ) );
    colnames( .logical.mat ) <- rownames( .mat.per.ppt[[1]] );
    return( .logical.mat );
} # get.uses.by.ppt.for.evaluator (..)

all.evaluators <- colnames( .mat.per.ppt[[1]] );
uses.by.evaluator <- sapply( all.evaluators, function ( the.evaluator ) {
     .uses.for.evaluator <- get.uses.by.ppt.for.evaluator( the.evaluator );
     .avg.uses.for.evaluator <- apply( .uses.for.evaluator, 2, mean, na.rm = TRUE );
     return( .avg.uses.for.evaluator );
 } );

    results <- t( uses.by.evaluator );
    
    # Filter out anything that's impossible -- so that's all the deterministic bounds as well as any non-matching time.  Also we use 30weeks never 20weeks for the 6 month time point.
    if( the.time == "1m" ) {
        the.times.it.aint <- c( "20weeks", "30weeks", "1m5weeks_6m30weeks", "sixmonths", "1monemonth_6msixmonths" );
    } else if( the.time == "1m6m" ) {
        the.times.it.aint <- c( "20weeks", "1m5weeks_6m30weeks", "1monemonth_6msixmonths" );
    } else if( the.time == "6m" ) {
        the.times.it.aint <- c( "5weeks", "20weeks", "1m5weeks_6m30weeks", "onemonth", "1monemonth_6msixmonths" );
    } else if( the.time == "1m.6m" ) {
        the.times.it.aint <- c( "\\.5weeks", "\\.30weeks", "20weeks", "\\.onemonth", "\\.sixmonths" );
    }
    ## Also exclude deterministic bounds
    ## Also exclude "lower" and "upper" versions of results, which appear to be redundant.
    ## Also exclude "DS.Starphy" versions of results, which are redundant w PFitter (here, because it's just the est / mut rate coef)
    ## Also exclude all "Starphy" versions of results, which are close enough to redundant w PFitter that it's not worth cluttering the output.
    results.filtered <- results[ grep( paste( c( "(DS)?Star[Pp]hy(Test)?", "DS\\.Star[Pp]hy", "lower", "upper", "deterministic", the.times.it.aint ), collapse = "|" ), rownames( results ), invert = TRUE ), , drop = FALSE ];

    ## Maybe also exclude some columns.
    if( !is.null( column.pattern ) && !is.na( column.pattern ) && ( column.pattern != "" ) ) {
        results.filtered <- results.filtered[ , grep( column.pattern, colnames( results.filtered ) ), drop = FALSE ];
    }

    ## Maybe also rename some cols.
    if( !is.null( colname.pattern.map ) && !is.na( colname.pattern.map ) ) {
        for( .pattern in names( colname.pattern.map ) ) {
            print( paste( "renaming cols according to pattern '", .pattern, "' => '", unlist( colname.pattern.map[ .pattern ] ), "'", sep = "" ) );
            .colname.matches <- grep( .pattern, colnames( results.filtered ), value = TRUE );
            if( length( .colname.matches ) == 0 ) {
                next;
            }
            # Now fix 'em.
            colnames( results.filtered ) <- gsub( .pattern, unlist( colname.pattern.map[ .pattern ] ), colnames( results.filtered ) );
        } # End foreach .pattern
    }
    
    ## Maybe also rename some rows.
    if( !is.null( rowname.pattern.map ) && !is.na( rowname.pattern.map ) ) {
        for( .pattern in names( rowname.pattern.map ) ) {
            print( paste( "renaming rows according to pattern '", .pattern, "' => '", unlist( rowname.pattern.map[ .pattern ] ), "'", sep = "" ) );
            .rowname.matches <- grep( .pattern, rownames( results.filtered ), value = TRUE );
            if( length( .rowname.matches ) == 0 ) {
                next;
            }
            # Now fix 'em.
            rownames( results.filtered ) <- gsub( .pattern, unlist( rowname.pattern.map[ .pattern ] ), rownames( results.filtered ) );
        } # End foreach .pattern
    }
    
    return( repeatedRowsToColumns( results.filtered ) );
} # getFilteredLassoUsageTables (..)

        get.formulas <- function ( results.by.region.and.time, .varname = "none", model.type = "glm", withbounds = TRUE, regions = c( "nflg", "v3" ), times = c( "1m", "6m" ), the.bounds.type = "unbounded" ) {
            if( withbounds ) {
                .withbounds.string <- paste( model.type, "withbounds", sep = "." );
            } else {
                .withbounds.string <- model.type;
            }
            if( length( regions ) == 2 ) {
                .results.for.region <- results.by.region.and.time[[3]][[1]][[1]];
            } else {
                .results.for.region <- results.by.region.and.time[[ regions ]];
            }
            if( length( times ) == 2 ) {
                the.time <- "1m.6m";
            } else {
                the.time <- times;
            }
            .results.by.removed.ptid <-
                .results.for.region[[ the.time ]][[ "evaluated.results" ]][[the.bounds.type]][[ paste( model.type, "formulas", sep = "." ) ]][[ .withbounds.string ]];
            if( !( .varname %in% colnames( .results.by.removed.ptid ) ) ) {
                ## To avoid a major crash, just revert to using another one as a template.
                .formulas <- .results.by.removed.ptid[ , 1, drop = FALSE ];
                .rv <- table( .formulas );
                warning( "MISSING VAR WHEN RETRIEVING FORMULAS! USED ARBITRARY OTHER VAR AS TEMPLATE, WITH 0 COUNTS!" );
                .rv[ 1:length( .rv ) ] <- 0;
                return( .rv );
            }
            
            .formulas <- .results.by.removed.ptid[ , .varname, drop = FALSE ];
            table( .formulas );
        } # get.formulas (..)

        get.uses <- function ( results.by.region.and.time, .varname = "none", withbounds = TRUE, regions = c( "nflg", "v3" ), times = c( "1m", "6m" ) ) {
            if( withbounds ) {
                .withbounds.string <- "lasso.withbounds";
            } else {
                .withbounds.string <- "lasso";
            }
            if( length( regions ) == 2 ) {
                .results.for.region <- results.by.region.and.time[[3]][[1]][[1]];
            } else {
                .results.for.region <- results.by.region.and.time[[ regions ]];
            }
            if( length( times ) == 2 ) {
                the.time <- "1m.6m";
            } else {
                the.time <- times;
            }
            .results.by.removed.ptid <-
                .results.for.region[[ the.time ]][[ "evaluated.results" ]][["unbounded"]][[ "lasso.coefs" ]][[ .withbounds.string ]];
            .uses <- lapply( 1:length( .results.by.removed.ptid ), function( .i ) { .dgCMatrix <- .results.by.removed.ptid[[.i]][[.varname]]; .rv <- as.logical( .dgCMatrix ); names( .rv ) <- rownames( .dgCMatrix ); return( .rv ); } );
            table( names( which( unlist( .uses ) ) ) );
        } # get.uses (..)

        get.rmses <- function ( results.by.region.and.time, evaluate.regions = train.regions, evaluate.times = train.times, train.regions = c( "nflg", "v3" ), train.times = c( "1m", "6m" ), the.bound = "sampledwitdth_uniform", zeroNAs = TRUE ) {
            if( zeroNAs ) {
                rmse.stat <- "rmse.zeroNAs";
            } else {
                rmse.stat <- "rmse";
            }
            if( is.null( the.bound ) || is.na( the.bound ) ) {
                the.bound <- "unbounded";
            }
            if( length( train.times ) == 2 ) {
                train.time <- "1m.6m";
                if( the.bound != "unbounded" ) {
                    the.bound <- "sampledwidth_uniform_1mmtn003_6mhvtn502";
                }
            } else {
                train.time <- train.times;
                if( train.time == "6m" ) {
                    if( the.bound != "unbounded" ) {
                        the.bound <- "sampledwidth_uniform_hvtn502";
                    }
                    stopifnot( length( evaluate.times ) == 1 );
                    stopifnot( evaluate.times == "6m" );
                } else {
                    if( the.bound != "unbounded" ) {
                        the.bound <- "sampledwidth_uniform_mtn003";
                    }
                    stopifnot( length( evaluate.times ) == 1 );
                    stopifnot( evaluate.times == "1m" );
                }
            }
            if( length( evaluate.times ) == 2 ) {
                # Leave the.bound alone.
            } else if( length( train.times ) == 2 ) {
                the.bound <- paste( the.bound, evaluate.times, sep = "." );
            }
            if( length( train.regions ) == 2 ) {
                .results.for.region <- results.by.region.and.time[[3]][[1]][[1]];
              if( length( evaluate.regions ) == 2 ) {
                  # Leave the.bound alone, then.
              } else {
                the.bound <- paste( the.bound, evaluate.regions, sep = "." );
              }
            } else {
              if( length( evaluate.regions ) == 2 ) {
                  stop( "can't evaluate more regions than trained" );
              } else if( evaluate.regions != train.regions ) {
                  stop( "can't evaluate a different region than trained" );
              }
              .results.for.region <- results.by.region.and.time[[ train.regions ]];
            }
            .lst <-
                sort( unlist( .results.for.region[[ train.time ]][[ "evaluated.results" ]][[ the.bound ]][[ rmse.stat ]] ), decreasing = T );
            ## TODO: REMOVE. Temporary.
            .lst <- .lst[ grep( "(one|six)month", names( .lst ), value = TRUE, invert = TRUE ) ];
            return( .lst );
        } # get.rmses (..)

        # This returns biases in order of rmses.
        get.biases <- function ( results.by.region.and.time, evaluate.regions = train.regions, evaluate.times = train.times, train.regions = c( "nflg", "v3" ), train.times = c( "1m", "6m" ), the.bound = "sampledwitdth_uniform", zeroNAs = TRUE ) {
            if( zeroNAs ) {
                bias.stat <- "bias.zeroNAs";
                rmse.stat <- "rmse.zeroNAs";
            } else {
                bias.stat <- "bias";
                rmse.stat <- "rmse";
            }
            if( is.null( the.bound ) || is.na( the.bound ) ) {
                the.bound <- "unbounded";
            }
            if( length( train.times ) == 2 ) {
                train.time <- "1m.6m";
                if( the.bound != "unbounded" ) {
                    the.bound <- "sampledwidth_uniform_1mmtn003_6mhvtn502";
                }
            } else {
                train.time <- train.times;
                if( train.time == "6m" ) {
                    if( the.bound != "unbounded" ) {
                        the.bound <- "sampledwidth_uniform_hvtn502";
                    }
                    stopifnot( length( evaluate.times ) == 1 );
                    stopifnot( evaluate.times == "6m" );
                } else {
                    if( the.bound != "unbounded" ) {
                        the.bound <- "sampledwidth_uniform_mtn003";
                    }
                    stopifnot( length( evaluate.times ) == 1 );
                    stopifnot( evaluate.times == "1m" );
                }
            }
            if( length( evaluate.times ) == 2 ) {
                # Leave the.bound alone.
            } else if( length( train.times ) == 2 ) {
                the.bound <- paste( the.bound, evaluate.times, sep = "." );
            }
            if( length( train.regions ) == 2 ) {
                .results.for.region <- results.by.region.and.time[[3]][[1]][[1]];
              if( length( evaluate.regions ) == 2 ) {
                  # Leave the.bound alone, then.
              } else {
                the.bound <- paste( the.bound, evaluate.regions, sep = "." );
              }
            } else {
              if( length( evaluate.regions ) == 2 ) {
                  stop( "can't evaluate more regions than trained" );
              } else if( evaluate.regions != train.regions ) {
                  stop( "can't evaluate a different region than trained" );
              }
              .results.for.region <- results.by.region.and.time[[ train.regions ]];
            }
            .lst <- unlist( .results.for.region[[ train.time ]][[ "evaluated.results" ]][[ the.bound ]][[ "bias.stat" ]] )[ order( unlist( .results.for.region[[ train.time ]][[ "evaluated.results" ]][[ the.bound ]][[ rmse.stat ]] ), decreasing = T ) ];
            ## TODO: REMOVE. Temporary.
            .lst <- .lst[ grep( "(one|six)month", names( .lst ), value = TRUE, invert = TRUE ) ];
            return( .lst );
        } # get.biases (..)

        get.bias.and.rmse <- function ( ... ) { cbind( bias = get.biases( ... ), rmse = get.rmses( ... ) ) }
        
        evaluate.specific.timings.model <-
          function ( results.by.region.and.time, model.vars, .include.intercept = FALSE, step = FALSE, train.regions = c( "nflg", "v3" ), train.times = c( "1m", "6m" ) ) {

            if( include.intercept ) {
                .formula <- as.formula( paste( "days.since.infection ~ ", paste( model.vars, collapse = "+" ) ) );
            } else {
                .formula <- as.formula( paste( "days.since.infection ~ 0 + ", paste( model.vars, collapse = "+" ) ) );
            }
            evaluate.specific.timings.model.formula( results.by.region.and.time, .formula, step = step, train.regions = train.regions, train.times = train.times );
        } # evaluate.specific.timings.model (..)
        evaluate.specific.timings.model.formula <-
          function ( results.by.region.and.time, .formula, step = FALSE, regions = c( "nflg", "v3" ), times = c( "1m", "6m" ) ) {
            if( length( times ) == 2 ) {
                the.time <- "1m.6m";
            } else {
                the.time <- times;
            }
            if( length( regions ) == 2 ) {
                .results.for.region <- results.by.region.and.time[[3]][[1]][[1]];
            } else {
                .results.for.region <- results.by.region.and.time[[ regions ]];
            }
            results.covars.per.person.df <-
                data.frame( .results.for.region[[ the.time ]][[ "results.covars.per.person.with.extra.cols" ]] );
            ## DO NOT Undo conversion of the colnames (X is added before "6m.not.1m").  We want it to be called "X6m.not.1m" so it can work in the regression formulas.
            #colnames( results.covars.per.person.df ) <- colnames( results.covars.per.person.with.extra.cols );

            regression.df <- cbind( data.frame( days.since.infection = .results.for.region[[ the.time ]][["days.since.infection" ]][ rownames( results.covars.per.person.df ) ] ), .results.for.region[[ the.time ]][["results.per.person"]][ rownames( results.covars.per.person.df ), , drop = FALSE ], results.covars.per.person.df, .results.for.region[[ the.time ]][[ "bounds" ]] );

            .lm <-
                suppressWarnings( lm( .formula, data = regression.df ) );
            if( step ) {
                .step.rv <- step( .lm ); # Stepwise regression, both forward and backward.
                return( .step.rv );            
            }
            return( .lm );
        } # evaluate.specific.timings.model.formula (..)
        
    evaluate.specific.isMultiple.model <-
        function ( results.by.region.and.time, model.vars, include.intercept = TRUE, step = FALSE, regions = c( "nflg", "v3" ), times = c( "1m", "6m" ) ) {
        if( include.intercept ) {
            .formula <- as.formula( paste( "is.one.founder ~ ", paste( model.vars, collapse = "+" ) ) );
        } else {
            .formula <- as.formula( paste( "is.one.founder ~ 0 + ", paste( model.vars, collapse = "+" ) ) );
        }
        return( evaluate.specific.isMultiple.model.formula( results.by.region.and.time, .formula, step = step, regions = regions, times = times ) );
        } # evaluate.specific.isMultiple.model (..)

        evaluate.specific.isMultiple.model.formula <-
          function ( results.by.region.and.time, .formula, step = FALSE, regions = c( "nflg", "v3" ), times = c( "1m", "6m" ) ) {
            if( length( times ) == 2 ) {
                the.time <- "1m.6m";
            } else {
                the.time <- times;
            }
            if( length( regions ) == 2 ) {
                .results.for.region <- results.by.region.and.time[[3]][[1]][[1]];
            } else {
                .results.for.region <- results.by.region.and.time[[ regions ]];
            }
            results.covars.per.person.df <-
                data.frame( .results.for.region[[ the.time ]][[ "results.covars.per.person.with.extra.cols" ]] );
            ## DO NOT Undo conversion of the colnames (X is added before "6m.not.1m").  We want it to be called "X6m.not.1m" so it can work in the regression formulas.
            #colnames( results.covars.per.person.df ) <- colnames( results.covars.per.person.with.extra.cols );

        regression.df <- cbind( data.frame( is.one.founder = results.by.region.and.time[[3]][[1]][[1]][[1]][["gold.is.one.founder.per.person" ]][ rownames( results.covars.per.person.df ) ] ), results.covars.per.person.df, results.by.region.and.time[[3]][[1]][[1]][[1]][["bounds" ]] );
        .lm <-
            suppressWarnings( glm( .formula, family = "binomial", data = regression.df ) );
        if( step ) {
            .step.rv <- step( .lm ); # Stepwise regression, both forward and backward.
            return( .step.rv );            
        }
       return( .lm );
        } # evaluate.specific.isMultiple.model.formula (..)

        compute.pearson.R.of.specific.isMultiple.predictor.with.gold.standard <-
            function ( .lm, use.residuals = FALSE ) {
                if( use.residuals ) {
                    cor( .lm$model$is.one.founder, residuals( .lm ) )
                } else {
                    .var <- .lm$model[ , setdiff( names( .lm$model ), "is.one.founder" ) ];
                    if( ( length( .var ) == 1 ) && is.na( .var  ) ) {
                        return( NA );
                    } else {
                        return( cor( .lm$model$is.one.founder, .var ) );
                    }
                }
        } # compute.pearson.R.of.specific.isMultiple.predictor.with.gold.standard (..)

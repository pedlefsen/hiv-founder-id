repeatedRowsToColumns <- function ( the.matrix, pattern = "(?:glm|lasso).*.validation.results." ) {
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
getFilteredResultsTables <- function (
    out.tab.file.suffix, the.region, the.time, the.bounds.type = "unbounded", to.region = NULL, results.dirname = "raw_edited_20160216", zeroNAs = TRUE, sort.column = "rmse", column.pattern = NA, rowname.pattern.map = list( "\\.(days|time)\\.est" = "", "\\.mut\\.rate\\.coef" = "", "multifounder\\." = "(w/in clusts) ", "Synonymous\\." = "(syn) ", "is\\.poisson" = "fits", "is\\.starlike" = "star-like", "is.one.founder" = "single-founder", "\\." = " " )
) {
    ## HACK: the isMultiple results don't have a zeroNAs option.
    if( length( grep( "sMultiple", out.tab.file.suffix ) ) > 0 ) {
        zeroNAs <- FALSE;
    }
    
    if( is.null( to.region ) || is.na( to.region ) ) {
        ## if to.region is not defined then it means we should use the single-region results.
        infile.results <- paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/", the.bounds.type, out.tab.file.suffix, sep = "" );
    } else { # is.null( to.region ) .. else ..
        ## if to.region is defined then it means we should use the pooled-over-multiple-regions results.
        from.region <- the.region;
        infile.results <- paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", from.region, "_and_", to.region, "_", the.time, "_", the.bounds.type, out.tab.file.suffix, sep = "" )
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
    results.filtered <- results[ grep( paste( c( "Star[Pp]hy", "DS\\.Star[Pp]hy", "lower", "upper", "deterministic", the.times.it.aint ), collapse = "|" ), rownames( results ), invert = TRUE ), , drop = FALSE ];

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

### Get uses of parameters over lasso runs.
## out.file.prefix should be "isMultiple" or "Timings"/
getFilteredLassoUsageTables <- function (
    out.file.prefix, the.region, the.time, the.bounds.type = "unbounded", to.region = NULL, results.dirname = "raw_edited_20160216", column.pattern = NA, rowname.pattern.map = list( "\\.(days|time)\\.est" = "", "\\.mut\\.rate\\.coef" = "", "multifounder\\." = "(w/in clusts) ", "Synonymous\\." = "(syn) ", "is\\.poisson" = "fits", "is\\.starlike" = "star-like", "is.one.founder" = "single-founder", "\\." = " " ), colname.pattern.map = list( "inf\\.sites" = "InSites", "multifounder\\." = "(w/in clusts) ", "Synonymous\\." = "(syn) ", "is\\.poisson" = "fits", "is\\.starlike" = "star-like", "is.one.founder" = "single-founder", "\\.hd" = " HD", "\\." = " " )
) {
    RESULTS.DIR <- "/fh/fast/edlefsen_p/bakeoff_analysis_results/";
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
    results.filtered <- results[ grep( paste( c( "Star[Pp]hy", "DS\\.Star[Pp]hy", "lower", "upper", "deterministic", the.times.it.aint ), collapse = "|" ), rownames( results ), invert = TRUE ), , drop = FALSE ];

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

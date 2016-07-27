writeResultsTables <- function ( results.by.region.and.time, out.tab.file.suffix, regions, results.are.bounded = TRUE ) {
    ## Note that there are now special entries in results.by.region.and.time that are not regions (under "results.across.regions.by.time") -- these are comparisons across the two (main) regions.

  if( results.are.bounded ) {
    getEvaluatedResults <- function ( results.by.region.and.time, the.region, the.time, bounds.names ) {
      results.by.region.and.time[[ the.region ]][[ the.time ]][[ "evaluated.results" ]][ bounds.names ]
    }
    getEvaluatedResultsAcrossRegions <- function ( results.by.region.and.time, from.region, to.region, the.time, bounds.names ) {
      results.by.region.and.time[[ "results.across.regions.by.time" ]][[ from.region ]][[ to.region ]][[ the.time ]][[ "evaluated.results" ]][ bounds.names ]
    }
  } else {
    getEvaluatedResults <- function ( results.by.region.and.time, the.region, the.time ) {
      results.by.region.and.time[[ the.region ]][[ the.time ]][[ "evaluated.results" ]]
    }
    getEvaluatedResultsAcrossRegions <- function ( results.by.region.and.time, from.region, to.region, the.time ) {
      results.by.region.and.time[[ "results.across.regions.by.time" ]][[ from.region ]][[ to.region ]][[ the.time ]][[ "evaluated.results" ]]
    }
  }
  
    # Make a table out of it. (one per study).  See below for making tables from results.across.regions.by.time.
    results.table.by.region.and.time.and.bounds.type <-
        lapply( regions, function( the.region ) {
            .rv <- 
                lapply( names( results.by.region.and.time[[ the.region ]] ), function( the.time ) {
                  # Results are bounded for the timings but not the is.single results.
                  if( results.are.bounded ) {
                    # the "bounds" here actually also may include "glm.fit.statistics", which is not really a bound.
                    ..relevant.bounds <- setdiff( names( results.by.region.and.time[[ the.region ]][[ the.time ]][[ "evaluated.results" ]] ), "glm.fit.statistics" );
                    ..rv <- 
                      lapply( getEvaluatedResults( results.by.region.and.time, the.region, the.time, ..relevant.bounds ), function( results.for.bounds.type ) {
                        ## TODO: REMOVE
                        ## Special workaround for bug fixed 19 March 2016
                                                if( class( results.for.bounds.type ) == "matrix" ) {
                                                  stopifnot( ncol( results.for.bounds.type ) == 1 );
                                                  results.for.bounds.type <- list( AUC = results.for.bounds.type[ , 1 ] );
                                                }
                                                if( "lasso.coefs" %in% names( results.for.bounds.type ) ) {
                                                  results.for.bounds.type <- results.for.bounds.type[ names( results.for.bounds.type ) != "lasso.coefs" ];
                                                }
                                                if( "glm.formulas" %in% names( results.for.bounds.type ) ) {
                                                  results.for.bounds.type <- results.for.bounds.type[ names( results.for.bounds.type ) != "glm.formulas" ];
                                                }
                                                if( "lasso.formulas" %in% names( results.for.bounds.type ) ) {
                                                  results.for.bounds.type <- results.for.bounds.type[ names( results.for.bounds.type ) != "lasso.formulas" ];
                                                }
                                                sapply( results.for.bounds.type, function( results.list ) { results.list } );
                      } );
                    names( ..rv ) <- ..relevant.bounds;
                    return( ..rv );
                  } else { # if results.are.bounded .. else ..
                    .evaluated.results <-
                      getEvaluatedResults( results.by.region.and.time, the.region, the.time );
                    ..rv <-
                      sapply( .evaluated.results, function( results.list ) { results.list } );
                    return( list( unbounded = as.matrix( ..rv ) ) );
                  } # End if results.are.bounded .. else ..
                } );
            names( .rv ) <- names( results.by.region.and.time[[ the.region ]] );
            return( .rv );
        } );
    names( results.table.by.region.and.time.and.bounds.type ) <- regions;
    
    ## Write these out.
    .result.ignored <- sapply( regions, function ( the.region ) {
        ..result.ignored <- 
        sapply( names( results.by.region.and.time[[ the.region ]] ), function ( the.time ) {
          .bounds.types <- names( results.table.by.region.and.time.and.bounds.type[[ the.region ]][[ the.time ]] );
          ...result.ignored <- 
              sapply( .bounds.types, function ( the.bounds.type ) {
                  out.dir <- paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/", sep = "" );
                  dir.create( out.dir, recursive = TRUE, showWarnings = FALSE );
              out.file <- paste( out.dir, the.bounds.type, out.tab.file.suffix, sep = "" );
              ## TODO: REMOVE
              print( paste( the.bounds.type, out.file ) );
              ..tbl <- results.table.by.region.and.time.and.bounds.type[[ the.region ]][[ the.time ]][[ the.bounds.type ]];
              #print( ..tbl );
              .tbl <-
                apply( ..tbl, 1:2, function( .x ) { sprintf( "%0.2f", .x ) } );
              #print( .tbl );
              write.table( .tbl, quote = FALSE, file = out.file, sep = "\t" );
              return( NULL );
            } );
          return( NULL );
        } );
        return( NULL );
    } );

    ## Next: write out the regions tables.
    # Make a table out of each result in results.across.regions.by.time, too.
    results.table.across.regions.by.time.and.bounds.type <- 
        lapply( regions[ -length( regions ) ], function( from.region ) {
            .rv <- 
                lapply( names( results.by.region.and.time[[ "results.across.regions.by.time" ]][[ from.region ]] ), function( to.region ) {
                    ..rv <- 
                        lapply( names( results.by.region.and.time[[ "results.across.regions.by.time" ]][[ from.region ]][[ to.region ]] ), function( the.time ) {
                          # Results are bounded for the timings but not the is.single results.
                          if( results.are.bounded ) {
                            # the "bounds" here actually also may include "glm.fit.statistics", which is not really a bound.
                            ..relevant.bounds <- setdiff( names( results.by.region.and.time[[ "results.across.regions.by.time" ]][[ from.region ]][[ to.region ]][[ the.time ]][[ "evaluated.results" ]] ), "glm.fit.statistics" );
                  
                            ..rv <- 
                              lapply( getEvaluatedResultsAcrossRegions( results.by.region.and.time, from.region, to.region, the.time, ..relevant.bounds ), function( results.for.bounds.type ) {
                        ## TODO: REMOVE
                        ## Special workaround for bug fixed 19 March 2016
                                                if( class( results.for.bounds.type ) == "matrix" ) {
                                                  stopifnot( ncol( results.for.bounds.type ) == 1 );
                                                  results.for.bounds.type <- list( AUC = results.for.bounds.type[ , 1 ] );
                                                }
                                                if( "lasso.coefs" %in% names( results.for.bounds.type ) ) {
                                                  results.for.bounds.type <- results.for.bounds.type[ names( results.for.bounds.type ) != "lasso.coefs" ];
                                                }
                                                if( "glm.formulas" %in% names( results.for.bounds.type ) ) {
                                                  results.for.bounds.type <- results.for.bounds.type[ names( results.for.bounds.type ) != "glm.formulas" ];
                                                }
                                                if( "lasso.formulas" %in% names( results.for.bounds.type ) ) {
                                                  results.for.bounds.type <- results.for.bounds.type[ names( results.for.bounds.type ) != "lasso.formulas" ];
                                                }
                                sapply( results.for.bounds.type, function( results.list ) { results.list } );
                              } );
                            names( ..rv ) <- ..relevant.bounds;
                            return( ..rv );
                          } else { # if results.are.bounded .. else ..
                            .evaluated.results <-
                              getEvaluatedResultsAcrossRegions( results.by.region.and.time, from.region, to.region, the.time );
                            ..rv <-
                              sapply( .evaluated.results, function( results.list ) { results.list } );
                            return( list( unbounded = as.matrix( ..rv ) ) );
                          } # End if results.are.bounded .. else ..

                        } );
                  names( ..rv ) <- names( results.by.region.and.time[[ "results.across.regions.by.time" ]][[ from.region ]][[ to.region ]] )
                  return( ..rv );
                } );
            names( .rv ) <- names( results.by.region.and.time[[ "results.across.regions.by.time" ]][[ from.region ]] );
            return( .rv );
        } );
    names( results.table.across.regions.by.time.and.bounds.type ) <- regions[ -length( regions ) ];
    
    ## Write these pooled-over-regions results out, too.
    .result.ignored <- sapply( regions[ -length( regions ) ], function ( from.region ) {
        ..result.ignored <- sapply( names( results.table.across.regions.by.time.and.bounds.type[[ from.region ]] ), function ( to.region ) {
            ...result.ignored <- 
        sapply( names( results.table.across.regions.by.time.and.bounds.type[[ from.region ]][[ to.region ]] ), function ( the.time ) {
          .bounds.types <- names( results.table.across.regions.by.time.and.bounds.type[[ from.region ]][[ to.region ]][[ the.time ]] );
          ....result.ignored <- 
            sapply( .bounds.types, function ( the.bounds.type ) {
              out.file <- paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", from.region, "_and_", to.region, "_", the.time, "_", the.bounds.type, out.tab.file.suffix, sep = "" );
              ## TODO: REMOVE
              print( paste( the.bounds.type, out.file ) );
              .tbl <-
                apply( results.table.across.regions.by.time.and.bounds.type[[ from.region ]][[ to.region ]][[ the.time ]][[ the.bounds.type ]], 1:2, function( .x ) { sprintf( "%0.2f", .x ) } );
              #print( .tbl );
              write.table( .tbl, quote = FALSE, file = out.file, sep = "\t" );
              return( NULL );
            } );
              return( NULL );
            } );
          return( NULL );
        } );
        return( NULL );
    } );
    
    return( NULL );
} # writeResultsTables (..)

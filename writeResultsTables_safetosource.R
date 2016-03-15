writeResultsTables <- function ( results.by.region.and.time, out.tab.file.suffix, regions, times ) {
    ## Note that there are now special entries in results.by.region.and.time that are not regions (under "results.across.regions.by.time") -- these are comparisons across the two (main) regions.
    
    # Make a table out of it. (one per study).  See below for making tables from results.across.regions.by.time.
    results.table.by.region.and.time.and.bounds.type <-
        lapply( regions, function( the.region ) {
            .rv <- 
                lapply( names( results.by.region.and.time[[ the.region ]] ), function( the.time ) {
                  ..rv <- 
                    lapply( results.by.region.and.time[[ the.region ]][[ the.time ]][[ "evaluated.results" ]], function( results.by.bounds.type ) {
                      sapply( results.by.bounds.type, function( results.list ) { results.list } );
                    } );
                  names( ..rv ) <- names( results.by.region.and.time[[ the.region ]][[ the.time ]][[ "evaluated.results" ]] )
                  return( ..rv );
                } );
            names( .rv ) <- names( results.by.region.and.time[[ the.region ]] );
            return( .rv );
        } );
    names( results.table.by.region.and.time.and.bounds.type ) <- regions;
    
    ## Write these out.
    .result.ignored <- sapply( regions, function ( the.region ) {
        ..result.ignored <- 
        sapply( times, function ( the.time ) {
          .bounds.types <- names( results.table.by.region.and.time.and.bounds.type[[ the.region ]][[ the.time ]] );
          ...result.ignored <- 
            sapply( .bounds.types, function ( the.bounds.type ) {
              out.file <- paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/", the.bounds.type, "_evaluateIsMultiple.tab", sep = "" );
              ## TODO: REMOVE
              print( the.bounds.type );
              .tbl <-
                apply( results.table.by.region.and.time.and.bounds.type[[ the.region ]][[ the.time ]][[ the.bounds.type ]], 1:2, function( .x ) { sprintf( "%0.2f", .x ) } );
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
                  ...rv <- 
                    lapply( results.by.region.and.time[[ "results.across.regions.by.time" ]][[ from.region ]][[ to.region ]][[ the.time ]][[ "evaluated.results" ]], function( results.by.bounds.type ) {
                      sapply( results.by.bounds.type, function( results.list ) { results.list } );
                    } );
                  names( ...rv ) <- names( results.by.region.and.time[[ "results.across.regions.by.time" ]][[ from.region ]][[ to.region ]][[ the.time ]][[ "evaluated.results" ]] )
                  return( ...rv );
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
        sapply( times, function ( the.time ) {
          .bounds.types <- names( results.table.across.regions.by.time.and.bounds.type[[ from.region ]][[ to.region ]][[ the.time ]] );
          ....result.ignored <- 
            sapply( .bounds.types, function ( the.bounds.type ) {
              out.file <- paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", from.region, "_and_", to.region, "_", the.time, "_", the.bounds.type, "_evaluateIsMultiple.tab", sep = "" );
              ## TODO: REMOVE
              print( the.bounds.type );
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

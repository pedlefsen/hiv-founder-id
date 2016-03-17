### Read results tables.
## setting the.time to "1m.6m" will return pooled results over those times (but note current limitation that the bounds will be for either 1m (5weeks) or 6m (30weeks) but presently not use-the-right-bounds-matching-the-data-timepoints).
getFilteredResultsTables <- function (
    out.tab.file.suffix, the.region, the.time, the.bounds.type = "unbounded", to.region = NULL, results.dirname = "raw_edited_20160216", zeroNAs = FALSE, sort.column = "rmse" ) {
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
        colnames( results.zeroNAs ) <-
            gsub( "\\.zeroNAs$", "", colnames( results.zeroNAs ) );
    } else {
        results <- results.in[ , grep( "zeroNAs$", colnames( results.in ), invert = TRUE ), drop = FALSE ];
    }

    # Filter out anything that's impossible -- so that's all the deterministic bounds as well as any non-matching time.  Also we use 30weeks never 20weeks for the 6 month time point.
    if( the.time == "1m" ) {
        the.times.it.aint <- c( "20weeks", "30weeks", "1m5weeks_6m30weeks" );
    } else if( the.time == "1m6m" ) {
        the.times.it.aint <- c( "20weeks", "1m5weeks_6m30weeks" );
    } else if( the.time == "6m" ) {
        the.times.it.aint <- c( "5weeks", "20weeks", "1m5weeks_6m30weeks" );
    } else if( the.time == "1m.6m" ) {
        the.times.it.aint <- c( "\\.5weeks", "\\.30weeks", "20weeks", "1m5weeks_6m30weeks" );
    }
    results.filtered <- results[ grep( paste( c( "deterministic", the.times.it.aint ), collapse = "|" ), rownames( results ), invert = TRUE ), , drop = FALSE ];
    
    if( !is.null( sort.column ) && ( length( sort.column ) > 0 ) && !is.na( sort.column ) ) {
        stopifnot( length( sort.column ) == 1 );
        stopifnot( sort.column %in% names( results.filtered ) );
        results.filtered.sorted <-
            results.filtered[ order( results.filtered[[ sort.column ]] ), , drop = FALSE ];
        return( results.filtered.sorted );
    }
    
    return( results.filtered );
} # getFilteredResultsTables (..)

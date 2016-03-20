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
## setting the.time to "1m.6m" will return pooled results over those times (but note current limitation that the bounds will be for either 1m (5weeks) or 6m (30weeks) but presently not use-the-right-bounds-matching-the-data-timepoints).
getFilteredResultsTables <- function (
    out.tab.file.suffix, the.region, the.time, the.bounds.type = "unbounded", to.region = NULL, results.dirname = "raw_edited_20160216", zeroNAs = FALSE, sort.column = "rmse", column.pattern = NA, rowname.pattern.map = list( "\\.(days|time)\\.est" = "", "\\.mut\\.rate\\.coef" = "", "multifounder\\." = "(w/in clusts) ", "Synonymous\\." = "(syn) ", "is\\.poisson" = "FITS", "is\\.starlike" = "star-like", "is.one.founder" = "single-founder", "\\." = " " )
) {
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
        the.times.it.aint <- c( "\\.5weeks", "\\.30weeks", "20weeks" );
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

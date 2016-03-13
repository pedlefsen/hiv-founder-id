# If partition.size is not NA, then the results will be subsetted to just partitions of that size, meaning fasta files ending in _p[partition.size]_[\d]+.fasta
readIdentifyFounders <- function ( identify.founders.tab.filename, partition.size = NA ) {
             results.in <- read.delim( identify.founders.tab.filename, sep = "\t" );

             results <- as.matrix( results.in );
             suppressWarnings( mode( results ) <- "numeric" );

             # Exclude the ill-fated "multi-region" results, for now.
             # multiregion.results <- results[ is.na( results.in[ , 2 ] ), , drop = FALSE ];
             # rownames( multiregion.results ) <- gsub( ".*caprisa002_(\\d+)_.*", "\\1", gsub( ".*rv217_(\\d+)_.*", "\\1", as.character( results.in[ is.na( results.in[ , 2 ] ), 1 ] ) ) );
             results <- results[ !is.na( results.in[ , 2 ] ), , drop = FALSE ];
             rownames( results ) <- gsub( ".*caprisa002_(\\d+)_.*", "\\1", gsub( ".*rv217_(\\d+)_.*", "\\1", as.character( results.in[ !is.na( results.in[ , 2 ] ), 1 ] ) ) );

             if( !is.na( partition.size ) ) {
                 results <- results[ grep( paste( "_p", partition.size, "_\\d+\\.fasta$", sep = "" ), results.in[ , 1 ] ), , drop = FALSE ];
                 .in.file.name <- as.character( unlist( results.in[ grep( paste( "_p", partition.size, "_\\d+\\.fasta$", sep = "" ), results.in[ , 1 ] ), 1, drop = FALSE ] ) );
                 results.partition.id <- as.numeric( gsub( paste( "^.+_p", partition.size, "_(\\d+)\\.fasta$", sep = "" ), "\\1", .in.file.name ) );
                 # Replace the unused "infile" and "file" columns with a numeric "partition.size" column.
                 results <-
                     cbind( partition.id = results.partition.id, results[ , 3:ncol( results ) ] );
             }

             lambda.colnames <- grep( "lambda", colnames( results ), value = T );
             lambda <- results[ , lambda.colnames, drop = FALSE ];
             lambda.colnames.nb <- gsub( "[^\\.]+\\.Star[Pp]hy", "PFitter", gsub( "lambda.*$", "nbases", lambda.colnames ) );
             lambda.nb <- results[ , lambda.colnames.nb, drop = FALSE ];
             
             days.colnames <- c( grep( "time", colnames( results ), value = T ), grep( "days", colnames( results ), value = T ) );
             days <- results[ , days.colnames, drop = FALSE ];
             days.colnames.nb <- gsub( "[^\\.]+\\.Star[Pp]hy", "PFitter", gsub( "days.*$", "nbases", days.colnames ) );
             days.nb <- results[ , days.colnames.nb, drop = FALSE ];

             single.colnames <- grep( "\\.is\\.|fits", colnames( results ), perl = TRUE, value = TRUE );
             single.exceptInSites.colnames <- grep( "InSites", single.colnames, value = T, invert = T );
             single.exceptInSites <- results[ , single.exceptInSites.colnames, drop = FALSE ];
             single.exceptInSites.colnames.nb <-
                 gsub( "^Star[Pp]hy", "PFitter", gsub( "[^\\.]+\\.Star[Pp]hy", "PFitter", gsub( "\\....er\\.PFitter\\.", "\\.", gsub( "(?:is\\.|fits).*$", "nbases", single.exceptInSites.colnames ) ) ) );
             ## Special case: the one called "StarPhy.is.one.founder" is actually using the synonymous PFitter results, so its "nbases" should be the synonymous one.
             single.exceptInSites.colnames.nb[ single.exceptInSites.colnames == "StarPhy.is.one.founder" ] <- "Synonymous.PFitter.nbases";
             single.exceptInSites.nb <- results[ , single.exceptInSites.colnames.nb, drop = FALSE ];
             
             ## Results for which nbases == 0 [should] reflect no variation in the input sequences.  Use lambda estimate 0, days est 0, and isSingle 1 (TRUE).
             results[ , lambda.colnames ] <-
                 t( apply( results, 1, function( .row ) {
                     ifelse( .row[ lambda.colnames.nb ] == 0, 0, .row[ lambda.colnames ] )
                 } ) );
             results[ , days.colnames ] <-
                 t( apply( results, 1, function( .row ) {
                     ifelse( .row[ days.colnames.nb ] == 0, 0, .row[ days.colnames ] )
                 } ) );
             results[ , single.exceptInSites.colnames ] <-
                 t( apply( results, 1, function( .row ) {
                     ifelse( .row[ single.exceptInSites.colnames.nb ] == 0, 1, .row[ single.exceptInSites.colnames ] )
                 } ) );
             
             return( results );
} # readIdentifyFounders (..)

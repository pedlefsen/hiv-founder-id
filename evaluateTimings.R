
## First do all the stuff in README.postprocessing.txt.

source( "readIdentifyFounders_safetosource.R" );
source( "getDaysSinceInfection_safetosource.R" );

use.infer <- TRUE;
use.anchre <- TRUE;
use.partitions <- TRUE;
results.dirname <- "raw_edited_20160216";
#results.dirname <- "raw";

force.recomputation <- FALSE;
timings.results.by.region.and.time.Rda.filename <-
    paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/timings.results.by.region.and.time.Rda", sep = "" );

identify.founders.date.estimates <-
    c( "PFitter.time.est", "Synonymous.PFitter.time.est", "multifounder.PFitter.time.est", "multifounder.Synonymous.PFitter.time.est" ); 

rv217.gold.standard.infection.dates.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/rv217_gold_standard_timings.csv" );
rv217.gold.standard.infection.dates <- as.Date( as.character( rv217.gold.standard.infection.dates.in[,2] ), "%m/%d/%y" );
names( rv217.gold.standard.infection.dates ) <- as.character( rv217.gold.standard.infection.dates.in[,1] );

caprisa002.gold.standard.infection.dates.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/caprisa_002_gold_standard_timings.csv" );
caprisa002.gold.standard.infection.dates <- as.Date( as.character( caprisa002.gold.standard.infection.dates.in[,2] ), "%Y/%m/%d" );
names( caprisa002.gold.standard.infection.dates ) <- as.character( caprisa002.gold.standard.infection.dates.in[,1] );

regions <- c( "nflg", "v3", "rv217_v3" );
times <- c( "1m", "6m", "1m6m" );
#times <- c( "1m6m" );

rmse <- function( x, na.rm = FALSE ) {
    if( na.rm ) {
        x <- x[ !is.na( x ) ];
    }
    return( sqrt( mean( x ** 2 ) ) );
}

## This is from PFitter.  Epsilon is the per position mutation rate, per generation -- but I think that the generations / day is fixed at 2, so really this is the mutation rate position per half-day.
default.epsilon <- 2.16e-05;
phi <- sqrt(1+4/3);
daysFromLambda <- function ( lambda, nb, epsilon = default.epsilon ) {
    1.5 * ((phi)/(1+phi)) * (lambda/(epsilon*nb) - (1-phi)/(phi^2) )
}

compute.results.one.per.ppt <- function ( results, weights ) {
    apply( results, 2, function ( .column ) {
        .rv <- 
        sapply( unique( rownames( results ) ), function( .ppt ) {
            .ppt.cells <- .column[ rownames( results ) == .ppt ];
            .ppt.weights <- weights[ rownames( results ) == .ppt ];
            .ppt.weights <- .ppt.weights / sum( .ppt.weights, na.rm = TRUE );
            sum( .ppt.cells * .ppt.weights, na.rm = TRUE );
        } );
        names( .rv ) <- unique( rownames( results ) );
        return( .rv );
    } );
} # compute.results.one.per.ppt

compute.diffs.by.stat <- function ( results.one.per.ppt ) {
    diffs.by.stat <- 
        lapply( colnames( results.one.per.ppt ), function( .stat ) {
            .rv <- ( as.numeric( results.one.per.ppt[ , .stat ] ) - as.numeric( days.since.infection[ rownames( results.one.per.ppt ) ] ) );
            names( .rv ) <- rownames( results.one.per.ppt );
            return( .rv );
        } );
    names( diffs.by.stat ) <- colnames( results.one.per.ppt );
    return( diffs.by.stat );
} # compute.diffs.by.stat ( results.one.per.ppt )

getTimingsResultsByRegionAndTime <- function ( partition.size = NA ) {
    if( !is.na( partition.size ) ) {
        regions <- "v3"; # Only v3 has partition results at this time.
    }
    timings.results.by.region.and.time <-
        lapply( regions, function( the.region ) {
            ## TODO: REMOVE
            cat( the.region, fill = T );
       timings.results.by.time <- 
           lapply( times, function( the.time ) {
               ## TODO: REMOVE
               cat( the.time, fill = T );

               .days.since.infection.filename <- paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/sampleDates.tbl", sep = "" );
               if( the.region == "v3" ) {
                   days.since.infection <-
                       getDaysSinceInfection(
                           .days.since.infection.filename,
                           caprisa002.gold.standard.infection.dates
                       );
               } else {
                   stopifnot( ( the.region == "nflg" ) || ( length( grep( "rv217", the.region ) ) > 0 ) );
                   days.since.infection <-
                       getDaysSinceInfection(
                           .days.since.infection.filename,
                           rv217.gold.standard.infection.dates
                       );
               }
               days.since.infection <- sapply( 1:nrow( sample.dates.in ), function( .i ) { as.numeric( as.Date( as.character( sample.dates.in[ .i, 2 ] ) ) - ifelse( the.region == "v3", caprisa002.gold.standard.infection.dates[ as.character( sample.dates.in[ .i, 1 ] ) ], rv217.gold.standard.infection.dates[ as.character( sample.dates.in[ .i, 1 ] ) ] ) ) } );
               names( days.since.infection ) <- sample.dates.in[ , "ptid" ];
                   
               ## identify-founders results
               if( is.na( partition.size ) ) {
                   results <- readIdentifyFounders( paste( paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/identify_founders.tab", sep = "" ) ) );
               } else {
                   results <- readIdentifyFounders( paste( paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/partitions/identify_founders.tab", sep = "" ) ), partition.size = partition.size );
               }
               days.colnames <- c( grep( "time", colnames( results ), value = T ), grep( "days", colnames( results ), value = T ) );
  
               ### TODO: HERE IS WHERE I CAN DO SOME PLAYING AROUND WITH RECALIBRATION.
               ## The current problem is that I can't reproduce the days perfectly -- in some cases the recomputed days are way off, so something is wrong with the identify_founders table info about the synonymous pfitter, perhaps.
               days.est.colnames <- grep( "est", days.colnames, value = TRUE );
               days.est <- results[ , days.est.colnames, drop = FALSE ];
               lambda.est.colnames <-
                   gsub( "PFitter\\.lambda\\.est", "PFitter.lambda", gsub( "(?:days|time|fits)", "lambda", days.est.colnames, perl = TRUE ) );
               stopifnot( all( lambda.est.colnames %in% colnames( results ) ) );
               days.est.colnames.nb <- gsub( "[^\\.]+\\.Star[Pp]hy", "PFitter", gsub( "(?:days|time|fits).*$", "nbases", days.est.colnames, perl = TRUE ) );
               days.est.nb <- results[ , days.est.colnames.nb, drop = FALSE ];
               days.est.colnames.nseq <- gsub( "[^\\.]+\\.Star[Pp]hy", "PFitter", gsub( "(?:days|time|fits).*$", "nseq", days.est.colnames, perl = TRUE ) );
               days.est.nseq <- results[ , days.est.colnames.nseq, drop = FALSE ];
               
               ## proof of concept:
               create.function.to.be.optimized <- # defaults to use rmse aka sqrt( sum( bias^2 ) ), but you can change that to abs( bias ).
                   function ( days.est.col.i, for.bias = FALSE, weights = days.est.nseq*days.est.nb ) {

                       return( function( log.epsilon, return.days.est.new = FALSE ) {
                       .days.est.new <- 
                         apply( results, 1, function( .row ) {
                             if( is.na( .row[ days.est.colnames.nb[ days.est.col.i ] ] ) || .row[ days.est.colnames.nb[ days.est.col.i ] ] == 0 ) {
                                 return( NA );
                             }
                             return( daysFromLambda( .row[ lambda.est.colnames[ days.est.col.i ] ], .row[ days.est.colnames.nb[ days.est.col.i ] ], epsilon = exp( log.epsilon ) ) );
                         } );
                       .mat <- matrix( .days.est.new, ncol = 1 );
                       colnames( .mat ) <- days.est.colnames[ days.est.col.i ];
                       rownames( .mat ) <- rownames( days.est );
                       if( return.days.est.new ) {
                           return( .mat );
                       }
                       .results.one.per.ppt <-
                           compute.results.one.per.ppt( .mat, weights );
                       .diffs.by.stat <- compute.diffs.by.stat( .results.one.per.ppt );
                       if( for.bias ) {
                         .bias <- mean( .diffs.by.stat[[1]], na.rm = TRUE );
                         return( abs( .bias ) );
                       } else {
                           .rmse <- rmse( .diffs.by.stat[[1]], na.rm = TRUE );
                           return( .rmse );
                       }
                     } );
                   } # create.function.to.be.optimized (..)
               optimal.log.epsilons.and.rmses <- sapply( 1:length( days.est.colnames ), function ( .col.i ) {
                   function.to.be.optimized <- create.function.to.be.optimized( .col.i );
                   optimize( function.to.be.optimized, interval = c( log( default.epsilon ) - log( 1000 ), log( default.epsilon ) + log( 1000 ) ) );
               } );
               colnames( optimal.log.epsilons.and.rmses ) <- days.est.colnames;
               rownames( optimal.log.epsilons.and.rmses ) <- c( "log.epsilon", "rmse" );
               optimal.epsilons <-
                   exp( unlist( optimal.log.epsilons.and.rmses[ "log.epsilon", ] ) );
               optimized.results <- 
                   sapply( 1:length( days.est.colnames ), function ( .col.i ) {
                   .optimal.log.epsilon <- optimal.log.epsilons.and.rmses[[ "log.epsilon", .col.i ]];
                   .new.column <- create.function.to.be.optimized( .col.i )( .optimal.log.epsilon, return.days.est.new = TRUE );
                   return( .new.column );
               } );
               colnames( optimized.results ) <-
                   paste( "optimized", days.est.colnames, sep = "." );
               rownames( optimized.results ) <-
                   rownames( results );

               ## Also do it for for.bias results (objective fn to minimize is abs( bias ))
               optimal.log.epsilons.and.biases <- sapply( 1:length( days.est.colnames ), function ( .col.i ) {
                   function.to.be.optimized.for.bias <- create.function.to.be.optimized( .col.i, for.bias = TRUE );

                   optimize( function.to.be.optimized.for.bias, interval = c( log( default.epsilon ) - log( 1000 ), log( default.epsilon ) + log( 1000 ) ) );
               } );
               colnames( optimal.log.epsilons.and.biases ) <- days.est.colnames;
               rownames( optimal.log.epsilons.and.biases ) <- c( "log.epsilon", "bias" );
               optimal.for.bias.epsilons <-
                   exp( unlist( optimal.log.epsilons.and.biases[ "log.epsilon", ] ) );
               optimized.for.bias.results <- 
                   sapply( 1:length( days.est.colnames ), function ( .col.i ) {
                   .optimal.for.bias.log.epsilon <- optimal.log.epsilons.and.biases[[ "log.epsilon", .col.i ]];
                   .new.column <- create.function.to.be.optimized( .col.i, for.bias = TRUE )( .optimal.for.bias.log.epsilon, return.days.est.new = TRUE );
                   return( .new.column );
               } );
               colnames( optimized.for.bias.results ) <-
                   paste( "optimized.for.bias", days.est.colnames, sep = "." );
               rownames( optimized.for.bias.results ) <-
                   rownames( results );
               
               results <-
                   cbind(
                       results[ , days.est.colnames, drop = FALSE ],
                       optimized.results,
                       optimized.for.bias.results
                   );
               
               if( use.infer && is.na( partition.size ) ) { ## TODO: Handle the infer results for the partitions
                 ## Add to results: "infer" results.
                 if( ( the.region == "v3" ) || ( the.region == "rv217_v3" ) ) {
                     the.region.dir <- "v3_edited_20160216";
                 } else {
                     the.region.dir <- "nflg_copy_20160222";
                 }
                 infer.results.directories <- dir( paste( "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/", results.dirname, "/", the.region.dir, "/", the.time, sep = "" ), "founder-inference-bakeoff_", full.name = TRUE );
                 # Special: for v3, separate out the caprisa seqs from the rv217 seqs
                 if( the.region == "v3" ) {
                      infer.results.directories <- grep( "_100\\d\\d\\d$", infer.results.directories, value = TRUE );
                 } else if( the.region == "rv217_v3" ) {
                      infer.results.directories <- grep( "_100\\d\\d\\d$", infer.results.directories, value = TRUE, invert = TRUE );
                 }
                   
                 infer.results.files <- sapply( infer.results.directories, dir, "outtoi.csv", full.name = TRUE );
                 infer.results.list <-
                     lapply( unlist( infer.results.files ), function( .file ) {
                         .rv <- as.matrix( read.csv( .file, header = FALSE ), nrow = 1 );
                         stopifnot( ncol( .rv ) == 3 );
                         return( .rv );
                     } );
                 names( infer.results.list ) <- unlist( infer.results.files );
  
                 if( length( infer.results.list ) > 0 ) {
                     infer.results <- do.call( rbind, infer.results.list );
                     colnames( infer.results ) <- c( "Infer", "Infer.CI.low", "Infer.CI.high" );
                     .tmp <- gsub( "^.+_(\\d+)/.+$", "\\1", names( infer.results.list ) );
                     rownames( infer.results ) <- .tmp;
                     
                     # Add just the estimate from infer.
                     new.results.columns <- matrix( NA, nrow = nrow( results ), ncol = 1 );
                     rownames( new.results.columns ) <- rownames( results );
                     colnames( new.results.columns ) <- "Infer.time.est";
                     .shared.ptids <-
                         intersect( rownames( new.results.columns ), rownames( infer.results ) );
                     .result.ignored <- sapply( .shared.ptids, function( .ptid ) {
                         .infer.subtable <-
                             infer.results[ rownames( infer.results ) == .ptid, 1, drop = FALSE ];
                         stopifnot( sum( rownames( infer.results ) == .ptid ) == nrow( .infer.subtable ) );
                         new.results.columns[ rownames( infer.results ) == .ptid, ] <<-
                             .infer.subtable;
                         return( NULL );
                     } );
                       
                     results <- cbind( results, new.results.columns );
                 }
                   
               } # End if use.infer
  
#                if( use.anchre && ( the.time == "1m6m" ) && is.na( partition.size ) ) {  ## TODO: Handle the anchre results for the partitions
#                  ## Add to results: "anchre" results. (only at 1m6m)
#                  anchre.results.directories <- dir( paste( "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/", results.dirname, "/", the.region.dir, "/1m6m", sep = "" ), "anchre", full.name = TRUE );
#                  if( length( anchre.results.directories ) > 0 ) {
#                    anchre.results.files <-
#                        sapply( anchre.results.directories, dir, "mrca.csv", full.name = TRUE );
#                    anchre.results <- do.call( rbind,
#                        lapply( unlist( anchre.results.files ), function( .file ) {
#                            stopifnot( file.exists( .file ) );
#                            .file.short <-
#                                gsub( "^.*?\\/?([^\\/]+?)$", "\\1", .file, perl = TRUE );
#                            .file.short.nosuffix <-
#                                gsub( "^([^\\.]+)(\\..+)?$", "\\1", .file.short, perl = TRUE );
#                            .file.converted <-
#                                paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region.dir, "/1m6m/", .file.short.nosuffix, ".anc2tsv.tab", sep = "" );
#                            # convert it.
#                            system( paste( "./anc2tsv.sh", .file, ">", .file.converted ) );
#                            stopifnot( file.exists( .file.converted ) );
#                            .rv <- as.matrix( read.delim( .file.converted, header = TRUE, sep = "\t" ), nrow = 1 );
#                            ## No negative dates!  Just call it NA.
#                            .rv <- apply( .rv, 1:2, function( .str ) { if( length( grep( "^-", .str ) ) > 0 ) { NA } else { .str } } );
#                            stopifnot( ncol( .rv ) == 4 );
#                            return( .rv );
#                        } ) );
#                    colnames( anchre.results ) <- c( "Anchre.r2t.est", "Anchre.est", "Anchre.CI.low", "Anchre.CI.high" );
#                    rownames( anchre.results ) <-
#                        gsub( "^.+_(\\d+)$", "\\1", names( unlist( anchre.results.files ) ) );
#                    # Special: for v3, only use caprisa seqs (not rv217, for now).
#                    if( the.region == "v3" ) {
#                        anchre.results <-
#                            anchre.results[ grep( "^100\\d\\d\\d", rownames( anchre.results ) ), , drop = FALSE ];
#                    } else if( the.region == "rv217_v3" ) {
#                        anchre.results <-
#                            anchre.results[ grep( "^100\\d\\d\\d", rownames( anchre.results ), invert = TRUE ), , drop = FALSE ];
#                    }
#                    # Add just the estimate from anchre.
#                    sample.dates <- as.Date( as.character( sample.dates.in[ , 2 ] ) );
#                    names( sample.dates ) <- sample.dates.in[ , 1 ];
#                    anchre.r2t.days.before.sample <- sapply( 1:nrow( anchre.results ), function( .i ) { 0 - as.numeric( as.Date( anchre.results[ .i, 1 ] ) - sample.dates[ rownames( anchre.results )[ .i ] ] ) } );
#                    names( anchre.r2t.days.before.sample ) <- rownames( anchre.results );
#                    anchre.days.before.sample <- sapply( 1:nrow( anchre.results ), function( .i ) { 0 - as.numeric( as.Date( anchre.results[ .i, 2 ] ) - sample.dates[ rownames( anchre.results )[ .i ] ] ) } );
#                    names( anchre.days.before.sample ) <- rownames( anchre.results );
#       
#                    results <- cbind( results, anchre.r2t.days.before.sample, anchre.days.before.sample );
#                    colnames( results )[ ncol( results ) - 1:0 ] <- c( "Anchre.r2t.time.est", "Anchre.bst.time.est" ) ;
#                  } # End if there's any anchre results.
#                } # End if use.ancher and the.time is 1m6m, add anchre results too.

               if( is.na( partition.size ) ) {
                 ## Now the issue is that there are multiple input files per ppt, eg for the NFLGs ther are often "LH" and "RH" files.  What to do?  The number of sequences varies.  Do a weighted average.
                 .weights <- days.est.nseq*days.est.nb;
                 results.one.per.ppt <-
                     compute.results.one.per.ppt( results, .weights );
                 diffs.by.stat <- compute.diffs.by.stat( results.one.per.ppt );
    
                 return( list( bias = lapply( diffs.by.stat, mean, na.rm = T ), se = lapply( diffs.by.stat, sd, na.rm = T ), rmse = lapply( diffs.by.stat, rmse, na.rm = T ), n = lapply( diffs.by.stat, function( .vec ) { sum( !is.na( .vec ) ) } ), epsilon = c( rep( default.epsilon, length( days.est.colnames ) ), optimal.epsilons, optimal.for.bias.epsilons ) ) );
               } else {
                   ## Here the multiple results per participant come from the partitions.  We want to evaluate each one, and summarize them afterwards.
                   results.per.ppt <- apply( results, 2, function ( .column ) {
                     .rv <- 
                       lapply( unique( rownames( results ) ), function( .ppt ) {
                           .values <- .column[ rownames( results ) == .ppt ];
                           .diffs <- ( as.numeric( .values ) - as.numeric( days.since.infection[ .ppt ] ) );
                           return( list( values = .values, diffs = .diffs ) );
                     } );
                     names( .rv ) <- unique( rownames( results ) );
                     return( .rv );
                   } );
                   diffs.per.ppt <- lapply( results.per.ppt, function( .lst ) { lapply( .lst, function( ..lst ) { ..lst[[ "diffs" ]] } ) } );
                   mean.bias.per.ppt <- lapply( diffs.per.ppt, function( .lst ) { lapply( .lst, mean, na.rm = T ) } );
                   sd.bias.per.ppt <- lapply( diffs.per.ppt, function( .lst ) { lapply( .lst, sd, na.rm = T ) } );
                   median.mean.bias <- lapply( mean.bias.per.ppt, function( .lst ) { median( unlist( .lst ), na.rm = T ) } );
                   min.mean.bias <- lapply( mean.bias.per.ppt, function( .lst ) { min( unlist( .lst ), na.rm = T ) } );
                   max.mean.bias <- lapply( mean.bias.per.ppt, function( .lst ) { max( unlist( .lst ), na.rm = T ) } );
                   median.sd.bias <- lapply( sd.bias.per.ppt, function( .lst ) { median( unlist( .lst ), na.rm = T ) } );
                   min.sd.bias <- lapply( sd.bias.per.ppt, function( .lst ) { min( unlist( .lst ), na.rm = T ) } );
                   max.sd.bias <- lapply( sd.bias.per.ppt, function( .lst ) { max( unlist( .lst ), na.rm = T ) } );
                   return( list( median.mean.bias = median.mean.bias, min.mean.bias = min.mean.bias, max.mean.bias = max.mean.bias,  median.sd.bias = median.sd.bias, min.sd.bias = min.sd.bias, max.sd.bias = max.sd.bias ) );
               }
           } ); # End foreach the.time
       names( timings.results.by.time ) <- times;
       return( timings.results.by.time );
   } ); # End foreach the.region
  names( timings.results.by.region.and.time ) <- regions;
  return( timings.results.by.region.and.time );
} # getTimingsResultsByRegionAndTime ( partition.size )

if( force.recomputation || !file.exists( timings.results.by.region.and.time.Rda.filename ) ) {
    timings.results.by.region.and.time <- getTimingsResultsByRegionAndTime();
    save( timings.results.by.region.and.time, file = timings.results.by.region.and.time.Rda.filename );
} else {
    # loads timings.results.by.region.and.time
    load( file = timings.results.by.region.and.time.Rda.filename );
}

# Make a table out of it. (one per study).
results.table.by.region.and.time <-
    lapply( names( timings.results.by.region.and.time ), function( the.region ) {
        .rv <- 
            lapply( names( timings.results.by.region.and.time[[ the.region ]] ), function( the.time ) {
                sapply( timings.results.by.region.and.time[[ the.region ]][[ the.time ]], function( results.list ) { results.list } ) } );
        names( .rv ) <- times;
        return( .rv );
    } );
names( results.table.by.region.and.time ) <- regions;

## Write these out.
.result.ignored <- sapply( regions, function ( the.region ) {
    ..result.ignored <- 
    sapply( times, function ( the.time ) {
        out.file <- paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/evaluateTimings.tab", sep = "" );
        .tbl <-
            apply( results.table.by.region.and.time[[ the.region ]][[ the.time ]], 1:2, function( .x ) { sprintf( "%0.2f", .x ) } );
        write.table( .tbl, quote = FALSE, file = out.file, sep = "\t" );
        return( NULL );
    } );
    return( NULL );
} );

## For partition size == 10
timings.results.by.region.and.time.p10 <- getTimingsResultsByRegionAndTime( partition.size = 10 );

# Make a table out of it. (one per study).
results.table.by.region.and.time.p10 <-
    lapply( names( timings.results.by.region.and.time.p10 ), function( the.region ) {
        .rv <- 
            lapply( names( timings.results.by.region.and.time.p10[[ the.region ]] ), function( the.time ) {
                sapply( timings.results.by.region.and.time.p10[[ the.region ]][[ the.time ]], function( results.list ) { results.list } ) } );
        names( .rv ) <- times;
        return( .rv );
    } );
names( results.table.by.region.and.time.p10 ) <- "v3";

## Write these out.
.result.ignored <- sapply( "v3", function ( the.region ) {
    ..result.ignored <- 
    sapply( times, function ( the.time ) {
        out.file <- paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/evaluateTimings_p10.tab", sep = "" );
        .tbl <- apply( results.table.by.region.and.time.p10[[ the.region ]][[ the.time ]], 1:2, function( .x ) { sprintf( "%0.2f", .x ) } );
        write.table( .tbl, quote = FALSE, file = out.file, sep = "\t" );
        return( NULL );
    } );
    return( NULL );
} );

######################
# ( results.table.by.region.and.time.p10 )
# $v3
# $v3$`1m`
#                                          median.mean.bias min.mean.bias
# PFitter.time.est                         -3.128968        -41.66667    
# Synonymous.PFitter.time.est              -37.07143        -46.35714    
# multifounder.PFitter.time.est            -7.85989         -43.76471    
# multifounder.Synonymous.PFitter.time.est -35.96273        -45.33333    
#                                          max.mean.bias median.sd.bias
# PFitter.time.est                         1266.304      30.76341      
# Synonymous.PFitter.time.est              114.6522      14.98793      
# multifounder.PFitter.time.est            204.7143      30.05146      
# multifounder.Synonymous.PFitter.time.est 26.85714      12.76584      
#                                          min.sd.bias max.sd.bias
# PFitter.time.est                         8.253342    257.8342   
# Synonymous.PFitter.time.est              0           44.67235   
# multifounder.PFitter.time.est            4.176263    128.7112   
# multifounder.Synonymous.PFitter.time.est 0           61.42207   
# 
# $v3$`6m`
#                                          median.mean.bias min.mean.bias
# PFitter.time.est                         63.74904         -115.8333    
# Synonymous.PFitter.time.est              -123.781         -171         
# multifounder.PFitter.time.est            -22.55556        -143.6667    
# multifounder.Synonymous.PFitter.time.est -110.6875        -147.3333    
#                                          max.mean.bias median.sd.bias
# PFitter.time.est                         872.7         51.16241      
# Synonymous.PFitter.time.est              19.2          20.39825      
# multifounder.PFitter.time.est            243.3103      60.83165      
# multifounder.Synonymous.PFitter.time.est -61.33333     37.16665      
#                                          min.sd.bias max.sd.bias
# PFitter.time.est                         20.57507    174.5478   
# Synonymous.PFitter.time.est              0           39.92873   
# multifounder.PFitter.time.est            19.82811    172.6328   
# multifounder.Synonymous.PFitter.time.est 8.485281    77.33471   
# 
# $v3$`1m6m`
#                                          median.mean.bias min.mean.bias
# PFitter.time.est                         148.8333         -12.90909    
# Synonymous.PFitter.time.est              -27.71429        -56.33333    
# multifounder.PFitter.time.est            28.77778         -38.11538    
# multifounder.Synonymous.PFitter.time.est -18.16667        -55.25       
#                                          max.mean.bias median.sd.bias
# PFitter.time.est                         1280.304      31.00045      
# Synonymous.PFitter.time.est              141.2609      13.99429      
# multifounder.PFitter.time.est            221.2759      28.65699      
# multifounder.Synonymous.PFitter.time.est 31.6087       25.86964      
#                                          min.sd.bias max.sd.bias
# PFitter.time.est                         18.43685    132.4031   
# Synonymous.PFitter.time.est              4.618802    35.29038   
# multifounder.PFitter.time.est            15.96975    109.4992   
# multifounder.Synonymous.PFitter.time.est 4.349329    51.98283   


###########################################################################
# results.table.by.region.and.time
### Raw run, without recombination detection/removal nor hypermutation detection/removal:
# $nflg
# $nflg$`1m`
#                                          bias      se       rmse    
# PFitter.time.est                         33.853    102.5296 106.613 
# Synonymous.PFitter.time.est              -24.64259 24.75365 34.684  
# multifounder.PFitter.time.est            4.501963  54.06851 53.50202
# multifounder.Synonymous.PFitter.time.est -31.84253 12.5926  34.1777 
# Infer.time.est                           146.0784  112.2472 183.2712
# 
# $nflg$`6m`
#                                          bias      se       rmse    
# PFitter.time.est                         -23.03854 101.3019 102.4257
# Synonymous.PFitter.time.est              -141.2521 27.61474 143.8482
# multifounder.PFitter.time.est            -65.98357 74.86203 98.96122
# multifounder.Synonymous.PFitter.time.est -150.5546 24.89646 152.5395
# Infer.time.est                           137.9264  177.8726 223.0063
# 
# $nflg$`1m6m`
#                                          bias     se       rmse    
# PFitter.time.est                         90.79881 108.7833 139.947 
# Synonymous.PFitter.time.est              -9.74849 29.24275 30.24138
# multifounder.PFitter.time.est            27.89046 63.01514 67.7003 
# multifounder.Synonymous.PFitter.time.est -23.7626 19.66338 30.58103
# Infer.time.est                           335.0502 476.9423 574.6778
# 
# 
# $v3
# $v3$`1m`
#                                          bias   se       rmse    
# PFitter.time.est                         164.65 340.783  370.7237
# Synonymous.PFitter.time.est              -23    58.29869 61.3009 
# multifounder.PFitter.time.est            12     111.4757 109.3138
# multifounder.Synonymous.PFitter.time.est -41.3  31.20746 51.2923 
# Infer.time.est                           413.3  173.6164 446.6009
# 
# $v3$`6m`
#                                          bias      se       rmse    
# PFitter.time.est                         161.6667  253.3017 294.5052
# Synonymous.PFitter.time.est              -119.0556 76.703   140.4661
# multifounder.PFitter.time.est            30.22222  262.8089 257.1863
# multifounder.Synonymous.PFitter.time.est -141.8333 60.15592 153.4092
# Infer.time.est                           361.9444  278.6665 452.0449
# 
# $v3$`1m6m`
#                                          bias      se       rmse    
# PFitter.time.est                         249.2941  355.0214 425.1748
# Synonymous.PFitter.time.est              -13.47059 66.84788 66.23621
# multifounder.PFitter.time.est            132.0588  352.2764 366.3854
# multifounder.Synonymous.PFitter.time.est -33.23529 46.61616 56.12329
# Infer.time.est                           464.6471  313.3397 555.2506
# Anchre.r2t.time.est                      656.1176  1477.022 1575.994
# Anchre.bst.time.est                      82935.71  107414.3 133182.1
# 
# 
# $rv217_v3
# $rv217_v3$`1m`
#                                          bias      se       rmse    
# PFitter.time.est                         41.55556  144.4919 148.4077
# Synonymous.PFitter.time.est              -28.27778 37.03096 46.18261
# multifounder.PFitter.time.est            -2.305556 68.9268  68.00184
# multifounder.Synonymous.PFitter.time.est -32.5     24.94623 40.75878
# Infer.time.est                           143.1389  205.5429 248.1191
# 
# $rv217_v3$`6m`
#                                          bias      se       rmse    
# PFitter.time.est                         40.09091  184.6078 186.1574
# Synonymous.PFitter.time.est              -125      59.86443 138.2033
# multifounder.PFitter.time.est            -64       77.77612 99.80891
# multifounder.Synonymous.PFitter.time.est -142.8182 30.8063  146.0045
# Infer.time.est                           247.3636  270.7396 363.686 
# 
# $rv217_v3$`1m6m`
#                                          bias      se       rmse    
# PFitter.time.est                         161       210.6419 262.5766
# Synonymous.PFitter.time.est              0.1515152 54.44844 53.61733
# multifounder.PFitter.time.est            30.33333  43.68328 52.63568
# multifounder.Synonymous.PFitter.time.est -16.75758 27.95424 32.22694
# Infer.time.est                           542.7273  545.3223 763.4906
# Anchre.r2t.time.est                      165.7273  568.4986 583.8344
# Anchre.bst.time.est                      527.7879  2035.306 2072.559





### Original run, using RAPBeta and my implementation of Hypermut 2.0.  Also the computation of these values was not correctly excluding multi-region results nor accounting for the sometimes multiple estimates per ppt.
# $nflg
# $nflg$`1m`
# $nflg$`1m`$rv217
#                                          bias      se       rmse    
# PFitter.time.est                         39.39216  115.1498 120.6285
# Synonymous.PFitter.time.est              -23.13725 28.43942 36.4455 
# multifounder.PFitter.time.est            6.137255  61.02033 60.73004
# multifounder.Synonymous.PFitter.time.est -30.84314 16.37116 34.84335
# 
# 
# $nflg$`6m`
# $nflg$`6m`$rv217
#                                          bias      se       rmse    
# PFitter.time.est                         -23.27778 124.5701 125.5874
# Synonymous.PFitter.time.est              -142.7407 31.42117 146.0956
# multifounder.PFitter.time.est            -66.92593 91.80863 112.9239
# multifounder.Synonymous.PFitter.time.est -151.3889 28.42296 153.9854
# 
# 
# $nflg$`1m6m`
# $nflg$`1m6m`$rv217
#                                          bias      se       rmse    
# PFitter.time.est                         92.3      130.6546 158.18  
# Synonymous.PFitter.time.est              -8.533333 34.35889 34.8425 
# multifounder.PFitter.time.est            21.83333  70.30431 72.48885
# multifounder.Synonymous.PFitter.time.est -25.23333 20.79072 32.47409
# 
# 
# 
# $v3
# $v3$`1m`
# $v3$`1m`$caprisa002
#                                          bias   se       rmse    
# PFitter.time.est                         154    341.5987 366.8395
# Synonymous.PFitter.time.est              -21.9  61.98209 64.25963
# multifounder.PFitter.time.est            -1     71.1248  69.33109
# multifounder.Synonymous.PFitter.time.est -41.05 31.57693 51.30643
# 
# $v3$`1m`$rv217
#                                          bias      se       rmse    
# PFitter.time.est                         30.94444  131.9155 133.7007
# Synonymous.PFitter.time.est              -30.16667 35.37998 46.11941
# multifounder.PFitter.time.est            -12.69444 32.68973 34.64222
# multifounder.Synonymous.PFitter.time.est -36.19444 15.79117 39.40142
# 
# 
# $v3$`6m`
# $v3$`6m`$caprisa002
#                                          bias      se       rmse    
# PFitter.time.est                         139.2778  203.4195 241.8242
# Synonymous.PFitter.time.est              -120.7778 73.17309 140.1575
# multifounder.PFitter.time.est            -25.88889 156.6712 154.4424
# multifounder.Synonymous.PFitter.time.est -144.7778 55.60528 154.5341
# 
# $v3$`6m`$rv217
#                                          bias      se       rmse    
# PFitter.time.est                         39.24242  183.8551 185.2521
# Synonymous.PFitter.time.est              -125.1515 59.594   138.2271
# multifounder.PFitter.time.est            -67.60606 75.71325 100.6447
# multifounder.Synonymous.PFitter.time.est -143.0909 30.37923 146.1846
# 
# 
# $v3$`1m6m`
# $v3$`1m6m`$caprisa002
#                                          bias      se       rmse    
# PFitter.time.est                         259.2105  339.6857 420.1232
# Synonymous.PFitter.time.est              6.210526  96.55947 94.18906
# multifounder.PFitter.time.est            112.8421  303.213  315.963 
# multifounder.Synonymous.PFitter.time.est -22.57895 69.40807 71.23017
# 
# $v3$`1m6m`$rv217
#                                          bias      se       rmse    
# PFitter.time.est                         131.9394  160.9126 206.1947
# Synonymous.PFitter.time.est              -5.333333 43.23603 42.90864
# multifounder.PFitter.time.est            27.81818  39.11877 47.51587
# multifounder.Synonymous.PFitter.time.est -16.78788 27.88229 32.18225

#### OLD, with Infer and anchre
# $nflg
# $nflg$`1m`
#                                          bias      rmse    
# PFitter.time.est                         39.39216  115.1498
# Synonymous.PFitter.time.est              -23.13725 28.43942
# multifounder.PFitter.time.est            6.137255  61.02033
# multifounder.Synonymous.PFitter.time.est -30.84314 16.37116
# Infer                                    86.92593  5.128547
# 
# $nflg$`6m`
#                                          bias      rmse    
# PFitter.time.est                         -23.27778 124.5701
# Synonymous.PFitter.time.est              -142.7407 31.42117
# multifounder.PFitter.time.est            -66.92593 91.80863
# multifounder.Synonymous.PFitter.time.est -151.3889 28.42296
# Infer                                    85.53333  13.15222
# 
# $nflg$`1m6m`
#                                          bias      rmse    
# PFitter.time.est                         92.3      130.6546
# Synonymous.PFitter.time.est              -8.533333 34.35889
# multifounder.PFitter.time.est            21.83333  70.30431
# multifounder.Synonymous.PFitter.time.est -25.23333 20.79072
# Infer                                    73.6875   25.60916
# Anchre.r2t                               53.4375   140.108 
# Anchre.bst                               291.9375  757.6369
# 
# 
# $v3
# $v3$`1m`
#                                          bias   rmse    
# PFitter.time.est                         154    341.5987
# Synonymous.PFitter.time.est              -21.9  61.98209
# multifounder.PFitter.time.est            -1     71.1248 
# multifounder.Synonymous.PFitter.time.est -41.05 31.57693
# Infer                                    77     6.24921 
# 
# $v3$`6m`
#                                          bias      rmse    
# PFitter.time.est                         139.2778  203.4195
# Synonymous.PFitter.time.est              -120.7778 73.17309
# multifounder.PFitter.time.est            -25.88889 156.6712
# multifounder.Synonymous.PFitter.time.est -144.7778 55.60528
# Infer                                    93.38889  12.71469
# 
# $v3$`1m6m`
#                                          bias      rmse    
# PFitter.time.est                         259.2105  339.6857
# Synonymous.PFitter.time.est              6.210526  96.55947
# multifounder.PFitter.time.est            112.8421  303.213 
# multifounder.Synonymous.PFitter.time.est -22.57895 69.40807
# Infer                                    88.84211  31.68817
# Anchre.r2t                               620.6842  1530.37 
# Anchre.bst                               325964    337789.5


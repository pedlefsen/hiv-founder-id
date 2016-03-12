## First do all the stuff in README.postprocessing.txt.

source( "readIdentifyFounders_safetosource.R" );
source( "getDaysSinceInfection_safetosource.R" );
source( "summarizeCovariatesOnePerParticipant_safetosource.R" );

#' Evaluate timings estimates and produce results tables.
#'
#' This function runs the BakeOff results analysis for the timings results.
#'
#' The "center of bounds" (COB) approach is the way that we do it at the
#' VTN, and the way it was done in RV144, etc: use the midpoint
#' between the bounds on the actual infection time computed from the
#' dates and results of the HIV positivity tests (antibody or PCR).
#' The typical approach is to perform antibody testing every X days
#' (historically this is 6 months in most HIV vaccine trials, except
#' during the vaccination phase there are more frequent visits and on
#' every visit HIV testing is conducted).  The (fake) bounds used here
#' are calculated in the createArtificialBoundsOnInfectionDate.R file.
#' The actual bounds would be too tight, since the participants were
#' detected HIV+ earlier in these people than what we expect to see in
#' a trial in which testing is conducted every X days.  For the center
#' of bounds approach we load the bounds files in subdirs of the
#' "bounds" subdirectory eg at
#' /fh/fast/edlefsen_p/bakeoff/analysis_sequences/bounds/nflg/1m/.
#' These files have names beginning with "artificialBounds_" and
#' ending with ".tab".
#'
#' @param use.center.of.bounds compute results for the COB approach, and also return evaluations of bounded versions of the other results.
#' @param use.infer compute results for the PREAST approach.
#' @param use.anchre compute results for the anchre approach (presently disabled).
#' @param results.dirname the subdirectory of "/fh/fast/edlefsen_p/bakeoff/analysis_sequences" and also of "/fh/fast/edlefsen_p/bakeoff_analysis_results"
#' @param force.recomputation if FALSE (default) and if there is a saved version called timings.results.by.region.and.time.Rda (under bakeoff_analysis_results/results.dirname), then that file will be loaded; otherwise the results will be recomputed and saved in that location.
#'
#' @return NULL
#' @export

evaluateTimings <- function (
                             use.center.of.bounds = TRUE,
                             use.infer = TRUE,
                             use.anchre = TRUE,
                             results.dirname = "raw_edited_20160216",
                             force.recomputation = FALSE
                            )
{
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
    
    compute.results.one.per.ppt <- function ( results, weights ) {
        apply( results, 2, function ( .column ) {
            .rv <- 
            sapply( unique( rownames( results ) ), function( .ppt ) {
                .ppt.cells <- .column[ rownames( results ) == .ppt ];
                if( all( is.na( .ppt.cells ) ) ) {
                  return( NA );
                }
                .ppt.weights <- weights[ rownames( results ) == .ppt ];
                .ppt.weights <- .ppt.weights / sum( .ppt.weights, na.rm = TRUE );
                sum( .ppt.cells * .ppt.weights, na.rm = TRUE );
            } );
            names( .rv ) <- unique( rownames( results ) );
            return( .rv );
        } );
    } # compute.results.one.per.ppt
    
    compute.diffs.by.stat <- function ( results.one.per.ppt, days.since.infection ) {
        diffs.by.stat <- 
            lapply( colnames( results.one.per.ppt ), function( .stat ) {
                .rv <- ( as.numeric( results.one.per.ppt[ , .stat ] ) - as.numeric( days.since.infection[ rownames( results.one.per.ppt ) ] ) );
                names( .rv ) <- rownames( results.one.per.ppt );
                return( .rv );
            } );
        names( diffs.by.stat ) <- colnames( results.one.per.ppt );
        return( diffs.by.stat );
    } # compute.diffs.by.stat ( results.one.per.ppt, days.since.infection )
    
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
                       
                   ## identify-founders results; we always get and use these.
                   if( is.na( partition.size ) ) {
                       results.in <- readIdentifyFounders( paste( paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/identify_founders.tab", sep = "" ) ) );
                   } else {
                       results.in <- readIdentifyFounders( paste( paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/partitions/identify_founders.tab", sep = "" ) ), partition.size = partition.size );
                   }
                   
                   results <- results.in;
                   
                   days.colnames <- c( grep( "time", colnames( results ), value = T ), grep( "days", colnames( results ), value = T ) );
    
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
                       function ( days.est.col.i, for.bias = FALSE, weights = days.est.nseq*days.est.nb, held.out.ptids = c() ) {
                           .retained.rows <- which( !( rownames( results ) %in% held.out.ptids ) );
                           return( function( log.epsilon, return.days.est.new = FALSE, retained.rows = .retained.rows ) {
                           .days.est.new <- 
                             apply( results[ retained.rows, , drop = FALSE ], 1, function( .row ) {
                                 if( is.na( .row[ days.est.colnames.nb[ days.est.col.i ] ] ) ) {
                                     return( NA );
                                 }
                                 return( daysFromLambda( .row[ lambda.est.colnames[ days.est.col.i ] ], .row[ days.est.colnames.nb[ days.est.col.i ] ], inverse.epsilon = ( 1 / exp( log.epsilon ) ) ) );
                             } );
                           .mat <- matrix( .days.est.new, ncol = 1 );
                           colnames( .mat ) <- days.est.colnames[ days.est.col.i ];
                           rownames( .mat ) <- rownames( days.est )[ retained.rows ];
                           if( return.days.est.new ) {
                               return( .mat );
                           }
                           .results.one.per.ppt <-
                               compute.results.one.per.ppt( .mat, weights[ retained.rows, ] );
                           .diffs.by.stat <- compute.diffs.by.stat( .results.one.per.ppt, days.since.infection );
                           if( for.bias ) {
                             .bias <- mean( .diffs.by.stat[[1]], na.rm = TRUE );
                             return( abs( .bias ) );
                           } else {
                               .rmse <- rmse( .diffs.by.stat[[1]], na.rm = TRUE );
                               return( .rmse );
                           }
                         } );
                       } # create.function.to.be.optimized (..)

                   ## Optimized for minimal RMSE, holding each ppt out in turn
                   optimal.for.rmse.log.epsilons <-
                       sapply( 1:length( days.est.colnames ), function( .col.i ) {
                           ## Hold each ptid out, one at a time.
                           sapply( rownames( results ), function( .ptid ) {
                               function.to.be.optimized <- create.function.to.be.optimized( .col.i, held.out.ptids = .ptid );
                               optimize( function.to.be.optimized, interval = c( log( default.epsilon ) - log( 1000 ), log( default.epsilon ) + log( 1000 ) ) )[ 1 ];
                           } );
                       } );
                   colnames( optimal.for.rmse.log.epsilons ) <- days.est.colnames;
                   rownames( optimal.for.rmse.log.epsilons ) <- rownames( results );
                   optimal.for.rmse.validation.results <- matrix( NA, nrow = nrow( results ), ncol = length( days.est.colnames ) );
                   for( .row.i in 1:nrow( results ) ) {
                       for( .col.i in 1:length( days.est.colnames ) ) {
                           .optimal.for.rmse.log.epsilon <- optimal.for.rmse.log.epsilons[[ .row.i, .col.i ]];
                           optimal.for.rmse.validation.results[ .row.i, .col.i ] <- 
                               create.function.to.be.optimized( .col.i )( .optimal.for.rmse.log.epsilon, return.days.est.new = TRUE, retained.rows = .row.i );
                       }
                   }
                   colnames( optimal.for.rmse.validation.results ) <-
                       paste( "optimal.for.rmse.validation", days.est.colnames, sep = "." );
                   rownames( optimal.for.rmse.validation.results ) <-
                       rownames( results );
    
                   ## Also do it for for.bias results (objective fn to minimize is abs( bias ))
                   ## Optimized for minimal bias, holding each ppt out in turn
                   optimal.for.bias.log.epsilons <-
                       sapply( 1:length( days.est.colnames ), function( .col.i ) {
                           ## Hold each ptid out, one at a time.
                           sapply( rownames( results ), function( .ptid ) {
                               function.to.be.optimized <- create.function.to.be.optimized( .col.i, for.bias = TRUE, held.out.ptids = .ptid );
                               optimize( function.to.be.optimized, interval = c( log( default.epsilon ) - log( 1000 ), log( default.epsilon ) + log( 1000 ) ) )[ 1 ];
                           } );
                       } );
                   colnames( optimal.for.bias.log.epsilons ) <- days.est.colnames;
                   rownames( optimal.for.bias.log.epsilons ) <- rownames( results );
                   optimal.for.bias.validation.results <- matrix( NA, nrow = nrow( results ), ncol = length( days.est.colnames ) );
                   for( .row.i in 1:nrow( results ) ) {
                       for( .col.i in 1:length( days.est.colnames ) ) {
                           .optimal.for.bias.log.epsilon <- optimal.for.bias.log.epsilons[[ .row.i, .col.i ]];
                           optimal.for.bias.validation.results[ .row.i, .col.i ] <- 
                               create.function.to.be.optimized( .col.i, for.bias = TRUE )( .optimal.for.bias.log.epsilon, return.days.est.new = TRUE, retained.rows = .row.i );
                       }
                   }
                   colnames( optimal.for.bias.validation.results ) <-
                       paste( "optimal.for.bias.validation", days.est.colnames, sep = "." );
                   rownames( optimal.for.bias.validation.results ) <-
                       rownames( results );
                   
                   results <-
                       cbind(
                           results[ , days.est.colnames, drop = FALSE ],
                           optimal.for.rmse.validation.results,
                           optimal.for.bias.validation.results
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
                         colnames( infer.results ) <- c( "Infer", "Infer.CI.high", "Infer.CI.low" );
                         ## reorder them
                         infer.results <- infer.results[ , c( "Infer", "Infer.CI.low", "Infer.CI.high" ), drop = FALSE ];
                         infer.results.bounds.ptid <- gsub( "^.+_(\\d+)/.+$", "\\1", names( infer.results.list ) );
                         rownames( infer.results ) <- infer.results.bounds.ptid;
                         
                         ## Separate it into separate tables by bounds type.
                         infer.results.bounds.type <-
                           gsub( "^.+_artificialBounds_(.+)_\\d+/.+$", "\\1", names( infer.results.list ) );
                         infer.results.bounds.type[ grep( "csv$", infer.results.bounds.type ) ] <- NA;
                         infer.results.nobounds.table <- infer.results[ is.na( infer.results.bounds.type ), , drop = FALSE ];
                         infer.results.bounds.types <- setdiff( unique( infer.results.bounds.type ), NA );
                         infer.results.bounds.tables <- lapply( infer.results.bounds.types, function ( .bounds.type ) {
                           return( infer.results[ !is.na( infer.results.bounds.type ) & ( infer.results.bounds.type == .bounds.type ), , drop = FALSE ] );
                         } );
                         names( infer.results.bounds.tables ) <- infer.results.bounds.types;
                         
                         # Add just the estimates from infer (not the CIs) to the results table.
                         new.results.columns <-
                           matrix( NA, nrow = nrow( results ), ncol = 1 + length( infer.results.bounds.tables ) );
                         rownames( new.results.columns ) <- rownames( results );
                         colnames( new.results.columns ) <- c( "Infer.time.est", infer.results.bounds.types );
                         
                         .shared.ptids.nobounds <-
                             intersect( rownames( new.results.columns ), rownames( infer.results.nobounds.table ) );
                         .result.ignored <- sapply( .shared.ptids.nobounds, function( .ptid ) {
                             .infer.subtable <-
                                 infer.results.nobounds.table[ rownames( infer.results.nobounds.table ) == .ptid, 1, drop = FALSE ];
                             stopifnot( sum( rownames( infer.results.nobounds.table ) == .ptid ) == nrow( .infer.subtable ) );
                             # But there might be fewer of these than there are rows in the results.table (if eg there are 3 input fasta files and infer results for only 2 of them).
                             stopifnot( nrow( .infer.subtable ) <= sum( rownames( new.results.columns ) == .ptid ) );
                             if( nrow( .infer.subtable ) < sum( rownames( new.results.columns ) == .ptid ) ) {
                               .infer.subtable <- rbind( .infer.subtable, matrix( NA, nrow = ( sum( rownames( new.results.columns ) == .ptid ) - nrow( .infer.subtable ) ), ncol = ncol( .infer.subtable ) ) );
                             }
                             new.results.columns[ rownames( new.results.columns ) == .ptid, 1 ] <<-
                                 .infer.subtable;
                             return( NULL );
                         } );
                         .result.ignored <- sapply( names( infer.results.bounds.tables ), function ( .bounds.type ) {
                           .shared.ptids <-
                               intersect( rownames( new.results.columns ), rownames( infer.results.bounds.tables[[ .bounds.type ]] ) );
                           ..result.ignored <- sapply( .shared.ptids, function( .ptid ) {
                               .infer.subtable <-
                                 infer.results.bounds.tables[[ .bounds.type ]][ rownames( infer.results.bounds.tables[[ .bounds.type ]] ) == .ptid, 1, drop = FALSE ];
                               stopifnot( sum( rownames( infer.results.bounds.tables[[ .bounds.type ]] ) == .ptid ) == nrow( .infer.subtable ) );
                               # But there might be fewer of these than there are rows in the results.table (if eg there are 3 input fasta files and infer results for only 2 of them).
                               stopifnot( nrow( .infer.subtable ) <= sum( rownames( new.results.columns ) == .ptid ) );
                               if( nrow( .infer.subtable ) < sum( rownames( new.results.columns ) == .ptid ) ) {
                                 .infer.subtable <- rbind( .infer.subtable, matrix( NA, nrow = ( sum( rownames( new.results.columns ) == .ptid ) - nrow( .infer.subtable ) ), ncol = ncol( .infer.subtable ) ) );
                               }
                               new.results.columns[ rownames( new.results.columns ) == .ptid, .bounds.type ] <<-
                                   .infer.subtable;
                               return( NULL );
                           } );                       
                           return( NULL );
                         } );
                           
                         colnames( new.results.columns ) <- c( "Infer.time.est", paste( "Infer", gsub( "_", ".", infer.results.bounds.types ), "time.est", sep = "." ) ); 
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

                   results.covars.one.per.ppt.with.extra.cols <-
                       summarizeCovariatesOnePerParticipant( results.in );
                     
                   .keep.cols <-
                       grep( "num.*\\.seqs|totalbases", colnames( results.covars.one.per.ppt.with.extra.cols ), value = TRUE, perl = TRUE, invert = TRUE );
                   ### TODO: Something else.  Just trying to get down to a reasonable set; basically there are very highly clustered covariates here and it screws up the inference.
                   ## Also remove all of the mut.rate.coef except for multifounder.Synonymous.PFitter.mut.rate.coef.
                   
                   # .keep.cols <- c( "multifounder.Synonymous.PFitter.mut.rate.coef",
                   #                 grep( "mut\\.rate\\.coef", .keep.cols, invert = TRUE, value = TRUE ) )
                   #.keep.cols <- c( "multifounder.Synonymous.PFitter.mut.rate.coef", "inf.to.priv.ratio", "priv.sites", "inf.sites.clusters", "InSites.founders", "multifounder.Synonymous.PFitter.is.poisson" );
                   ## Keep only the mut.rate.coef cols and priv.sites and multifounder.Synonymous.PFitter.is.poisson.
                   helpful.additional.cols <- c( "priv.sites","multifounder.Synonymous.PFitter.is.poisson" );
                   mutation.rate.coefs <- grep( "mut\\.rate\\.coef", .keep.cols, value = TRUE );
                   .keep.cols <- c( mutation.rate.coefs, helpful.additional.cols );
                   results.covars.one.per.ppt <-
                       results.covars.one.per.ppt.with.extra.cols[ , .keep.cols, drop = FALSE ];
                   results.covars.one.per.ppt.df <-
                       data.frame( results.covars.one.per.ppt );

                   regression.df <- cbind( data.frame( days.since.infection = days.since.infection[ rownames( results.covars.one.per.ppt.df ) ] ), results.covars.one.per.ppt.df );
                   
                   ## ERE I AM...
                   # library( "glmnet" )
                   # cv.glmnet.fit <- cv.glmnet( results.covars.one.per.ppt, days.since.infection, nfolds = nrow( results.covars.one.per.ppt ), type.measure = "mae", grouped = FALSE, intercept = FALSE ); # mean absolute error, corresponding to the "for.bias" version of the other exploration.
                   
                   #library( "boot" );
                   ## ENDMARK
                   # gaussian.fit.formula <- as.formula( paste( "days.since.infection ~ 0 + ", paste( colnames( results.covars.one.per.ppt ), collapse = "+" ) ) );
                   # gaussian.fit <- glm( gaussian.fit.formula, family = "gaussian", data = regression.df );
                   # summary( gaussian.fit );
                   #cv.glm( data = regression.df, glmfit = gaussian.fit, K = nrow( regression.df ) );
                   ## new proof of concept:
                   helpful.additional.parameters.validation.results.one.per.ppt <- matrix( NA, nrow = nrow( results.covars.one.per.ppt.df ), ncol = length( days.est.colnames ) );
                   for( .row.i in 1:nrow( regression.df ) ) {
                       for( .col.i in 1:length( days.est.colnames ) ) {
                           .mut.rate.coef.colname <- colnames( mutation.rate.coefs )[ .col.i ];
                           ## Ok build a regression model with no intercept, including only the helpful.additional.cols
                           .formula <- as.formula( paste( "days.since.infection ~ 0 + ", paste( c( helpful.additional.cols, .mut.rate.coef.colname ), collapse = "+" ) ) );
                           .pred.value <- predict( lm( .formula, data = regression.df[ -.row.i, ] ), regression.df[ .row.i, , drop = FALSE ] );
                           helpful.additional.parameters.validation.results.one.per.ppt[ .row.i, .col.i ] <- 
                               .pred.value;
                       } # End foreach .col.i
                   } # End foreach .row.i
                   colnames( helpful.additional.parameters.validation.results.one.per.ppt ) <-
                       paste( "helpful.additional.cols.validation", days.est.colnames, sep = "." );
                   rownames( helpful.additional.parameters.validation.results.one.per.ppt ) <-
                       rownames( regression.df );
                   
                     results.one.per.ppt <-
                         cbind( results.one.per.ppt,
                               helpful.additional.parameters.validation.results.one.per.ppt );
                     
                     if( use.center.of.bounds ) {
                       ## The "center of bounds" approach is the way that
                       ## we do it at the VTN, and the way it was done in
                       ## RV144, etc: use the midpoint between the bounds
                       ## on the actual infection time computed from the
                       ## dates and results of the HIV positivity tests
                       ## (antibody or PCR).  The typical approach is to
                       ## perform antibody testing every X days
                       ## (historically this is 6 months in most HIV
                       ## vaccine trials, except during the vaccination
                       ## phase there are more frequent visits and on
                       ## every visit HIV testing is conducted).  The
                       ## (fake) bounds used here are calculated in the
                       ## createArtificialBoundsOnInfectionDate.R file.
                       ## The actual bounds would be too tight, since the
                       ## participants were detected HIV+ earlier in these
                       ## people than what we expect to see in a trial in
                       ## which testing is conducted every X days.  For
                       ## the center of bounds approach we load the bounds
                       ## files in subdirs of the "bounds" subdirectory eg
                       ## at
                       ## /fh/fast/edlefsen_p/bakeoff/analysis_sequences/bounds/nflg/1m/.
                       ## These files have names beginning with
                       ## "artificialBounds_" and ending with ".tab".
                       bounds.subdirname <- "bounds";
                       .artificial.bounds.dirname <-
                           paste( "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/", results.dirname, "/", bounds.subdirname, "/", the.region, "/", the.time, "/", sep = "" );
                       artificial.bounds.filenames <-
                           dir( .artificial.bounds.dirname, pattern = "artificialBounds_.*.tab", recursive = FALSE, full.names = TRUE );
                       names( artificial.bounds.filenames ) <- gsub( "^.*artificialBounds_(.*).tab$", "\\1", artificial.bounds.filenames );
                       the.artificial.bounds <- lapply( names( artificial.bounds.filenames ), function ( .artificial.bounds.name ) {
                         return( read.table( artificial.bounds.filenames[[ .artificial.bounds.name ]], header = TRUE, sep = "\t" ) );
                       } );
                       names( the.artificial.bounds ) <- names( artificial.bounds.filenames );
                       
                       center.of.bounds.table <- sapply( names( the.artificial.bounds ), function ( .artificial.bounds.name ) {
                           round( apply( the.artificial.bounds[[ .artificial.bounds.name ]], 1, mean ) )
                       } );
                       colnames( center.of.bounds.table ) <-
                         paste( "COB", gsub( "_", ".", colnames( center.of.bounds.table ) ), "time.est", sep = "." );
    
                       ## Ok, well, now we can also evaluate versions of variants of each
                       ## method, but bounding the results.  Ie if the
                       ## estimated time is within the bounds, that time
                       ## is used, otherwise, the boundary.
                       
                       results.one.per.ppt.bounded <-
                         lapply( names( the.artificial.bounds ), function ( .artificial.bounds.name ) {
                           .mat <-
                           apply( results.one.per.ppt, 2, function ( .results.column ) {
                           sapply( names( .results.column ), function ( .ppt ) {
                             .value <- .results.column[ .ppt ];
                             if( !is.na( .value ) ) {
                               .value.is.below.lb <- ( .value < the.artificial.bounds[[ .artificial.bounds.name ]][ .ppt, "lower" ] );
                               if( .value.is.below.lb ) {
                                 return( the.artificial.bounds[[ .artificial.bounds.name ]][ .ppt, "lower" ] );
                               }
                               .value.is.above.ub <- ( .value > the.artificial.bounds[[ .artificial.bounds.name ]][ .ppt, "upper" ] );
                               if( .value.is.above.ub ) {
                                 return( the.artificial.bounds[[ .artificial.bounds.name ]][ .ppt, "upper" ] );
                               }
                             }
                             return( .value );
                           } );
                         } );
                           rownames( .mat ) <- rownames( results.one.per.ppt );
                           return( .mat );
                       } );
                       names( results.one.per.ppt.bounded ) <- names( the.artificial.bounds );
    
                       ## NOTE: we add the COB bounds only to the
                       ## original, unbounded results, after computing the
                       ## bounded ones.
                       results.one.per.ppt <-
                         cbind( results.one.per.ppt, center.of.bounds.table[ rownames( results.one.per.ppt ), , drop = FALSE ] );
                     } # End if use.center.of.bounds
    
                     ## unbounded results:
                     diffs.by.stat <- compute.diffs.by.stat( results.one.per.ppt );
    
                     unbounded.results <- 
                       list( bias = lapply( diffs.by.stat, mean, na.rm = T ), se = lapply( diffs.by.stat, sd, na.rm = T ), rmse = lapply( diffs.by.stat, rmse, na.rm = T ), n = lapply( diffs.by.stat, function( .vec ) { sum( !is.na( .vec ) ) } ) );
                     
                     ## bounded results:
                     if( use.center.of.bounds ) {
                       bounded.results.by.bound.type <- lapply( results.one.per.ppt.bounded, function ( .results.one.per.ppt ) {
                         .diffs.by.stat <- compute.diffs.by.stat( .results.one.per.ppt );
        
                         return( list( bias = lapply( .diffs.by.stat, mean, na.rm = T ), se = lapply( .diffs.by.stat, sd, na.rm = T ), rmse = lapply( .diffs.by.stat, rmse, na.rm = T ), n = lapply( .diffs.by.stat, function( .vec ) { sum( !is.na( .vec ) ) } ) ) );
                       } );
                       return( c( list( unbounded = unbounded.results ), bounded.results.by.bound.type ) );
                     } else {
                       return( unbounded.results );
                     }
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
               } );
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
    
    return( invisible( NULL ) );
} # evaluateTimings (..)

## Here is where the action is.
evaluateTimings();

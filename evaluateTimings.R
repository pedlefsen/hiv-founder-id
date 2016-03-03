## First do all the stuff in README.postprocessing.txt.

use.infer <- TRUE;
use.anchre <- FALSE;
results.dirname <- "raw_edited_20160216";
#results.dirname <- "raw";

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

timings.results.by.region.and.time <- 
 lapply( regions, function( the.region ) {
             ## TODO: REMOVE
             cat( the.region, fill = T );
     timings.results.by.time <- 
         lapply( times, function( the.time ) {
             ## TODO: REMOVE
             cat( the.time, fill = T );
             sample.dates.in <- read.delim( paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/sampleDates.tbl", sep = "" ), sep = " ", header = F, fill = T );
             colnames( sample.dates.in ) <- c( "ptid", "date" );
             ## Remove anything that's not really a ptid/date combo
             sample.dates.in <- sample.dates.in[ grep( "^\\d+$", as.character( sample.dates.in[ , 1 ] ) ), , drop = FALSE ];
             # remove anything with a missing date.
             sample.dates.in <-
                 sample.dates.in[ sample.dates.in[ , 2 ] != "", , drop = FALSE ];
             # Special: for v3, separate out the caprisa seqs from the rv217 seqs
#              if( the.region == "v3" ) {
#                  sample.dates.in.rv217 <- sample.dates.in[ grep( "^100\\d\\d\\d", as.character( sample.dates.in[ , 1 ] ), invert = TRUE ), , drop = FALSE ];
#                  sample.dates.in <- sample.dates.in[ grep( "^100\\d\\d\\d", as.character( sample.dates.in[ , 1 ] ) ), , drop = FALSE ];
#              }
             days.since.infection <- sapply( 1:nrow( sample.dates.in ), function( .i ) { as.numeric( as.Date( as.character( sample.dates.in[ .i, 2 ] ) ) - ifelse( the.region == "v3", caprisa002.gold.standard.infection.dates[ as.character( sample.dates.in[ .i, 1 ] ) ], rv217.gold.standard.infection.dates[ as.character( sample.dates.in[ .i, 1 ] ) ] ) ) } );
             names( days.since.infection ) <- sample.dates.in[ , "ptid" ];
#              if( the.region == "v3" ) {
#                  days.since.infection.rv217 <- sapply( 1:nrow( sample.dates.in.rv217 ), function( .i ) { as.numeric( as.Date( as.character( sample.dates.in.rv217[ .i, 2 ] ) ) - rv217.gold.standard.infection.dates[ as.character( sample.dates.in.rv217[ .i, 1 ] ) ] ) } );
#                  names( days.since.infection.rv217 ) <- sample.dates.in.rv217[ , "ptid" ];
#              }
                 
             ## identify-founders results
             #results.in <- read.delim( paste( paste( "/fh/fast/edlefsen_p/bakeoff_merged_analysis_sequences/raw_fixed/", the.region, "/", the.time, "/identify_founders.tab", sep = "" ) ), sep = "\t" );
             results.in <- read.delim( paste( paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/identify_founders.tab", sep = "" ) ), sep = "\t" );
             # Hack to fix a bug in which the colnames don't have "multifounder." on them (but instead, because of how read.delim works, they have ".1" at the end).
             #colnames( results.in ) <- gsub( "^(.+)\\.1$", "multifounder.\\1", colnames( results.in ) );

             results <- as.matrix( results.in );
             suppressWarnings( mode( results ) <- "numeric" );
             
             rownames( results ) <- gsub( ".*caprisa002_(\\d+)_.*", "\\1", gsub( ".*rv217_(\\d+)_.*", "\\1", as.character( results.in[ , 1 ] ) ) );


             ## ERE I AM, trying something about recalibrating the rate.
             epsilon <- 2.16e-05;
             phi <- sqrt(1+4/3)
             days <- function(l,nb,epsilon) 1.5*((phi)/(1+phi))*(l/(epsilon*nb) - (1-phi)/(phi^2))
             
             lambda.colnames <- grep( "lambda", colnames( results ), value = T );
             lambda <- results[ , lambda.colnames, drop = FALSE ];
             lambda.colnames.nb <- gsub( "[^\\.]+\\.Star[Pp]hy", "PFitter", gsub( "lambda.*$", "nbases", lambda.colnames ) );
             lambda.nb <- results[ , lambda.colnames.nb, drop = FALSE ];
             
             days.colnames <- grep( "days", colnames( results ), value = T );
             days <- results[ , days.colnames, drop = FALSE ];
             days.colnames.nb <- gsub( "[^\\.]+\\.Star[Pp]hy", "PFitter", gsub( "days.*$", "nbases", days.colnames ) );
             days.nb <- results[ , days.colnames.nb, drop = FALSE ];

             single.colnames <- grep( "\\.is\\.", colnames( results ), value = T );
             single.exceptInSites.colnames <- grep( "InSites", single.colnames, value = T, invert = T );
             single.exceptInSites <- results[ , single.exceptInSites.colnames, drop = FALSE ];
             single.exceptInSites.colnames.nb <-
                 gsub( "^Star[Pp]hy", "PFitter", gsub( "[^\\.]+\\.Star[Pp]hy", "PFitter", gsub( "\\....er\\.PFitter\\.", "\\.", gsub( "is\\..*$", "nbases", single.exceptInSites.colnames ) ) ) );
             ## Special case: the one called "StarPhy.is.one.founder" is actually using the synonymous PFitter results, so its "nbases" should be the synonymous one.
             single.exceptInSites.colnames.nb[ single.exceptInSites.colnames == "StarPhy.is.one.founder" ] <- "Synonymous.PFitter.nbases";
             single.exceptInSites.nb <- results[ , single.exceptInSites.colnames.nb, drop = FALSE ];
             
             ## Results for which nbases == 0 should be NAs!  Remove the time and lambda estimates.
             results[ , lambda.colnames ] <-
                 t( apply( results, 1, function( .row ) {
                     ifelse( .row[ lambda.colnames.nb ] == 0, NA, .row[ lambda.colnames ] )
                 } ) );
             results[ , days.colnames ] <-
                 t( apply( results, 1, function( .row ) {
                     ifelse( .row[ days.colnames.nb ] == 0, NA, .row[ days.colnames ] )
                 } ) );
             ## TODO: [in evaluateIsMultiple] Fix the single ones there
             results[ , single.exceptInSites.colnames ] <-
                 t( apply( results, 1, function( .row ) {
                     ifelse( .row[ single.exceptInSites.colnames.nb ] == 0, NA, .row[ single.exceptInSites.colnames ] )
                 } ) );

             ### ERE I AM, working on this.
             
             lambda.est.colnames <- grep( "\\.ci\\.", lambda.colnames, value = T, invert = TRUE );
             lambda.est <- results[ , lambda.est.colnames, drop = FALSE ];
             nb.est.colnames <- gsub( "[^\\.]+\\.Star[Pp]hy", "PFitter", gsub( "lambda(\\.est)?", "nbases", lambda.est.colnames ) );
             nb.est <- results[ , nb.est.colnames, drop = FALSE ];

             results <- results[ , identify.founders.date.estimates, drop = FALSE ];
#              if( the.region == "v3" ) {
#                  results.rv217 <- results.rv217[ , identify.founders.date.estimates, drop = FALSE ];
#              }

             if( use.infer ) {
               ## Add to results: "infer" results.
                                        #infer.results.directories <- dir( paste( "/fh/fast/edlefsen_p/bakeoff_merged_analysis_sequences/raw_fixed/", the.region, "/", the.time, sep = "" ), "founder-inference-bakeoff_", full.name = TRUE );
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

               if( length( infer.results.list ) > 0 ) {
                   infer.results <- do.call( rbind, infer.results.list );
                   colnames( infer.results ) <- c( "Infer", "Infer.CI.low", "Infer.CI.high" );
                   rownames( infer.results ) <-
                       gsub( "^.+_(\\d+)$", "\\1", names( unlist( infer.results.files ) ) );
                   
                   # Add just the estimate from infer.
                   new.results.columns <- matrix( nrow = nrow( results ), ncol = 1 );
                   rownames( new.results.columns ) <- rownames( results );
                   colnames( new.results.columns ) <- "Infer.time.est";
                   new.results.columns[ rownames( infer.results ), 1 ] <- infer.results[ , 1 ];
                     
                   results <- cbind( results, new.results.columns );
               }
                 
             } # End if use.infer

             if( use.anchre && ( the.time == "1m6m" ) ) {
               ## Add to results: "anchre" results. (only at 1m6m)
               anchre.results.directories <- dir( paste( "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw/", the.region, "/1m6m", sep = "" ), "anchre", full.name = TRUE );
               anchre.results.files <-
                   sapply( anchre.results.directories, dir, "mrca.csv", full.name = TRUE );
               anchre.results <- do.call( rbind,
                   lapply( unlist( anchre.results.files ), function( .file ) {
                       stopifnot( file.exists( .file ) );
                       .file.short <-
                           gsub( "^.*?\\/?([^\\/]+?)$", "\\1", .file, perl = TRUE );
                       .file.short.nosuffix <-
                           gsub( "^([^\\.]+)(\\..+)?$", "\\1", .file.short, perl = TRUE );
                       .file.converted <-
                           paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw_fixed/", the.region, "/1m6m/", .file.short.nosuffix, ".anc2tsv.tab", sep = "" );
                       # convert it.
                       system( paste( "./anc2tsv.sh", .file, ">", .file.converted ) );
                       stopifnot( file.exists( .file.converted ) );
                       .rv <- as.matrix( read.delim( .file.converted, header = TRUE, sep = "\t" ), nrow = 1 );
                       ## No negative dates!
                       .rv <- apply( .rv, 1:2, function( .str ) { if( length( grep( "^-", .str ) ) > 0 ) { "0001-01-01" } else { .str } } );
                       stopifnot( ncol( .rv ) == 4 );
                       return( .rv );
                   } ) );
               colnames( anchre.results ) <- c( "Anchre.r2t.est", "Anchre.est", "Anchre.CI.low", "Anchre.CI.high" );
               rownames( anchre.results ) <-
                   gsub( "^.+_(\\d+)$", "\\1", names( unlist( anchre.results.files ) ) );
               # Special: for v3, only use caprisa seqs (not rv217, for now).
               if( the.region == "v3" ) {
                   anchre.results <-
                       anchre.results[ grep( "^100\\d\\d\\d", rownames( anchre.results ) ), , drop = FALSE ];
               }
               # Add just the estimate from anchre.
               sample.dates <- as.Date( as.character( sample.dates.in[ , 2 ] ) );
               names( sample.dates ) <- sample.dates.in[ , 1 ];
               anchre.r2t.days.before.sample <- sapply( 1:nrow( anchre.results ), function( .i ) { 0 - as.numeric( as.Date( anchre.results[ .i, 1 ] ) - sample.dates[ rownames( anchre.results )[ .i ] ] ) } );
               names( anchre.r2t.days.before.sample ) <- rownames( anchre.results );
               anchre.days.before.sample <- sapply( 1:nrow( anchre.results ), function( .i ) { 0 - as.numeric( as.Date( anchre.results[ .i, 2 ] ) - sample.dates[ rownames( anchre.results )[ .i ] ] ) } );
               names( anchre.days.before.sample ) <- rownames( anchre.results );
  
               results <- cbind( results, anchre.r2t.days.before.sample, anchre.days.before.sample );
               colnames( results )[ ncol( results ) - 1:0 ] <- c( "Anchre.r2t", "Anchre.bst" ) ;
             } # End if use.ancher and the.time is 1m6m, add anchre results too.
             
             diffs.by.stat <-
                 lapply( colnames( results ), function( .stat ) {
                     as.numeric( results[ , .stat ] ) - as.numeric( days.since.infection[ rownames( results ) ] );
                 } );
             names( diffs.by.stat ) <- colnames( results );

             # Special: for v3, separate out the caprisa seqs from the rv217 seqs
             if( the.region == "v3" ) {
                 diffs.by.stat.rv217 <-
                     lapply( colnames( results.rv217 ), function( .stat ) {
                         as.numeric( results.rv217[ , .stat ] ) - as.numeric( days.since.infection.rv217[ rownames( results.rv217 ) ] );
                     } );
                 names( diffs.by.stat.rv217 ) <- colnames( results.rv217 );
                 return( list( caprisa002 = list( bias = lapply( diffs.by.stat, mean, na.rm = T ), se = lapply( diffs.by.stat, sd, na.rm = T ), rmse = lapply( diffs.by.stat, rmse, na.rm = T ) ), rv217 = list( bias = lapply( diffs.by.stat.rv217, mean, na.rm = T ), se = lapply( diffs.by.stat.rv217, sd, na.rm = T ), rmse = lapply( diffs.by.stat.rv217, rmse, na.rm = T ) ) ) );
             } else {
                 return( list( rv217 = list( bias = lapply( diffs.by.stat, mean, na.rm = T ), se = lapply( diffs.by.stat, sd, na.rm = T ), rmse = lapply( diffs.by.stat, rmse, na.rm = T ) ) ) );
             }

         } ); # End foreach the.time
     names( timings.results.by.time ) <- times;
     return( timings.results.by.time );
 } ); # End foreach the.region
names( timings.results.by.region.and.time ) <- regions;

# Make a table out of it. (one per study).
results.tables <-
    lapply( names( timings.results.by.region.and.time ), function( the.region ) {
        .rv <- 
            lapply( names( timings.results.by.region.and.time[[ the.region ]] ), function( the.time ) {
                ..rv <- 
                    lapply( names( timings.results.by.region.and.time[[ the.region ]][[ the.time ]] ), function( the.study ) {
                        sapply( timings.results.by.region.and.time[[ the.region ]][[ the.time ]][[ the.study ]], function( results.list ) { results.list } ) }
                        );
                names( ..rv ) <- names( timings.results.by.region.and.time[[ the.region ]][[ the.time ]] );
                return( ..rv );
            } );
        names( .rv ) <- times;
        return( .rv );
    } );
names( results.tables ) <- regions;

## ERE I AM.  Should write these out somehow .
# results.tables


results.tables


### Raw run, without recombination detection/removal nor hypermutation detection/removal:
# $nflg
# $nflg$`1m`
# $nflg$`1m`$rv217
#                                          bias      se       rmse    
# PFitter.time.est                         43.9375   126.8828 133.3349
# Synonymous.PFitter.time.est              -23.5     28.77058 36.97381
# multifounder.PFitter.time.est            4.666667  55.83858 55.48503
# multifounder.Synonymous.PFitter.time.est -31.29412 14.19337 34.30486
# 
# 
# $nflg$`6m`
# $nflg$`6m`$rv217
#                                          bias      se       rmse    
# PFitter.time.est                         -28.87671 117.9385 120.6351
# Synonymous.PFitter.time.est              -145.3699 29.46491 148.2858
# multifounder.PFitter.time.est            -71.72222 83.44774 109.447 
# multifounder.Synonymous.PFitter.time.est -152.7778 27.21034 155.1378
# 
# 
# $nflg$`1m6m`
# $nflg$`1m6m`$rv217
#                                          bias      se       rmse    
# PFitter.time.est                         87.36842  122.4774 149.1282
# Synonymous.PFitter.time.est              -11.84211 29.62698 31.54195
# multifounder.PFitter.time.est            26        65.86704 69.81774
# multifounder.Synonymous.PFitter.time.est -24.3871  19.71409 31.15828
# 
# 
# 
# $v3
# $v3$`1m`
# $v3$`1m`$caprisa002
#                                          bias   se       rmse    
# PFitter.time.est                         164.65 340.783  370.7237
# Synonymous.PFitter.time.est              -23    58.29869 61.3009 
# multifounder.PFitter.time.est            12     111.4757 109.3138
# multifounder.Synonymous.PFitter.time.est -41.3  31.20746 51.2923 
# 
# $v3$`1m`$rv217
#                                          bias      sd       rmse    
# PFitter.time.est                         41.55556  144.4919 148.4077
# Synonymous.PFitter.time.est              -28.27778 37.03096 46.18261
# multifounder.PFitter.time.est            -2.305556 68.9268  68.00184
# multifounder.Synonymous.PFitter.time.est -32.5     24.94623 40.75878
# 
# 
# $v3$`6m`
# $v3$`6m`$caprisa002
#                                          bias      se       rmse    
# PFitter.time.est                         161.6667  253.3017 294.5052
# Synonymous.PFitter.time.est              -119.0556 76.703   140.4661
# multifounder.PFitter.time.est            30.22222  262.8089 257.1863
# multifounder.Synonymous.PFitter.time.est -141.8333 60.15592 153.4092
# 
# $v3$`6m`$rv217
#                                          bias      sd       rmse    
# PFitter.time.est                         40.09091  184.6078 186.1574
# Synonymous.PFitter.time.est              -125      59.86443 138.2033
# multifounder.PFitter.time.est            -64       77.77612 99.80891
# multifounder.Synonymous.PFitter.time.est -142.8182 30.8063  146.0045
# 
# 
# $v3$`1m6m`
# $v3$`1m6m`$caprisa002
#                                          bias      se       rmse    
# PFitter.time.est                         249.2941  355.0214 425.1748
# Synonymous.PFitter.time.est              -13.47059 66.84788 66.23621
# multifounder.PFitter.time.est            132.0588  352.2764 366.3854
# multifounder.Synonymous.PFitter.time.est -33.23529 46.61616 56.12329
# 
# $v3$`1m6m`$rv217
#                                          bias      sd       rmse    
# PFitter.time.est                         161       210.6419 262.5766
# Synonymous.PFitter.time.est              0.1515152 54.44844 53.61733
# multifounder.PFitter.time.est            30.33333  43.68328 52.63568
# multifounder.Synonymous.PFitter.time.est -16.75758 27.95424 32.22694


### Original run, using RAPBeta and my implementation of Hypermut 2.0.
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


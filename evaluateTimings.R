identify.founders.date.estimates <- c( "PFitter.time.est", "Synonymous.PFitter.time.est", "multifounder.PFitter.time.est", "multifounder.Synonymous.PFitter.time.est" ); 

rv217.gold.standard.infection.dates.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/rv217_gold_standard_timings.csv" );
rv217.gold.standard.infection.dates <- as.Date( as.character( rv217.gold.standard.infection.dates.in[,2] ), "%m/%d/%y" );
names( rv217.gold.standard.infection.dates ) <- as.character( rv217.gold.standard.infection.dates.in[,1] );

caprisa002.gold.standard.infection.dates.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/caprisa_002_gold_standard_timings.csv" );
caprisa002.gold.standard.infection.dates <- as.Date( as.character( caprisa002.gold.standard.infection.dates.in[,2] ), "%Y/%m/%d" );
names( caprisa002.gold.standard.infection.dates ) <- as.character( caprisa002.gold.standard.infection.dates.in[,1] );

studies <- c( "nflg", "v3" );
times <- c( "1m", "6m", "1m6m" );

timings.results.by.study.and.time <- 
 lapply( studies, function( the.study ) {
     timings.results.by.time <- 
         lapply( times, function( the.time ) {
             sample.dates.in <- read.delim( paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw/", the.study, "/", the.time, "/sampleDates.tbl", sep = "" ), sep = " ", header = F, fill = T );
             colnames( sample.dates.in ) <- c( "ptid", "date" );
             # remove anything with a missing date.
             sample.dates.in <-
                 sample.dates.in[ sample.dates.in[ , 2 ] != "", , drop = FALSE ];
             # Special: for v3, only use caprisa seqs (not rv144, for now).
             if( the.study == "v3" ) {
                 sample.dates.in <- sample.dates.in[ grep( "^100\\d\\d\\d", as.character( sample.dates.in[ , 1 ] ) ), , drop = FALSE ];
             }
             days.since.infection <- sapply( 1:nrow( sample.dates.in ), function( .i ) { as.numeric( as.Date( as.character( sample.dates.in[ .i, 2 ] ) ) - ifelse( the.study == "v3", caprisa002.gold.standard.infection.dates[ as.character( sample.dates.in[ .i, 1 ] ) ], rv217.gold.standard.infection.dates[ as.character( sample.dates.in[ .i, 1 ] ) ] ) ) } );
             names( days.since.infection ) <- sample.dates.in[ , "ptid" ];

             ## identify-founder results
             results.in <- read.delim( paste( paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw/", the.study, "/", the.time, "/identify_founders.tab", sep = "" ) ), sep = "\t" );
             # Hack to fix a bug in which the colnames don't have "multifounder." on them (but instead, because of how read.delim works, they have ".1" at the end).
             colnames( results.in ) <- gsub( "^(.+)\\.1$", "multifounder.\\1", colnames( results.in ) );
             results <- as.matrix( results.in );
             suppressWarnings( mode( results ) <- "numeric" );
             
             rownames( results ) <- gsub( ".*caprisa002_(\\d+)_.*", "\\1", gsub( ".*rv217_(\\d+)_.*", "\\1", as.character( results.in[ , 1 ] ) ) );
             # Special: for v3, only use caprisa seqs (not rv144, for now).
             if( the.study == "v3" ) {
                 results <- results[ grep( "^100\\d\\d\\d", rownames( results ) ), , drop = FALSE ];
             }

             results <- results[ , identify.founders.date.estimates, drop = FALSE ];

             ## Add to results: get the anchre, infer results.
             # ERE I AM.  founder-inference-bakeoff_10066/bakeoff_analysis_rv217_10066_1M_RH.outtoi.csv
             infer.results.directories <- dir( paste( "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw/", the.study, "/", the.time, sep = "" ), "founder-inference-bakeoff_", full.name = TRUE );
             infer.results.files <- sapply( infer.results.directories, dir, "outtoi.csv", full.name = TRUE );
             infer.results <- do.call( rbind,
                 lapply( unlist( infer.results.files ), function( .file ) {
                     return( as.matrix( read.csv( .file, header = FALSE ), nrow = 1 ) );
                 } ) );
             colnames( infer.results ) <- c( "MatsenWarth", "MatsenWarth.CI.low", "MatsenWarth.CI.high" );
             rownames( infer.results ) <-
                 gsub( "^.+_(\\d+)$", "\\1", names( unlist( infer.results.files ) ) );
             # Special: for v3, only use caprisa seqs (not rv144, for now).
             if( the.study == "v3" ) {
                 infer.results <-
                     infer.results[ grep( "^100\\d\\d\\d", rownames( infer.results ) ), , drop = FALSE ];
             }
             # Add just the estimate from infer.
             sample.dates <- as.Date( as.character( sample.dates.in[ , 2 ] ) );
             names( sample.dates ) <- sample.dates.in[ , 1 ];
             infer.days.before.sample <- sapply( 1:nrow( infer.results ), function( .i ) { 0 - as.numeric( as.Date( infer.results[ .i, 1 ] ) - sample.dates[ rownames( infer.results )[ .i ] ] ) } );
             names( infer.days.before.sample ) <- rownames( infer.results );

             results <- cbind( results, infer.days.before.sample );
             colnames( results )[ ncol( results ) ] <- "MatsenWarth";
             
             diffs.by.stat <-
                 lapply( colnames( results ), function( .stat ) {
                     as.numeric( results[ , .stat ] ) - as.numeric( days.since.infection[ rownames( results ) ] );
                 } );
             names( diffs.by.stat ) <- colnames( results );
             
             return( list( bias = lapply( diffs.by.stat, mean, na.rm = T ), rmse = lapply( diffs.by.stat, sd, na.rm = T ) ) );
         } ); # End foreach the.time
     names( timings.results.by.time ) <- times;
     return( timings.results.by.time );
 } ); # End foreach the.study
names( timings.results.by.study.and.time ) <- studies;




timings.results.by.study.and.time
# $nflg
# $nflg$`1m`
# $nflg$`1m`$bias
# $nflg$`1m`$bias$PFitter.time.est
# [1] 43.25
# 
# $nflg$`1m`$bias$Synonymous.PFitter.time.est
# [1] -23.625
# 
# $nflg$`1m`$bias$multifounder.PFitter.time.est
# [1] 6.137255
# 
# $nflg$`1m`$bias$multifounder.Synonymous.PFitter.time.est
# [1] -30.84314
# 
# 
# $nflg$`1m`$rmse
# $nflg$`1m`$rmse$PFitter.time.est
# [1] 120.3937
# 
# $nflg$`1m`$rmse$Synonymous.PFitter.time.est
# [1] 28.13319
# 
# $nflg$`1m`$rmse$multifounder.PFitter.time.est
# [1] 61.02033
# 
# $nflg$`1m`$rmse$multifounder.Synonymous.PFitter.time.est
# [1] 16.37116
# 
# 
# 
# $nflg$`6m`
# $nflg$`6m`$bias
# $nflg$`6m`$bias$PFitter.time.est
# [1] -26.36986
# 
# $nflg$`6m`$bias$Synonymous.PFitter.time.est
# [1] -145.2329
# 
# $nflg$`6m`$bias$multifounder.PFitter.time.est
# [1] -66.92593
# 
# $nflg$`6m`$bias$multifounder.Synonymous.PFitter.time.est
# [1] -151.3889
# 
# 
# $nflg$`6m`$rmse
# $nflg$`6m`$rmse$PFitter.time.est
# [1] 121.6416
# 
# $nflg$`6m`$rmse$Synonymous.PFitter.time.est
# [1] 29.92886
# 
# $nflg$`6m`$rmse$multifounder.PFitter.time.est
# [1] 91.80863
# 
# $nflg$`6m`$rmse$multifounder.Synonymous.PFitter.time.est
# [1] 28.42296
# 
# 
# 
# $nflg$`1m6m`
# $nflg$`1m6m`$bias
# $nflg$`1m6m`$bias$PFitter.time.est
# [1] 91.62162
# 
# $nflg$`1m6m`$bias$Synonymous.PFitter.time.est
# [1] -11.10811
# 
# $nflg$`1m6m`$bias$multifounder.PFitter.time.est
# [1] 21.83333
# 
# $nflg$`1m6m`$bias$multifounder.Synonymous.PFitter.time.est
# [1] -25.23333
# 
# 
# $nflg$`1m6m`$rmse
# $nflg$`1m6m`$rmse$PFitter.time.est
# [1] 135.8231
# 
# $nflg$`1m6m`$rmse$Synonymous.PFitter.time.est
# [1] 32.30822
# 
# $nflg$`1m6m`$rmse$multifounder.PFitter.time.est
# [1] 70.30431
# 
# $nflg$`1m6m`$rmse$multifounder.Synonymous.PFitter.time.est
# [1] 20.79072
# 
# 
# 
# 
# $v3
# $v3$`1m`
# $v3$`1m`$bias
# $v3$`1m`$bias$PFitter.time.est
# [1] 179.2
# 
# $v3$`1m`$bias$Synonymous.PFitter.time.est
# [1] -4.2
# 
# $v3$`1m`$bias$multifounder.PFitter.time.est
# [1] 5.9
# 
# $v3$`1m`$bias$multifounder.Synonymous.PFitter.time.est
# [1] -40.9
# 
# 
# $v3$`1m`$rmse
# $v3$`1m`$rmse$PFitter.time.est
# [1] 351.5259
# 
# $v3$`1m`$rmse$Synonymous.PFitter.time.est
# [1] 113.6942
# 
# $v3$`1m`$rmse$multifounder.PFitter.time.est
# [1] 76.81481
# 
# $v3$`1m`$rmse$multifounder.Synonymous.PFitter.time.est
# [1] 35.93693
# 
# 
# 
# $v3$`6m`
# $v3$`6m`$bias
# $v3$`6m`$bias$PFitter.time.est
# [1] 139.2778
# 
# $v3$`6m`$bias$Synonymous.PFitter.time.est
# [1] -120.7778
# 
# $v3$`6m`$bias$multifounder.PFitter.time.est
# [1] -25.88889
# 
# $v3$`6m`$bias$multifounder.Synonymous.PFitter.time.est
# [1] -144.7778
# 
# 
# $v3$`6m`$rmse
# $v3$`6m`$rmse$PFitter.time.est
# [1] 203.4195
# 
# $v3$`6m`$rmse$Synonymous.PFitter.time.est
# [1] 73.17309
# 
# $v3$`6m`$rmse$multifounder.PFitter.time.est
# [1] 156.6712
# 
# $v3$`6m`$rmse$multifounder.Synonymous.PFitter.time.est
# [1] 55.60528
# 
# 
# 
# $v3$`1m6m`
# $v3$`1m6m`$bias
# $v3$`1m6m`$bias$PFitter.time.est
# [1] 259.4118
# 
# $v3$`1m6m`$bias$Synonymous.PFitter.time.est
# [1] 2.705882
# 
# $v3$`1m6m`$bias$multifounder.PFitter.time.est
# [1] 123.4118
# 
# $v3$`1m6m`$bias$multifounder.Synonymous.PFitter.time.est
# [1] -18.58824
# 
# 
# $v3$`1m6m`$rmse
# $v3$`1m6m`$rmse$PFitter.time.est
# [1] 360.1496
# 
# $v3$`1m6m`$rmse$Synonymous.PFitter.time.est
# [1] 96.10981
# 
# $v3$`1m6m`$rmse$multifounder.PFitter.time.est
# [1] 319.6836
# 
# $v3$`1m6m`$rmse$multifounder.Synonymous.PFitter.time.est
# [1] 72.51643

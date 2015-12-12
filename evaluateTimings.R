## First gather the tab files eg:
# cat /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw/nflg/1m/hiv_founder_id_processed_*/identify_founders.tab > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw/nflg/1m/identify_founders.tab 
# cat /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw/v3/1m/hiv_founder_id_processed_*/identify_founders.tab > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw/v3/1m/identify_founders.tab
## Also run the evaluateTimings.sh script eg
# ./evaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw/v3/1m/ > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw/v3/1m/sampleDates.tbl
# ./evaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw/nflg/1m/ > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw/nflg/1m/sampleDates.tbl

identify.founders.date.estimates <- c( "PFitter.time.est", "Synonymous.PFitter.time.est", "multifounder.PFitter.time.est", "multifounder.Synonymous.PFitter.time.est" ); 

rv217.gold.standard.infection.dates.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/rv217_gold_standard_timings.csv" );
rv217.gold.standard.infection.dates <- as.Date( as.character( rv217.gold.standard.infection.dates.in[,2] ), "%m/%d/%y" );
names( rv217.gold.standard.infection.dates ) <- as.character( rv217.gold.standard.infection.dates.in[,1] );

caprisa002.gold.standard.infection.dates.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/caprisa_002_gold_standard_timings.csv" );
caprisa002.gold.standard.infection.dates <- as.Date( as.character( caprisa002.gold.standard.infection.dates.in[,2] ), "%Y/%m/%d" );
names( caprisa002.gold.standard.infection.dates ) <- as.character( caprisa002.gold.standard.infection.dates.in[,1] );

studies <- c( "nflg", "v3" );
times <- c( "1m", "6m", "1m6m" );
#times <- c( "1m6m" );

timings.results.by.study.and.time <- 
 lapply( studies, function( the.study ) {
             ## TODO: REMOVE
             cat( the.study, fill = T );
     timings.results.by.time <- 
         lapply( times, function( the.time ) {
             ## TODO: REMOVE
             cat( the.time, fill = T );
             sample.dates.in <- read.delim( paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw/", the.study, "/", the.time, "/sampleDates.tbl", sep = "" ), sep = " ", header = F, fill = T );
             colnames( sample.dates.in ) <- c( "ptid", "date" );
             ## Remove anything that's not really a ptid/date combo
             sample.dates.in <- sample.dates.in[ grep( "^\\d+$", as.character( sample.dates.in[ , 1 ] ) ), , drop = FALSE ];
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

             ## Add to results: "infer" results.
             infer.results.directories <- dir( paste( "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw/", the.study, "/", the.time, sep = "" ), "founder-inference-bakeoff_", full.name = TRUE );
             infer.results.files <- sapply( infer.results.directories, dir, "outtoi.csv", full.name = TRUE );
             infer.results.list <-
                 lapply( unlist( infer.results.files ), function( .file ) {
                     .rv <- as.matrix( read.csv( .file, header = FALSE ), nrow = 1 );
                     stopifnot( ncol( .rv ) == 3 );
                     return( .rv );
                 } );
             #print( infer.results.list ); ## TODO REMOVE
             infer.results <- do.call( rbind, infer.results.list );
             colnames( infer.results ) <- c( "Infer", "Infer.CI.low", "Infer.CI.high" );
             rownames( infer.results ) <-
                 gsub( "^.+_(\\d+)$", "\\1", names( unlist( infer.results.files ) ) );
             # Special: for v3, only use caprisa seqs (not rv144, for now).
             if( the.study == "v3" ) {
                 infer.results <-
                     infer.results[ grep( "^100\\d\\d\\d", rownames( infer.results ) ), , drop = FALSE ];
             }
             # Add just the estimate from infer.
             sample.dates.char <- as.character( sample.dates.in[ , 2 ] );
             sample.dates <- as.Date( sample.dates.char );
             names( sample.dates ) <- sample.dates.in[ , 1 ];
             infer.days.before.sample <- sapply( 1:nrow( infer.results ), function( .i ) { 0 - as.numeric( as.Date( infer.results[ .i, 1 ] ) - sample.dates[ rownames( infer.results )[ .i ] ] ) } );
             names( infer.days.before.sample ) <- rownames( infer.results );

             results <- cbind( results, infer.days.before.sample );
             colnames( results )[ ncol( results ) ] <- "Infer";

             if( the.time == "1m6m" ) {
               ## Add to results: "anchre" results. (only at 1m6m)
               anchre.results.directories <- dir( paste( "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw/", the.study, "/1m6m", sep = "" ), "anchre", full.name = TRUE );
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
                           paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw/", the.study, "/1m6m/", .file.short.nosuffix, ".anc2tsv.tab", sep = "" );
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
               # Special: for v3, only use caprisa seqs (not rv144, for now).
               if( the.study == "v3" ) {
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
             } # End if 1m6m, add anchre results too.
             
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

# Make a table out of it. (one per study).
results.tables <- 
lapply( studies, function( the.study ) {
.tbl <- do.call( cbind, lapply( timings.results.by.study.and.time[[ the.study ]], function( .lst.by.time ) {
  unlist( list( Bias = .lst.by.time[[ "bias" ]], RMSE = .lst.by.time[[ "rmse" ]] ) )
} ) );
return( .tbl );
} );
names( results.tables ) <- studies;
results.tables
# $nflg
#                                                       1m         6m      1m6m
# Bias.PFitter.time.est                          43.250000  -26.36986  91.62162
# Bias.Synonymous.PFitter.time.est              -23.625000 -145.23288 -11.10811
# Bias.multifounder.PFitter.time.est              6.137255  -66.92593  21.83333
# Bias.multifounder.Synonymous.PFitter.time.est -30.843137 -151.38889 -25.23333
# Bias.Infer                                     88.531250   82.11111  72.65000
# Bias.Anchre.r2t                               120.393666  121.64161  75.26316
# Bias.Anchre.bst                                28.133187   29.92886 291.21053
# RMSE.PFitter.time.est                          61.020331   91.80863 135.82308
# RMSE.Synonymous.PFitter.time.est               16.371161   28.42296  32.30822
# RMSE.multifounder.PFitter.time.est              4.607038   12.27091  70.30431
# RMSE.multifounder.Synonymous.PFitter.time.est  43.250000  -26.36986  20.79072
# RMSE.Infer                                    -23.625000 -145.23288  24.31433
# RMSE.Anchre.r2t                                 6.137255  -66.92593 174.86371
# RMSE.Anchre.bst                               -30.843137 -151.38889 714.73403
# 
# $v3
#                                                       1m         6m
# Bias.PFitter.time.est                         179.200000  139.27778
# Bias.Synonymous.PFitter.time.est               -4.200000 -120.77778
# Bias.multifounder.PFitter.time.est              5.900000  -25.88889
# Bias.multifounder.Synonymous.PFitter.time.est -40.900000 -144.77778
# Bias.Infer                                     77.450000   93.38889
# Bias.Anchre.r2t                               351.525937  203.41953
# Bias.Anchre.bst                               113.694234   73.17309
# RMSE.PFitter.time.est                          76.814815  156.67121
# RMSE.Synonymous.PFitter.time.est               35.936933   55.60528
# RMSE.multifounder.PFitter.time.est              5.296225   12.71469
# RMSE.multifounder.Synonymous.PFitter.time.est 179.200000  139.27778
# RMSE.Infer                                     -4.200000 -120.77778
# RMSE.Anchre.r2t                                 5.900000  -25.88889
# RMSE.Anchre.bst                               -40.900000 -144.77778
#                                                        1m6m
# Bias.PFitter.time.est                            259.411765
# Bias.Synonymous.PFitter.time.est                   2.705882
# Bias.multifounder.PFitter.time.est               123.411765
# Bias.multifounder.Synonymous.PFitter.time.est    -18.588235
# Bias.Infer                                       101.352941
# Bias.Anchre.r2t                                  695.529412
# Bias.Anchre.bst                               321256.235294
# RMSE.PFitter.time.est                            360.149632
# RMSE.Synonymous.PFitter.time.est                  96.109810
# RMSE.multifounder.PFitter.time.est               319.683644
# RMSE.multifounder.Synonymous.PFitter.time.est     72.516428
# RMSE.Infer                                        24.264020
# RMSE.Anchre.r2t                                 1604.720736
# RMSE.Anchre.bst                               333778.988566

# ------- FOR THE RECORD ---------
# timings.results.by.study.and.time
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
# $nflg$`1m`$bias$Infer
# [1] 88.53125
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
# $nflg$`1m`$rmse$Infer
# [1] 4.607038
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
# $nflg$`6m`$bias$Infer
# [1] 82.11111
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
# $nflg$`6m`$rmse$Infer
# [1] 12.27091
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
# $nflg$`1m6m`$bias$Infer
# [1] 76.05
# 
# $nflg$`1m6m`$bias$Anchre.r2t
# [1] 12
# 
# $nflg$`1m6m`$bias$Anchre.bst
# [1] 100.7059
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
# $nflg$`1m6m`$rmse$Infer
# [1] 21.87157
# 
# $nflg$`1m6m`$rmse$Anchre.r2t
# [1] 32.39792
# 
# $nflg$`1m6m`$rmse$Anchre.bst
# [1] 323.2957
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
# $v3$`1m`$bias$Infer
# [1] 77.45
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
# $v3$`1m`$rmse$Infer
# [1] 5.296225
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
# $v3$`6m`$bias$Infer
# [1] 93.38889
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
# $v3$`6m`$rmse$Infer
# [1] 12.71469
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
# $v3$`1m6m`$bias$Infer
# [1] 101.3529
# 
# $v3$`1m6m`$bias$Anchre.r2t
# [1] 695.5294
# 
# $v3$`1m6m`$bias$Anchre.bst
# [1] 321256.2
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
# 
# $v3$`1m6m`$rmse$Infer
# [1] 24.26402
# 
# $v3$`1m6m`$rmse$Anchre.r2t
# [1] 1604.721
# 
# $v3$`1m6m`$rmse$Anchre.bst
# [1] 333779

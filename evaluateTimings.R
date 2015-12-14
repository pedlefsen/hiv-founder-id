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
             results.in <- read.delim( paste( paste( "/fh/fast/edlefsen_p/bakeoff_merged_analysis_sequences/raw_fixed/", the.study, "/", the.time, "/identify_founders.tab", sep = "" ) ), sep = "\t" );
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
             infer.results.directories <- dir( paste( "/fh/fast/edlefsen_p/bakeoff_merged_analysis_sequences/raw_fixed/", the.study, "/", the.time, sep = "" ), "founder-inference-bakeoff_", full.name = TRUE );
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
                           paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw_fixed/", the.study, "/1m6m/", .file.short.nosuffix, ".anc2tsv.tab", sep = "" );
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
    lapply( names( timings.results.by.study.and.time ), function( the.study ) {
        .rv <- 
        lapply( names( timings.results.by.study.and.time[[ the.study ]] ), function( the.time ) { sapply( timings.results.by.study.and.time[[ the.study ]][[ the.time ]], function( results.list ) { results.list } ) } );
        names( .rv ) <- times;
        return( .rv );
    } );
names( results.tables ) <- studies;
results.tables
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
# Infer                                    77.45  5.296225
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
# Infer                                    96.57895  28.13127
# Anchre.r2t                               620.6842  1530.37 
# Anchre.bst                               325964    337789.5


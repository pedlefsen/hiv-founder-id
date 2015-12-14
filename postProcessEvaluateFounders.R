### NOTE: This is assuming you have _already_ compiled the evaluateFounders.tbl results into a single file at the top-level dir (per time point, region): eg /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_fixed/nflg/1m/evaluateFounders.tbl
##  NOTE continued: for now I'm manually creating these files because I am lazy and in one motion within emacs I can get rid of the extra table headers.
# eg cat /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_fixed/nflg/1m/*/evaluateFounders.tbl > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_fixed/nflg/1m/evaluateFounders.tbl
# then open and remove the extra instances of the first (header) line.

# From this file we read in indicators of whether to use the multiple- or single-founder true profile.
rv217.gold.standards.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/RV217_gold_standards.csv" );
rv217.gold.is.multiple <- rv217.gold.standards.in[ , "gold.is.multiple" ];
names( rv217.gold.is.multiple ) <- rv217.gold.standards.in[ , "ptid" ];

# From this file we read in indicators of whether to use the multiple- or single-founder true profile.
caprisa002.gold.standards.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/caprisa_002_gold_standards.csv" );
caprisa002.gold.is.multiple <- caprisa002.gold.standards.in[ , "gold.is.multiple" ];
names( caprisa002.gold.is.multiple ) <- caprisa002.gold.standards.in[ , "ptid" ];

## Read in the evaluateFounders results for the identify-founders method, aggregate and return a matrix of results with ptids in rows.
postProcessEvaluateFounders <- function ( the.study, the.time, use.infer = FALSE ) {
    if( use.infer ) {
        maybe.subdir <- "/infer";
    } else {
        maybe.subdir <- "";
    }
    if( the.study == "nflg" ) {
        the.study.alt <- "rv217";
    } else {
        the.study.alt <- "caprisa002";
    }
    results.in <- read.delim( paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw_fixed/", the.study, maybe.subdir, "/", the.time, "/evaluateFounders.tbl", sep = "" ), sep = "\t" );
    results.in.estimates.ptid <-
        gsub( paste( ".*", the.study.alt, "_([^_]+)_.+", sep = "" ), "\\1", as.character( results.in[ , "estimates.file" ] ) );
    results.in.truths.ptid <-
        gsub( paste( ".*", the.study.alt, "_([^_]+)_.+", sep = "" ), "\\1", as.character( results.in[ , "truths.file" ] ) );
    if( !all( results.in.estimates.ptid == results.in.truths.ptid ) ) {
      warning( paste( "ptids mismatch?:", paste( results.in.estimates.ptid, collapse = ", " ), "not same as", paste( results.in.truths.ptid, collapse = ", " )  ) );
    }
    results.in.ptids <-
        unique( results.in.truths.ptid );
    
    ## Note that we do not subset to same-region results, since the regions overlap, we can aggregate HDs (ignoring gaps) by summing numerator and denominators.  This solidifies my thinking that the one that includes gaps does not make sense here, and gaps might be inserted by geneCutter to make codons work as reasonable AAs in translation, which is really no fault of the estimation, as it might be context-dependent (eg if a substitution in one homolgous pair results in that pair being not-aligned after codon-alignment because the placement leads to different best-guesses of the functional codon [though my thinking is that this should be exceedingly rare]).
    
    results.fromto <- results.in[ , 1:2, drop = FALSE ];
    
    ## Here we do the sanity check to ensure that every ( from, to ) pair appears at most once.
    .fromto <- apply( results.fromto, 1, paste, collapse = "---" );
    stopifnot( length( .fromto ) == length( unique( .fromto ) ) );
    
    if( the.study == "v3" ) {
        gold.is.multiple <- caprisa002.gold.is.multiple;
    } else {
        gold.is.multiple <- rv217.gold.is.multiple;
    }
    
    estimates.file.ptid <-
        gsub( paste( ".*", the.study.alt, "_([^_]+)_.+", sep = "" ), "\\1", as.character( results.fromto[ , "estimates.file" ] ) );
    truths.file.ptid <-
        gsub( paste( ".*", the.study.alt, "_([^_]+)_.+", sep = "" ), "\\1", as.character( results.fromto[ , "truths.file" ] ) );
    stopifnot( all( truths.file.ptid == estimates.file.ptid ) );
    stopifnot( all( truths.file.ptid %in% names( gold.is.multiple ) ) );
    
    estimates.file.is.multiple <-
        sapply( as.character( results.fromto[ , "estimates.file" ] ), function( .str ) { return( length( grep( "single", .str, invert = TRUE ) ) > 0 ) } );
    truths.file.is.multiple <-
        sapply( as.character( results.fromto[ , "truths.file" ] ), function( .str ) { return( length( grep( "single", .str, invert = TRUE ) ) > 0 ) } );
    
    ## Ok here is where we subset to use just the "single" or "multiple" results, depending on the value of gold.is.multiple.  Note that we prefer to use the matching result, but if it is not there our next best result is the one with the matching "truth" is.multiple.  That is, we can compare a multi-founder estimate to a single-founder truth, if there is only a multi-founder estimate.
    truths.matches.gold.is.multiple <-
        ( truths.file.is.multiple == gold.is.multiple[ truths.file.ptid ] );
    estimates.matches.gold.is.multiple <- ( estimates.file.is.multiple == gold.is.multiple[ truths.file.ptid ] );
    
    ### Some ptids are missing the gold-matching pair, so we choose a different match-state for these, preferring truths-matches over estimates-matches.
    .ptids.with.matching.estimates <- unique( estimates.file.ptid[ ( truths.matches.gold.is.multiple & estimates.matches.gold.is.multiple ) ] );
    .missing.ptids <- setdiff( results.in.ptids, .ptids.with.matching.estimates );
    .ptids.with.mismatching.estimates <- intersect( .missing.ptids, estimates.file.ptid[ ( truths.matches.gold.is.multiple & !estimates.matches.gold.is.multiple ) ] );
    .missing.ptids <- setdiff( .missing.ptids, .ptids.with.mismatching.estimates );
    .ptids.with.mismatching.truths <- intersect( .missing.ptids, estimates.file.ptid[ ( !truths.matches.gold.is.multiple & estimates.matches.gold.is.multiple ) ] );
    .missing.ptids <- setdiff( .missing.ptids, .ptids.with.mismatching.estimates );
    .ptids.with.mismatching.both <- intersect( .missing.ptids, estimates.file.ptid[ ( !truths.matches.gold.is.multiple & !estimates.matches.gold.is.multiple ) ] );
    .missing.ptids <- setdiff( .missing.ptids, .ptids.with.mismatching.estimates );
    ## All ptids should be in one of these groups (not missing the gold-matching valu of is.multiple; or one of the .ptids.with.missing* -- note that a ptid would possibly qualify for multiple groups but we have a hierarchy, prefering matches to truth, for instance).
    #stopifnot( length( .missing.ptids ) == 0 );
    if( length( .missing.ptids ) > 0 ) {
      warning( paste( "After filtering, we are missing ptids:", paste( .missing.ptids, collapse = ", " ) ) );
    }
    
    results.matches.gold.is.multiple <-
        results.in[ ( truths.matches.gold.is.multiple & estimates.matches.gold.is.multiple ), , drop = FALSE ];
    results.matches.gold.is.multiple.truths.ptid <-
        truths.file.ptid[ ( truths.matches.gold.is.multiple & estimates.matches.gold.is.multiple ) ];
    results.matches.gold.is.multiple.in.truths.but.not.estimates <-
        results.in[ ( results.in.truths.ptid %in% .ptids.with.mismatching.estimates ) & ( truths.matches.gold.is.multiple & !estimates.matches.gold.is.multiple ), , drop = FALSE ];
    results.matches.gold.is.multiple.in.truths.but.not.estimates.truths.ptid <- 
        truths.file.ptid[ ( results.in.truths.ptid %in% .ptids.with.mismatching.estimates ) & ( truths.matches.gold.is.multiple & !estimates.matches.gold.is.multiple ) ];
    results.matches.gold.is.multiple.in.estimates.but.not.truths <-
        results.in[ ( results.in.truths.ptid %in% .ptids.with.mismatching.truths ) & ( !truths.matches.gold.is.multiple & estimates.matches.gold.is.multiple ), , drop = FALSE ];
    results.matches.gold.is.multiple.in.estimates.but.not.truths.truths.ptid <- 
        truths.file.ptid[ ( results.in.truths.ptid %in% .ptids.with.mismatching.truths ) & ( !truths.matches.gold.is.multiple & estimates.matches.gold.is.multiple ) ];
    results.doesnt.match.gold.is.multiple.in.either <-
        results.in[ ( results.in.truths.ptid %in% .ptids.with.mismatching.both ) & ( !truths.matches.gold.is.multiple & !estimates.matches.gold.is.multiple ), , drop = FALSE ];
    results.doesnt.match.gold.is.multiple.in.either.truths.ptid <- 
        truths.file.ptid[ ( results.in.truths.ptid %in% .ptids.with.mismatching.both ) & ( !truths.matches.gold.is.multiple & !estimates.matches.gold.is.multiple ) ];
    
    ## The ptids won't have changed but their order may have.
    results.ptids <- c( results.matches.gold.is.multiple.truths.ptid, results.matches.gold.is.multiple.in.truths.but.not.estimates.truths.ptid, results.matches.gold.is.multiple.in.estimates.but.not.truths.truths.ptid, results.doesnt.match.gold.is.multiple.in.either.truths.ptid );
    stopifnot( length( unique( results.ptids ) ) == length( results.in.ptids ) );
        
    results.prep <- as.matrix( rbind( results.matches.gold.is.multiple, results.matches.gold.is.multiple.in.truths.but.not.estimates, results.matches.gold.is.multiple.in.estimates.but.not.truths, results.doesnt.match.gold.is.multiple.in.either ) )[ , 3:ncol( results.matches.gold.is.multiple ) ];
    suppressWarnings( mode( results.prep ) <- "numeric" );
    
    # ## TODO: Read this in, use all the "fits" and "is.starlike" etc options to choose whether to use single or multiple in the estimated set.  NOTE: FOR NOW we just use the "right" one for the job, ie the same one we use for the "truth".
    # identifyfounders.in <- read.delim( "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw_fixed/nflg/1m/identify_founders.tab", sep = "\t" );
    # identifyfounders <- as.matrix( identifyfounders.in );
    # suppressWarnings( mode( identifyfounders ) <- "numeric" );

    ## We remove the "includingGaps"; keeping only "ignoringGaps".
    ## We remove the "FULL_SEQUENCE" result if it is there, and we compute our own "nearly full length genome" result.
    cols.to.keep <-
        grep( "FULL_SEQUENCE", grep( "\\.sum\\.ignoringGaps", colnames( results.prep ), value = TRUE ), invert = TRUE, value = TRUE );
    cols.to.keep.prefixes <- gsub( "^(.+)\\.perspective.+", "\\1", cols.to.keep );
    results.prep2 <- results.prep[ , cols.to.keep, drop = FALSE ];
    rownames( results.prep2 ) <- results.ptids;
    
    ## Here we add in new aggregation over the whole genome (for multi-gene results only; for now it applies to nflg not v3).
    if( the.study == "nflg" ) {
        # gather columns with like suffixes.
        colnames.sans.gene <- gsub( "^[^\\.]+(\\..+)$", "\\1", colnames( results.prep2 ) );
        unique.colnames.sans.gene <- unique( colnames.sans.gene );
        results.prep.NFLG.list <- lapply( unique.colnames.sans.gene, function( .colname.sans.gene ) {
            apply( results.prep2[ , ( colnames.sans.gene == .colname.sans.gene ), drop = FALSE ], 1, sum, na.rm = T );
        } );
        names( results.prep.NFLG.list ) <- paste( "NFLG", unique.colnames.sans.gene, sep = "" );
        results.prep.NFLG <- do.call( cbind, results.prep.NFLG.list );
        results.prep.withNFLG <-
            cbind( results.prep.NFLG, results.prep );

        # Adopt this.
        results.prep2 <- results.prep.withNFLG;
        cols.to.keep <-
            grep( "FULL_SEQUENCE", grep( "\\.sum\\.ignoringGaps", colnames( results.prep.withNFLG ), value = TRUE ), invert = TRUE, value = TRUE );
        cols.to.keep.prefixes <- gsub( "^(.+)\\.perspective.+", "\\1", cols.to.keep );
    }
        
    ## (in addition to the includingGaps results; see above), we remove the interesting (but not for reporting) perspectives results, which are components of the result that we use: truths.perspective is the average over true founders of the HD of the smallest-HD estimate for that true founder (even if it is the same as another one).
    ## Ok so now we need to recompute per-ptid averages over the "truths" and "esimtates" 'perspectives'; first we sum like columns for a ptid.  There may be multiple rows for the rv217 "nflg" data (LH, RH, env, NFLG) but not for caprisa002 "v3".
    results.prep3 <- t(
        sapply( results.ptids, function( .ptid ) {
            .ptid.sums <- 
                apply( results.prep2[ results.ptids == .ptid, , drop = FALSE ], 2, sum, na.rm = T );
            .rv <- sapply( unique( cols.to.keep.prefixes ), function( .prefix ) {
                .ptid.sums[ paste( .prefix, "perspective.HD.sum.ignoringGaps", sep = "." ) ] /
                .ptid.sums[ paste( .prefix, "perspective.denominator.sum.ignoringGaps", sep = "." ) ]
            } );
            names( .rv ) <- unique( cols.to.keep.prefixes );
            return( .rv );
        } ) );
    
    ## And finally replace it all with averages (we average the truths-perspective and the estimates-perspective).
    results.prep3.prefixes <- gsub( "^(.+)\\.(truths|estimates)", "\\1", colnames( results.prep3 ) );
    
    results <- sapply( unique( results.prep3.prefixes ), function( .prefix ) {
        apply( results.prep3[ , paste( .prefix, c( "truths", "estimates" ), sep = "." ), drop = FALSE ], 1, mean, na.rm = T )
    } );
    
    return( results );
} # postProcessEvaluateFounders (..)

# These have only two results: V3 for AA and for NA.
caprisa002.v3.1m.results <- postProcessEvaluateFounders( "v3", "1m" );
caprisa002.v3.1m.results.infer <- postProcessEvaluateFounders( "v3", "1m", use.infer = TRUE );
caprisa002.v3.6m.results <- postProcessEvaluateFounders( "v3", "6m" );
caprisa002.v3.6m.results.infer <- postProcessEvaluateFounders( "v3", "6m", use.infer = TRUE );
caprisa002.v3.1m6m.results <- postProcessEvaluateFounders( "v3", "1m6m" );
caprisa002.v3.1m6m.results.infer <- postProcessEvaluateFounders( "v3", "1m6m", use.infer = TRUE );

pdf( file = "caprisa002.v3.1m.results.pdf" )
boxplot( caprisa002.v3.1m.results[,1], caprisa002.v3.1m.results[,2 ], names = colnames( caprisa002.v3.1m.results ) )
dev.off()

pdf( file = "caprisa002.v3.1m.results.infer.pdf" )
boxplot( caprisa002.v3.1m.results.infer[,1], caprisa002.v3.1m.results.infer[,2 ], names = colnames( caprisa002.v3.1m.results.infer ) )
dev.off()

pdf( file = "caprisa002.v3.6m.results.pdf" )
boxplot( caprisa002.v3.6m.results[,1], caprisa002.v3.6m.results[,2 ], names = colnames( caprisa002.v3.6m.results ) )
dev.off()

pdf( file = "caprisa002.v3.6m.results.infer.pdf" )
boxplot( caprisa002.v3.6m.results.infer[,1], caprisa002.v3.6m.results.infer[,2 ], names = colnames( caprisa002.v3.6m.results.infer ) )
dev.off()

pdf( file = "caprisa002.v3.1m6m.results.pdf" )
boxplot( caprisa002.v3.1m6m.results[,1], caprisa002.v3.1m6m.results[,2 ], names = colnames( caprisa002.v3.1m6m.results ) )
dev.off()

pdf( file = "caprisa002.v3.1m6m.results.infer.pdf" )
boxplot( caprisa002.v3.1m6m.results.infer[,1], caprisa002.v3.1m6m.results.infer[,2 ], names = colnames( caprisa002.v3.1m6m.results.infer ) )
dev.off()

pdf( file = "caprisa002.v3.NA.over.time.pdf" )
boxplot( caprisa002.v3.1m.results[,1 ], caprisa002.v3.6m.results[,1 ], caprisa002.v3.1m6m.results[,1 ], names = c( "1m.V3.NA", "6m.V3.NA", "1m6m.V3.NA" ) )
dev.off()

pdf( file = "caprisa002.v3.AA.over.time.pdf" )
boxplot( caprisa002.v3.1m.results[,2 ], caprisa002.v3.6m.results[,2 ], caprisa002.v3.1m6m.results[,2 ], names = c( "1m.V3.AA", "6m.V3.AA", "1m6m.V3.AA" ) )
dev.off()

## For like subjects, take diffs
.shared.ptids.6m <- intersect( rownames( caprisa002.v3.1m.results ), rownames( caprisa002.v3.6m.results ) );
.shared.ptids.1m6m <- intersect( rownames( caprisa002.v3.1m.results ), rownames( caprisa002.v3.1m6m.results ) );
pdf( file = "caprisa002.v3.NA.over.time.within.subjects.pdf" )
boxplot( caprisa002.v3.6m.results[.shared.ptids.6m,1 ] - caprisa002.v3.1m.results[.shared.ptids.6m,1 ], caprisa002.v3.1m6m.results[.shared.ptids.1m6m,1 ] - caprisa002.v3.1m.results[.shared.ptids.1m6m,1 ], names = c( "6m.V3.AA-1m.V3.AA", "1m6m.V3.AA-1m.V3.AA" ) )
dev.off()
.shared.ptids <- intersect( rownames( caprisa002.v3.1m.results ), rownames( caprisa002.v3.6m.results ) );
pdf( file = "caprisa002.v3.AA.over.time.within.subjects.pdf" )
boxplot( caprisa002.v3.6m.results[.shared.ptids,2 ] - caprisa002.v3.1m.results[.shared.ptids,2 ], caprisa002.v3.1m6m.results[.shared.ptids,2 ] - caprisa002.v3.1m.results[.shared.ptids,2 ], names = c( "6m.V3.AA-1m.V3.AA", "1m6m.V3.AA-1m.V3.AA" ) )
dev.off()


## RV217
rv217.nflg.1m.results <- postProcessEvaluateFounders( "nflg", "1m" );
rv217.nflg.1m.results.infer <- postProcessEvaluateFounders( "nflg", "1m", use.infer = TRUE );
rv217.nflg.6m.results <- postProcessEvaluateFounders( "nflg", "6m" );
rv217.nflg.6m.results.infer <- postProcessEvaluateFounders( "nflg", "6m", use.infer = TRUE );
rv217.nflg.1m6m.results <- postProcessEvaluateFounders( "nflg", "1m6m" );
rv217.nflg.1m6m.results.infer <- postProcessEvaluateFounders( "nflg", "1m6m", use.infer = TRUE );

pdf( file = "rv217.nflg.1m.NA.methods.scatter" )
boxplot( c( rv217.nflg.1m.results.infer[, 1:9 ] ), c( rv217.nflg.1m.results[, 1:9 ] ), names = c( "nflg.1m.NA.mean", "nflg.1m.NA.infer.mean" ) )
dev.off()
pdf( file = "rv217.nflg.1m.AA.methods.scatter" )
boxplot( c( rv217.nflg.1m.results.infer[, 10:18 ] ), c( rv217.nflg.1m.results[, 10:18 ] ), names = c( "nflg.1m.AA.infer.mean", "nflg.1m.AA.mean" ) )
dev.off()

pdf( file = "rv217.nflg.1m.NA.methods" )
boxplot( apply( rv217.nflg.1m.results.infer[, 1:9 ], 1, mean, na.rm = T ), apply( rv217.nflg.1m.results[, 1:9 ], 1, mean, na.rm = T ), names = c( "nflg.1m.NA.infer.mean", "nflg.1m.NA.mean" ) )
dev.off()
t.test( apply( rv217.nflg.1m.results.infer[, 10:18 ], 1, mean, na.rm = T ), apply( rv217.nflg.1m.results[, 10:18 ], 1, mean, na.rm = T ) )
pdf( file = "rv217.nflg.1m.AA.methods" )
boxplot( apply( rv217.nflg.1m.results.infer[, 10:18 ], 1, mean, na.rm = T ), apply( rv217.nflg.1m.results[, 10:18 ], 1, mean, na.rm = T ), names = c( "nflg.1m.AA.infer.mean", "nflg.1m.AA.mean" ) )
dev.off()
t.test( apply( rv217.nflg.1m.results.infer[, 10:18 ], 1, mean, na.rm = T ), apply( rv217.nflg.1m.results[, 10:18 ], 1, mean, na.rm = T ) )


pdf( file = "rv217.nflg.6m.NA.methods.scatter" )
boxplot( c( rv217.nflg.6m.results.infer[, 1:9 ] ), c( rv217.nflg.6m.results[, 1:9 ] ), names = c( "nflg.6m.NA.mean", "nflg.6m.NA.infer.mean" ) )
dev.off()
pdf( file = "rv217.nflg.6m.AA.methods.scatter" )
boxplot( c( rv217.nflg.6m.results.infer[, 10:18 ] ), c( rv217.nflg.6m.results[, 10:18 ] ), names = c( "nflg.6m.AA.infer.mean", "nflg.6m.AA.mean" ) )
dev.off()

pdf( file = "rv217.nflg.6m.NA.methods" )
boxplot( apply( rv217.nflg.6m.results.infer[, 1:9 ], 1, mean, na.rm = T ), apply( rv217.nflg.6m.results[, 1:9 ], 1, mean, na.rm = T ), names = c( "nflg.6m.NA.infer.mean", "nflg.6m.NA.mean" ) )
dev.off()
t.test( apply( rv217.nflg.6m.results.infer[, 10:18 ], 1, mean, na.rm = T ), apply( rv217.nflg.6m.results[, 10:18 ], 1, mean, na.rm = T ) )
pdf( file = "rv217.nflg.6m.AA.methods" )
boxplot( apply( rv217.nflg.6m.results.infer[, 10:18 ], 1, mean, na.rm = T ), apply( rv217.nflg.6m.results[, 10:18 ], 1, mean, na.rm = T ), names = c( "nflg.6m.AA.infer.mean", "nflg.6m.AA.mean" ) )
dev.off()
t.test( apply( rv217.nflg.6m.results.infer[, 10:18 ], 1, mean, na.rm = T ), apply( rv217.nflg.6m.results[, 10:18 ], 1, mean, na.rm = T ) )

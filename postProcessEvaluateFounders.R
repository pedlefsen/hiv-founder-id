### NOTE: This is assuming you have _already_ done everything in the README.postprocessing.txt file.

source( "createBoxplotShowingSignificance_safetosource.R" )
source( "prettyPrintPValues_safetosource.R" )

GOLD.STANDARD.DIR <- "/fh/fast/edlefsen_p/bakeoff/gold_standard";
#RESULTS.DIR <- "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw_fixed/";
#RESULTS.DIR <- "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/";
RESULTS.DIR <- "/fh/fast/edlefsen_p/bakeoff_merged_analysis_sequences_results/raw_fixed/";

# From this file we read in indicators of whether to use the multiple- or single-founder true profile.
rv217.gold.standards.in <- read.csv( paste( GOLD.STANDARD.DIR, "rv217/RV217_gold_standards.csv", sep = "/" ) );
rv217.gold.is.multiple <- rv217.gold.standards.in[ , "gold.is.multiple" ];
names( rv217.gold.is.multiple ) <- rv217.gold.standards.in[ , "ptid" ];

# From this file we read in indicators of whether to use the multiple- or single-founder true profile.
caprisa002.gold.standards.in <- read.csv( paste( GOLD.STANDARD.DIR, "caprisa_002/caprisa_002_gold_standards.csv", sep = "/" ) );
caprisa002.gold.is.multiple <- caprisa002.gold.standards.in[ , "gold.is.multiple" ];
names( caprisa002.gold.is.multiple ) <- caprisa002.gold.standards.in[ , "ptid" ];

evaluateResultsMatrix <- function( results.mat, the.name = "Result", the.title = NULL ) {
    .stats <- sapply( colnames( results.mat )[ 1:2 ], function( .colname ) {
        .data <- results.mat[ , .colname ];
        c( n = sum( !is.na( results.mat[ , .colname ] ) ), summary( .data ), sd = sd( .data, na.rm = TRUE ), frac.zero = mean( results.mat[ , .colname ] == 0, na.rm = T ) )
      } );
    # The first two are the main ones; the remainder if there are a breakdown by region.
    .name.A <- gsub( "NA", "RNA", gsub( "^[^\\.]+\\.", "", colnames( results.mat )[ 1 ] ) );
    .name.B <- gsub( "NA", "RNA", gsub( "^[^\\.]+\\.", "", colnames( results.mat )[ 2 ] ) );
    .var.name <- paste( "Mean Hamming distance of", the.name, "to the true founder(s)" );
    # Note that we don't actually show significance in the sense that we are not comparing the boxplots to each other.

    ## TODO: Also create plots for the other (non-first-two-column) results
    return( list( ggp = createBoxplotShowingSignificance( as.factor( c( rep( .name.A, nrow( results.mat ) ), rep( .name.B, nrow( results.mat ) ) ) ), c( results.mat[ , 1 ], results.mat[ , 2 ] ), the.var.name = .var.name, test.results.as.list.of.lists = list( .name.A = list( .name.B = NA  ) ), the.title = the.title ), stats = .stats ) );
} # evaluateResultsMatrix (..)

compareResultsMatrix <- function( results.mat.A, results.mat.B, name.A = "A", name.B = "B", the.title = NULL ) {
    stopifnot( ncol( results.mat.A ) == ncol( results.mat.B ) );
    .shared.ptids <- intersect( rownames( results.mat.A ), rownames( results.mat.B ) );
    .unique.ptids.A <- setdiff( .shared.ptids, rownames( results.mat.A ) );
    .unique.ptids.B <- setdiff( .shared.ptids, rownames( results.mat.B ) );
    ## TODO: Deal with these -- can we add them in another color to the boxplot?
    if( length( .unique.ptids.A ) > 0 ) {
        cat( "NOTE: There are", length( .unique.ptids.A ), "ptids unique to A", fill = TRUE );
    }
    if( length( .unique.ptids.B ) > 0 ) {
        cat( "NOTE: There are", length( .unique.ptids.B ), "ptids unique to B", fill = TRUE );
    }
    # 
    .within.ptid.diffs <- results.mat.A[ .shared.ptids, ] - results.mat.B[ .shared.ptids, ];
    .within.ptid.diffs.t.test.results <- apply( .within.ptid.diffs, 2, t.test );
    .within.ptid.diffs.t.test.results.p.values <- lapply( .within.ptid.diffs.t.test.results, function( .result ) { .result$p.value } )
    # The first two are the main ones; the remainder if there are a breakdown by region.
    .name.A <- paste( gsub( "NA", "RNA", gsub( "^[^\\.]+\\.", "", colnames( results.mat.A )[ 1 ] ) ), paste( "(p", prettyPrintPValuesTo4Digits( .within.ptid.diffs.t.test.results.p.values[[ 1 ]], " = ", " <= " ), ")", sep = "" ), sep = "\n" );
    .name.B <- paste( gsub( "NA", "RNA", gsub( "^[^\\.]+\\.", "", colnames( results.mat.A )[ 2 ] ) ), paste( "(p", prettyPrintPValuesTo4Digits( .within.ptid.diffs.t.test.results.p.values[[ 2 ]], " = ", " <= " ), ")", sep = "" ), sep = "\n" );
    .var.name <- paste( "Within-participant difference", paste( "(", name.A, " - ", name.B, ")", sep = "" ) );
    # Note that we don't actually show significance in the sense that we are not comparing the boxplots to each other.

    ## TODO: Also create plots for the other (non-first-two-column) results
    return( list( ggp = createBoxplotShowingSignificance( as.factor( c( rep( .name.A, nrow( .within.ptid.diffs ) ), rep( .name.B, nrow( .within.ptid.diffs ) ) ) ), c( .within.ptid.diffs[ , 1 ], .within.ptid.diffs[ , 2 ] ), the.var.name = .var.name, test.results.as.list.of.lists = list( .name.A = list( .name.B = NA  ) ), the.title = the.title ), p.values = unlist( .within.ptid.diffs.t.test.results.p.values ) ) );
} # compareResultsMatrix (..)


## Read in the evaluateFounders results for the identify-founders method, aggregate and return a matrix of results with ptids in rows.
postProcessEvaluateFounders <- function ( the.region, the.time, use.infer = FALSE ) {
    if( use.infer ) {
        maybe.subdir <- "/infer";
    } else {
        maybe.subdir <- "";
    }
    if( ( the.region == "nflg" ) || ( the.region == "rv217_v3" ) ) {
        the.study <- "rv217";
    } else {
        stopifnot( the.region == "v3" );
        the.study <- "caprisa002";
    }
    results.in <- read.delim( paste( RESULTS.DIR, the.region, "/", the.time, maybe.subdir, "/evaluateFounders.tbl", sep = "" ), sep = "\t" );
    results.in.estimates.ptid <-
        gsub( paste( ".*", the.study, "_([^_]+)_.+", sep = "" ), "\\1", as.character( results.in[ , "estimates.file" ] ) );
    results.in.truths.ptid <-
        gsub( paste( ".*", the.study, "_([^_]+)_.+", sep = "" ), "\\1", as.character( results.in[ , "truths.file" ] ) );
    if( !all( results.in.estimates.ptid == results.in.truths.ptid ) ) {
      warning( paste( "ptids mismatch?:", paste( results.in.estimates.ptid, collapse = ", " ), "not same as", paste( results.in.truths.ptid, collapse = ", " )  ) );
    }
    results.in.ptids <-
        unique( results.in.truths.ptid );
    
    ## Note that we do not subset to same-region results, since the regions do not overlap, we can aggregate HDs (ignoring gaps) by summing numerator and denominators.  This solidifies my thinking that the one that includes gaps does not make sense here, and gaps might be inserted by geneCutter to make codons work as reasonable AAs in translation, which is really no fault of the estimation, as it might be context-dependent (eg if a substitution in one homolgous pair results in that pair being not-aligned after codon-alignment because the placement leads to different best-guesses of the functional codon [though my thinking is that this should be exceedingly rare]).
    
    results.fromto <- results.in[ , 1:2, drop = FALSE ];
    
    ## Here we do the sanity check to ensure that every ( from, to ) pair appears at most once.
    .fromto <- apply( results.fromto, 1, paste, collapse = "---" );
    stopifnot( length( .fromto ) == length( unique( .fromto ) ) );
    
    if( the.region == "v3" ) {
        gold.is.multiple <- caprisa002.gold.is.multiple;
    } else {
        gold.is.multiple <- rv217.gold.is.multiple;
    }
    
    estimates.file.ptid <-
        gsub( paste( ".*", the.study, "_([^_]+)_.+", sep = "" ), "\\1", as.character( results.fromto[ , "estimates.file" ] ) );
    truths.file.ptid <-
        gsub( paste( ".*", the.study, "_([^_]+)_.+", sep = "" ), "\\1", as.character( results.fromto[ , "truths.file" ] ) );
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
    if( the.region == "nflg" ) {
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
        sapply( unique( results.ptids ), function( .ptid ) {
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
out.file <- paste( RESULTS.DIR, "v3/1m/postProcessEvaluateFounders.tab", sep = "" );
write.table( apply( caprisa002.v3.1m.results, 1:2, function( .x ) { sprintf( "%0.4f", .x ) } ), quote = FALSE, file = out.file, sep = "\t" );
caprisa002.v3.1m.results.infer <- postProcessEvaluateFounders( "v3", "1m", use.infer = TRUE );
out.file <- paste( RESULTS.DIR, "v3/1m/infer/postProcessEvaluateFounders.tab", sep = "" );
write.table( apply( caprisa002.v3.1m.results.infer, 1:2, function( .x ) { sprintf( "%0.4f", .x ) } ), quote = FALSE, file = out.file, sep = "\t" );

caprisa002.v3.6m.results <- postProcessEvaluateFounders( "v3", "6m" );
out.file <- paste( RESULTS.DIR, "v3/6m/postProcessEvaluateFounders.tab", sep = "" );
write.table( apply( caprisa002.v3.6m.results, 1:2, function( .x ) { sprintf( "%0.4f", .x ) } ), quote = FALSE, file = out.file, sep = "\t" );
caprisa002.v3.6m.results.infer <- postProcessEvaluateFounders( "v3", "6m", use.infer = TRUE );
out.file <- paste( RESULTS.DIR, "v3/6m/infer/postProcessEvaluateFounders.tab", sep = "" );
write.table( apply( caprisa002.v3.6m.results.infer, 1:2, function( .x ) { sprintf( "%0.4f", .x ) } ), quote = FALSE, file = out.file, sep = "\t" );

caprisa002.v3.1m6m.results <- postProcessEvaluateFounders( "v3", "1m6m" );
out.file <- paste( RESULTS.DIR, "v3/1m6m/postProcessEvaluateFounders.tab", sep = "" );
write.table( apply( caprisa002.v3.1m6m.results, 1:2, function( .x ) { sprintf( "%0.4f", .x ) } ), quote = FALSE, file = out.file, sep = "\t" );
caprisa002.v3.1m6m.results.infer <- postProcessEvaluateFounders( "v3", "1m6m", use.infer = TRUE );
out.file <- paste( RESULTS.DIR, "v3/1m6m/infer/postProcessEvaluateFounders.tab", sep = "" );
write.table( apply( caprisa002.v3.1m6m.results.infer, 1:2, function( .x ) { sprintf( "%0.4f", .x ) } ), quote = FALSE, file = out.file, sep = "\t" );

### Identify results:
caprisa002.v3.1m.Identify.results <-
    evaluateResultsMatrix( caprisa002.v3.1m.results, "Identify (1m) founders", "Founder Mean HD to True\nCaprisa 002 V3\nIdentify (1m)" )
pdf( file = paste( RESULTS.DIR, "caprisa002.v3.1m.Identify.pdf", sep = "" ) );
caprisa002.v3.1m.Identify.results$ggp;
dev.off();
caprisa002.v3.1m.Identify.results$stats
#                V3.AA       V3.NA
# n         20.00000000 20.00000000
# Min.       0.00000000  0.00000000
# 1st Qu.    0.00000000  0.00000000
# Median     0.00000000  0.00000000
# Mean       0.00593000  0.00302600
# 3rd Qu.    0.00000000  0.00000000
# Max.       0.10630000  0.05246000
# sd         0.02375953  0.01169363
# frac.zero  0.85000000  0.80000000

caprisa002.v3.6m.Identify.results <-
    evaluateResultsMatrix( caprisa002.v3.6m.results, "Identify (6m) founders", "Founder Mean HD to True\nCaprisa 002 V3\nIdentify (6m)" )
pdf( file = paste( RESULTS.DIR, "caprisa002.v3.6m.Identify.pdf", sep = "" ) );
caprisa002.v3.6m.Identify.results$ggp;
dev.off();
caprisa002.v3.6m.Identify.results$stats
#                V3.AA        V3.NA
# n         18.00000000 18.000000000
# Min.       0.00000000  0.000000000
# 1st Qu.    0.00000000  0.000000000
# Median     0.00000000  0.000000000
# Mean       0.00697700  0.003287000
# 3rd Qu.    0.00000000  0.000973000
# Max.       0.07160000  0.036230000
# sd         0.01809126  0.008702398
# frac.zero  0.77777778  0.722222222

caprisa002.v3.1m6m.Identify.results <-
    evaluateResultsMatrix( caprisa002.v3.1m6m.results, "Identify (1m6m) founders", "Founder Mean HD to True\nCaprisa 002 V3\nIdentify (1m6m)" )
pdf( file = paste( RESULTS.DIR, "caprisa002.v3.1m6m.Identify.pdf", sep = "" ) );
caprisa002.v3.1m6m.Identify.results$ggp;
dev.off();
caprisa002.v3.1m6m.Identify.results$stats
#                V3.AA        V3.NA
# n         17.00000000 17.000000000
# Min.       0.00000000  0.000000000
# 1st Qu.    0.00000000  0.000000000
# Median     0.00000000  0.000000000
# Mean       0.00568400  0.002924000
# 3rd Qu.    0.00000000  0.000000000
# Max.       0.07160000  0.036230000
# sd         0.01780483  0.008830657
# frac.zero  0.82352941  0.764705882

### Infer results:
caprisa002.v3.1m.Infer.results <-
    evaluateResultsMatrix( caprisa002.v3.1m.results.infer, "Infer (1m) founders", "Founder Mean HD to True\nCaprisa 002 V3\nInfer (1m)" )
pdf( file = paste( RESULTS.DIR, "caprisa002.v3.1m.Infer.pdf", sep = "" ) );
caprisa002.v3.1m.Infer.results$ggp;
dev.off();
caprisa002.v3.1m.Infer.results$stats
#                V3.AA       V3.NA
# n         19.00000000 19.00000000
# Min.       0.00000000  0.00000000
# 1st Qu.    0.00000000  0.00000000
# Median     0.00000000  0.00000000
# Mean       0.00825500  0.00423100
# 3rd Qu.    0.00000000  0.00143300
# Max.       0.11970000  0.05530000
# sd         0.02778673  0.01272688
# frac.zero  0.78947368  0.68421053

caprisa002.v3.6m.Infer.results <-
    evaluateResultsMatrix( caprisa002.v3.6m.results.infer, "Infer (6m) founders", "Founder Mean HD to True\nCaprisa 002 V3\nInfer (6m)" )
pdf( file = paste( RESULTS.DIR, "caprisa002.v3.6m.Infer.pdf", sep = "" ) );
caprisa002.v3.6m.Infer.results$ggp;
dev.off();
caprisa002.v3.6m.Infer.results$stats
#                V3.AA        V3.NA
# n         18.00000000 18.000000000
# Min.       0.00000000  0.000000000
# 1st Qu.    0.00000000  0.000000000
# Median     0.00000000  0.000000000
# Mean       0.00660100  0.002742000
# 3rd Qu.    0.00000000  0.000311600
# Max.       0.08165000  0.033790000
# sd         0.01993176  0.008093776
# frac.zero  0.77777778  0.722222222

caprisa002.v3.1m6m.Infer.results <-
    evaluateResultsMatrix( caprisa002.v3.1m6m.results.infer, "Infer (1m6m) founders", "Founder Mean HD to True\nCaprisa 002 V3\nInfer (1m6m)" )
pdf( file = paste( RESULTS.DIR, "caprisa002.v3.1m6m.Infer.pdf", sep = "" ) );
caprisa002.v3.1m6m.Infer.results$ggp;
dev.off();
caprisa002.v3.1m6m.Infer.results$stats
#                V3.AA       V3.NA
# n         16.00000000 16.00000000
# Min.       0.00000000  0.00000000
# 1st Qu.    0.00000000  0.00000000
# Median     0.00000000  0.00000000
# Mean       0.00866300  0.00410300
# 3rd Qu.    0.00031130  0.00092430
# Max.       0.10140000  0.05008000
# sd         0.02577062  0.01250597
# frac.zero  0.75000000  0.68750000

### Identify, over time:
caprisa002.v3.6m.vs.1m.Identify.results <-
    compareResultsMatrix( caprisa002.v3.6m.results, caprisa002.v3.1m.results, "6m", "1m", "Caprisa 002 V3\nIdentify: 6m vs 1m" );
pdf( file = paste( RESULTS.DIR, "caprisa002.v3.6m.vs.1m.Identify.pdf", sep = "" ) );
caprisa002.v3.6m.vs.1m.Identify.results$ggp;
dev.off();
## Conclusion: Using Identify with Caprisa 002 v3, we do aobut the same with only 6m data as with only 1m data.

caprisa002.v3.1m6m.vs.1m.Identify.results <-
    compareResultsMatrix( caprisa002.v3.1m6m.results, caprisa002.v3.1m.results, "1m6m", "1m", "Caprisa 002 V3\nIdentify: 1m6m vs 1m" );
pdf( file = paste( RESULTS.DIR, "caprisa002.v3.1m6m.vs.1m.Identify.pdf", sep = "" ) );
caprisa002.v3.1m6m.vs.1m.Identify.results$ggp;
dev.off();
## Conclusion: Using Identify with Caprisa 002 v3, we do aobut the same with both 1m and 6m data as with only 1m data.

### Infer, over time:
caprisa002.v3.6m.vs.1m.Infer.results <-
    compareResultsMatrix( caprisa002.v3.6m.results.infer, caprisa002.v3.1m.results.infer, "6m", "1m", "Caprisa 002 V3\nInfer: 6m vs 1m" );
pdf( file = paste( RESULTS.DIR, "caprisa002.v3.6m.vs.1m.Infer.pdf", sep = "" ) );
caprisa002.v3.6m.vs.1m.Infer.results$ggp;
dev.off();
## Conclusion: Using Infer, we do about the same with both 6m data than with only 1m data.

caprisa002.v3.1m6m.vs.1m.Infer.results <-
    compareResultsMatrix( caprisa002.v3.1m6m.results.infer, caprisa002.v3.1m.results.infer, "1m6m", "1m", "Caprisa 002 V3\nInfer: 1m6m vs 1m" );
pdf( file = paste( RESULTS.DIR, "caprisa002.v3.1m6m.vs.1m.Infer.pdf", sep = "" ) );
caprisa002.v3.1m6m.vs.1m.Infer.results$ggp;
dev.off();
## Conclusion: Using Infer, we do about the same with both 1m6m data than with only 1m data.

### IdentifyVsInfer:
caprisa002.v3.1m.IdentifyVsInfer.results <-
    compareResultsMatrix( caprisa002.v3.1m.results, caprisa002.v3.1m.results.infer, "Identify", "Infer", "Caprisa 002 V3 1m\nIdentify vs Infer" );
pdf( file = paste( RESULTS.DIR, "caprisa002.v3.1m.IdentifyVsInfer.pdf", sep = "" ) );
caprisa002.v3.1m.IdentifyVsInfer.results$ggp;
dev.off();
## Conclusion: Using 1m Caprisa 002 data, we do about the same when using Infer as when using Identify.

caprisa002.v3.6m.IdentifyVsInfer.results <-
    compareResultsMatrix( caprisa002.v3.6m.results, caprisa002.v3.6m.results.infer, "Identify", "Infer", "Caprisa 002 V3 6m\nIdentify vs Infer" );
pdf( file = paste( RESULTS.DIR, "caprisa002.v3.6m.IdentifyVsInfer.pdf", sep = "" ) );
caprisa002.v3.6m.IdentifyVsInfer.results$ggp;
dev.off();
## Conclusion: Using 6m Caprisa 002 data, we do about the same when using Infer as when using Identify (maybe slightly better using Identify in NAs)

caprisa002.v3.1m6m.IdentifyVsInfer.results <-
    compareResultsMatrix( caprisa002.v3.1m6m.results, caprisa002.v3.1m6m.results.infer, "Identify", "Infer", "Caprisa 002 V3 1m6m\nIdentify vs Infer" );
pdf( file = paste( RESULTS.DIR, "caprisa002.v3.1m6m.IdentifyVsInfer.pdf", sep = "" ) );
caprisa002.v3.1m6m.IdentifyVsInfer.results$ggp;
dev.off();
## Conclusion: Using 1m6m Caprisa 002 data, we do about the same when using Infer as when using Identify.


## RV217
rv217.nflg.1m.results <- postProcessEvaluateFounders( "nflg", "1m" );
out.file <- paste( RESULTS.DIR, "nflg/1m/postProcessEvaluateFounders.tab", sep = "" );
write.table( apply( rv217.nflg.1m.results, 1:2, function( .x ) { sprintf( "%0.4f", .x ) } ), quote = FALSE, file = out.file, sep = "\t" );
rv217.nflg.1m.results.infer <- postProcessEvaluateFounders( "nflg", "1m", use.infer = TRUE );
out.file <- paste( RESULTS.DIR, "nflg/1m/infer/postProcessEvaluateFounders.tab", sep = "" );
write.table( apply( rv217.nflg.1m.results.infer, 1:2, function( .x ) { sprintf( "%0.4f", .x ) } ), quote = FALSE, file = out.file, sep = "\t" );

rv217.nflg.6m.results <- postProcessEvaluateFounders( "nflg", "6m" );
out.file <- paste( RESULTS.DIR, "nflg/6m/postProcessEvaluateFounders.tab", sep = "" );
write.table( apply( rv217.nflg.6m.results, 1:2, function( .x ) { sprintf( "%0.4f", .x ) } ), quote = FALSE, file = out.file, sep = "\t" );
rv217.nflg.6m.results.infer <- postProcessEvaluateFounders( "nflg", "6m", use.infer = TRUE );
out.file <- paste( RESULTS.DIR, "nflg/6m/infer/postProcessEvaluateFounders.tab", sep = "" );
write.table( apply( rv217.nflg.6m.results.infer, 1:2, function( .x ) { sprintf( "%0.4f", .x ) } ), quote = FALSE, file = out.file, sep = "\t" );

rv217.nflg.1m6m.results <- postProcessEvaluateFounders( "nflg", "1m6m" );
out.file <- paste( RESULTS.DIR, "nflg/1m6m/postProcessEvaluateFounders.tab", sep = "" );
write.table( apply( rv217.nflg.1m6m.results, 1:2, function( .x ) { sprintf( "%0.4f", .x ) } ), quote = FALSE, file = out.file, sep = "\t" );
rv217.nflg.1m6m.results.infer <- postProcessEvaluateFounders( "nflg", "1m6m", use.infer = TRUE );
out.file <- paste( RESULTS.DIR, "nflg/1m6m/infer/postProcessEvaluateFounders.tab", sep = "" );
write.table( apply( rv217.nflg.1m6m.results.infer, 1:2, function( .x ) { sprintf( "%0.4f", .x ) } ), quote = FALSE, file = out.file, sep = "\t" );

### Identify results:
rv217.nflg.1m.Identify.results <-
    evaluateResultsMatrix( rv217.nflg.1m.results, "Identify (1m) founders", "Founder Mean HD to True\nRV217 NFLG\nIdentify (1m)" )
pdf( file = paste( RESULTS.DIR, "rv217.nflg.1m.Identify.pdf", sep = "" ) );
rv217.nflg.1m.Identify.results$ggp;
dev.off();
rv217.nflg.1m.Identify.results$stats
#               NFLG.AA      NFLG.NA
# n         36.000000000 36.000000000
# Min.       0.000000000  0.000000000
# 1st Qu.    0.000274700  0.000000000
# Median     0.000805500  0.000192500
# Mean       0.002535000  0.001467000
# 3rd Qu.    0.001570000  0.000659500
# Max.       0.025360000  0.013160000
# sd         0.004928037  0.003019404
# frac.zero  0.166666667  0.388888889

rv217.nflg.6m.Identify.results <-
    evaluateResultsMatrix( rv217.nflg.6m.results, "Identify (6m) founders", "Founder Mean HD to True\nRV217 NFLG\nIdentify (6m)" )
pdf( file = paste( RESULTS.DIR, "rv217.nflg.6m.Identify.pdf", sep = "" ) );
rv217.nflg.6m.Identify.results$ggp;
dev.off();
rv217.nflg.6m.Identify.results$stats
#               NFLG.AA      NFLG.NA
# n         34.000000000 34.000000000
# Min.       0.000000000  0.000000000
# 1st Qu.    0.001559000  0.000410300
# Median     0.003258000  0.001387000
# Mean       0.005606000  0.003064000
# 3rd Qu.    0.007899000  0.004123000
# Max.       0.021760000  0.015270000
# sd         0.005888187  0.003949968
# frac.zero  0.029411765  0.117647059

rv217.nflg.1m6m.Identify.results <-
    evaluateResultsMatrix( rv217.nflg.1m6m.results, "Identify (1m6m) founders", "Founder Mean HD to True\nRV217 NFLG\nIdentify (1m6m)" )
pdf( file = paste( RESULTS.DIR, "rv217.nflg.1m6m.Identify.pdf", sep = "" ) );
rv217.nflg.1m6m.Identify.results$ggp;
dev.off();
rv217.nflg.1m6m.Identify.results$stats
#               NFLG.AA     NFLG.NA
# n         24.000000000 24.00000000
# Min.       0.000000000  0.00000000
# 1st Qu.    0.000554600  0.00010860
# Median     0.000949400  0.00027530
# Mean       0.003868000  0.00126900
# 3rd Qu.    0.002519000  0.00138400
# Max.       0.046160000  0.01048000
# sd         0.009502685  0.00239576
# frac.zero  0.083333333  0.25000000

### Infer results:
rv217.nflg.1m.Infer.results <-
    evaluateResultsMatrix( rv217.nflg.1m.results.infer, "Infer (1m) founders", "Founder Mean HD to True\nRV217 NFLG\nInfer (1m)" )
pdf( file = paste( RESULTS.DIR, "rv217.nflg.1m.Infer.pdf", sep = "" ) );
rv217.nflg.1m.Infer.results$ggp;
dev.off();
rv217.nflg.1m.Infer.results$stats
#               NFLG.AA      NFLG.NA
# n         34.000000000 34.000000000
# Min.       0.000000000  0.000000000
# 1st Qu.    0.000000000  0.000000000
# Median     0.000532600  0.000039210
# Mean       0.003338000  0.002039000
# 3rd Qu.    0.001961000  0.000811700
# Max.       0.028870000  0.018290000
# sd         0.006942697  0.004343099
# frac.zero  0.382352941  0.500000000

rv217.nflg.6m.Infer.results <-
    evaluateResultsMatrix( rv217.nflg.6m.results.infer, "Infer (6m) founders", "Founder Mean HD to True\nRV217 NFLG\nInfer (6m)" )
pdf( file = paste( RESULTS.DIR, "rv217.nflg.6m.Infer.pdf", sep = "" ) );
rv217.nflg.6m.Infer.results$ggp;
dev.off();
rv217.nflg.6m.Infer.results$stats
#               NFLG.AA      NFLG.NA
# n         33.000000000 33.000000000
# Min.       0.000000000  0.000000000
# 1st Qu.    0.002145000  0.000806300
# Median     0.004231000  0.001411000
# Mean       0.007215000  0.003770000
# 3rd Qu.    0.008830000  0.003671000
# Max.       0.032580000  0.019270000
# sd         0.007983511  0.004920596
# frac.zero  0.121212121  0.121212121

rv217.nflg.1m6m.Infer.results <-
    evaluateResultsMatrix( rv217.nflg.1m6m.results.infer, "Infer (1m6m) founders", "Founder Mean HD to True\nRV217 NFLG\nInfer (1m6m)" )
pdf( file = paste( RESULTS.DIR, "rv217.nflg.1m6m.Infer.pdf", sep = "" ) );
rv217.nflg.1m6m.Infer.results$ggp;
dev.off();
rv217.nflg.1m6m.Infer.results$stats
#               NFLG.AA      NFLG.NA
# n         22.000000000 22.000000000
# Min.       0.000000000  0.000000000
# 1st Qu.    0.001422000  0.000501400
# Median     0.002569000  0.001087000
# Mean       0.004328000  0.002429000
# 3rd Qu.    0.004020000  0.001538000
# Max.       0.017440000  0.012530000
# sd         0.004894294  0.003625868
# frac.zero  0.045454545  0.090909091

### Identify, over time:
rv217.nflg.6m.vs.1m.Identify.results <-
    compareResultsMatrix( rv217.nflg.6m.results, rv217.nflg.1m.results, "6m", "1m", "RV217 NFLG\nIdentify: 6m vs 1m" );
pdf( file = paste( RESULTS.DIR, "rv217.nflg.6m.vs.1m.Identify.pdf", sep = "" ) );
rv217.nflg.6m.vs.1m.Identify.results$ggp;
dev.off();
## Conclusion: Using Identify, we do worse with only 6m data than with only 1m data.

rv217.nflg.1m6m.vs.1m.Identify.results <-
    compareResultsMatrix( rv217.nflg.1m6m.results, rv217.nflg.1m.results, "1m6m", "1m", "RV217 NFLG\nIdentify: 1m6m vs 1m" );
pdf( file = paste( RESULTS.DIR, "rv217.nflg.1m6m.vs.1m.Identify.pdf", sep = "" ) );
rv217.nflg.1m6m.vs.1m.Identify.results$ggp;
dev.off();
## Conclusion: Using Identify, we do about the same with both 1m6m data than with only 1m data.  Maybe in one case we do better with 1m6m on the AAs, and there is some trend towards doing better overall on the NAs -- if anything it is BETTER using 1m rather than 1m6m!

### Infer, over time:
rv217.nflg.6m.vs.1m.Infer.results <-
    compareResultsMatrix( rv217.nflg.6m.results.infer, rv217.nflg.1m.results.infer, "6m", "1m", "RV217 NFLG\nInfer: 6m vs 1m" );
pdf( file = paste( RESULTS.DIR, "rv217.nflg.6m.vs.1m.Infer.pdf", sep = "" ) );
rv217.nflg.6m.vs.1m.Infer.results$ggp;
dev.off();
## Conclusion: Using Infer, we do worse with only 6m data than with only 1m data.

rv217.nflg.1m6m.vs.1m.Infer.results <-
    compareResultsMatrix( rv217.nflg.1m6m.results.infer, rv217.nflg.1m.results.infer, "1m6m", "1m", "RV217 NFLG\nInfer: 1m6m vs 1m" );
pdf( file = paste( RESULTS.DIR, "rv217.nflg.1m6m.vs.1m.Infer.pdf", sep = "" ) );
rv217.nflg.1m6m.vs.1m.Infer.results$ggp;
dev.off();
## Conclusion: Suprisingly, using Infer, we do worse with both 1m6m data than with only 1m data.

### IdentifyVsInfer:
rv217.nflg.1m.IdentifyVsInfer.results <-
    compareResultsMatrix( rv217.nflg.1m.results, rv217.nflg.1m.results.infer, "Identify", "Infer", "RV217 NFLG 1m\nIdentify vs Infer" );
pdf( file = paste( RESULTS.DIR, "rv217.nflg.1m.IdentifyVsInfer.pdf", sep = "" ) );
rv217.nflg.1m.IdentifyVsInfer.results$ggp;
dev.off();
## Conclusion: Using 1m RV217 NFLG data, we do about the same when using Infer as when using Identify, but if anything the Infer results are better.

rv217.nflg.6m.IdentifyVsInfer.results <-
    compareResultsMatrix( rv217.nflg.6m.results, rv217.nflg.6m.results.infer, "Identify", "Infer", "RV217 NFLG 6m\nIdentify vs Infer" );
pdf( file = paste( RESULTS.DIR, "rv217.nflg.6m.IdentifyVsInfer.pdf", sep = "" ) );
rv217.nflg.6m.IdentifyVsInfer.results$ggp;
dev.off();
## Conclusion: Using 6m RV217 NFLG data, we do significantly worse when using Identify than when using Infer.

rv217.nflg.1m6m.IdentifyVsInfer.results <-
    compareResultsMatrix( rv217.nflg.1m6m.results, rv217.nflg.1m6m.results.infer, "Identify", "Infer", "RV217 NFLG 1m6m\nIdentify vs Infer" );
pdf( file = paste( RESULTS.DIR, "rv217.nflg.1m6m.IdentifyVsInfer.pdf", sep = "" ) );
rv217.nflg.1m6m.IdentifyVsInfer.results$ggp;
dev.off();
## Conclusion: Using 1m6m RV217 NFLG data, we do about the same when using Infer as when using Identify, with one case of Identify doing better on AAs -- however in NAs the results are significantly better using Infer.


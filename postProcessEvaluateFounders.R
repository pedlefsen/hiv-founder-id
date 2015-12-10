### NOTE: This is assuming you have _already_ compiled the evaluateFounders.tbl results into a single file at the top-level dir (per time point, region): eg /fh/fast/edlefsen_p/bakeoff_analysis_results/raw/nflg/1m/evaluateFounders.tbl
##  NOTE continued: for now I'm manually creating these files because I am lazy and in one motion within emacs I can get rid of the extra table headers.
# eg cat /fh/fast/edlefsen_p/bakeoff_analysis_results/raw/nflg/1m/*/evaluateFounders.tbl > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw/nflg/1m/evaluateFounders.tbl
# then open and remove the extra instances of the first (header) line.

## TODO: Rerun the evaluateFoundersFromInfer.bash script (and Cap variant). [RUNNING]
## TODO: Gather the results as described at the top of this file.

# From this file we read in indicators of whether to use the multiple- or single-founder true profile.
rv217.gold.standards.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/RV217_gold_standards.csv" );
rv217.gold.is.multiple <- rv217.gold.standards.in[ , "gold.is.multiple" ];
names( rv217.gold.is.multiple ) <- rv217.gold.standards.in[ , "ptid" ];

# From this file we read in indicators of whether to use the multiple- or single-founder true profile.
caprisa002.gold.standards.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/caprisa_002_gold_standards.csv" );
caprisa002.gold.is.multiple <- caprisa002.gold.standards.in[ , "gold.is.multiple" ];
names( caprisa002.gold.is.multiple ) <- caprisa002.gold.standards.in[ , "ptid" ];

## Read in the rv217 nflg evaluateFounders results for the identify-founders method.
rv217.nflg.1m.results.in <- read.delim( "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw/nflg/1m/evaluateFounders.tbl", sep = "\t" );
rv217.nflg.1m.results.in.estimates.ptid <-
    gsub( ".*rv217_([^_]+)_.+", "\\1", as.character( rv217.nflg.1m.results.in[ , "estimates.file" ] ) );
rv217.nflg.1m.results.in.truths.ptid <-
    gsub( ".*rv217_([^_]+)_.+", "\\1", as.character( rv217.nflg.1m.results.in[ , "truths.file" ] ) );
stopifnot( all( rv217.nflg.1m.results.in.estimates.ptid == rv217.nflg.1m.results.in.truths.ptid ) );
rv217.nflg.1m.results.in.ptids <-
    unique( rv217.nflg.1m.results.in.truths.ptid );

## Note that we do not subset to same-region results, since the regions overlap, we can aggregate HDs (ignoring gaps) by summing numerator and denominators.  This solidifies my thinking that the one that includes gaps does not make sense here, and gaps might be inserted by geneCutter to make codons work as reasonable AAs in translation, which is really no fault of the estimation, as it might be context-dependent (eg if a substitution in one homolgous pair results in that pair being not-aligned after codon-alignment because the placement leads to different best-guesses of the functional codon [though my thinking is that this should be exceedingly rare]).

rv217.nflg.1m.results.fromto <- rv217.nflg.1m.results.in[ , 1:2, drop = FALSE ];
###### ERE I AM UPDATING FIXING ETC.
# These should be removed methought, now not so sure: rv217.nflg.1m.results[ ( rv217.nflg.1m.truths.file.region == "LH" & rv217.nflg.1m.estimates.file.region == "RH" ) | ( rv217.nflg.1m.truths.file.region == "RH" & rv217.nflg.1m.estimates.file.region == "LH" ), ]
## THINKING KEEP.  JUSTIFICATION IS THAT IF THEY ALIGN THEY COUNT. WHY NOT?
### BUT DO WE NEED TO WEIGH THEM WHEN AVERAGING OVER REGIONS?  by number of seqs, I guess, to mimic what we would have gotten if we made one alignment and computed the HD.  yes?  yes.  so we need the number of seqs yarg.  that means an additional output or two in evaluateFounders.R.. Ok.
## NO we do not need to weigh them.  Oh wait, hmm.  ... ..  we need the denominators.

stopifnot( length( rv217.nflg.1m.truths.file.region ) == nrow( rv217.nflg.1m.results ) );
## Here we do the sanity check to ensure that every ( from, to ) pair appears at most once.
.fromto <- apply( rv217.nflg.1m.results.fromto, 1, paste, collapse = "---" );
stopifnot( length( .fromto ) == length( unique( .fromto ) ) );

rv217.nflg.1m.estimates.file.ptid <-
    gsub( ".*rv217_([^_]+)_.+", "\\1", as.character( rv217.nflg.1m.results.fromto[ , "estimates.file" ] ) );
rv217.nflg.1m.truths.file.ptid <-
    gsub( ".*rv217_([^_]+)_.+", "\\1", as.character( rv217.nflg.1m.results.fromto[ , "truths.file" ] ) );
stopifnot( all( rv217.nflg.1m.truths.file.ptid == rv217.nflg.1m.estimates.file.ptid ) );
stopifnot( all( rv217.nflg.1m.truths.file.ptid %in% names( rv217.gold.is.multiple ) ) );

rv217.nflg.1m.estimates.file.is.multiple <-
    sapply( as.character( rv217.nflg.1m.results.fromto[ , "estimates.file" ] ), function( .str ) { length( grep( ".*_single.+", .str, invert = TRUE ) ) > 0 } );
rv217.nflg.1m.truths.file.is.multiple <-
    sapply( as.character( rv217.nflg.1m.results.fromto[ , "truths.file" ] ), function( .str ) { length( grep( ".*_single.+", .str, invert = TRUE ) ) > 0 } );

## Ok here is where we subset to use just the "single" or "multiple" results, depending on the value of gold.is.multiple.  Note that we prefer to use the matching result, but if it is not there our next best result is the one with the matching "truth" is.multiple.  That is, we can compare a multi-founder estimate to a single-founder truth, if there is only a multi-founder estimate.
rv217.nflg.1m.truths.matches.gold.is.multiple <- ( rv217.nflg.1m.truths.file.is.multiple == rv217.gold.is.multiple[ rv217.nflg.1m.truths.file.ptid ] );
rv217.nflg.1m.estimates.matches.gold.is.multiple <- ( rv217.nflg.1m.estimates.file.is.multiple == rv217.gold.is.multiple[ rv217.nflg.1m.truths.file.ptid ] );

rv217.nflg.1m.results.matches.gold.is.multiple <-
    rv217.nflg.1m.results.in[ ( rv217.nflg.1m.truths.matches.gold.is.multiple & rv217.nflg.1m.estimates.matches.gold.is.multiple ), , drop = FALSE ];
rv217.nflg.1m.results.matches.gold.is.multiple.truths.ptid <-
    rv217.nflg.1m.truths.file.ptid[ ( rv217.nflg.1m.truths.matches.gold.is.multiple & rv217.nflg.1m.estimates.matches.gold.is.multiple ) ];

# Stop if we lost any ptids.
stopifnot( length( unique( rv217.nflg.1m.estimates.file.ptid[ ( rv217.nflg.1m.truths.matches.gold.is.multiple & rv217.nflg.1m.estimates.matches.gold.is.multiple ) ] ) ) == length( rv217.nflg.1m.results.in.ptids ) );

# ## TODO: Read this in, use all the "fits" and "is.starlike" etc options to choose whether to use single or multiple in the estimated set.  NOTE: FOR NOW we just use the "right" one for the job, ie the same one we use for the "truth".
# rv217.nflg.1m.identifyfounders.in <- read.delim( "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw/nflg/1m/identify_founders.tab", sep = "\t" );
# rv217.nflg.1m.identifyfounders <- as.matrix( rv217.nflg.1m.identifyfounders.in );
# suppressWarnings( mode( rv217.nflg.1m.identifyfounders ) <- "numeric" );

## Also (in addition to the includingGaps results; see above), we remove the interesting (but not for reporting) perspectives results, which are components of the result that we use: truths.perspective is the average over true founders of the HD of the smallest-HD estimate for that true founder (even if it is the same as another one).
rv217.nflg.1m.results.prep <- as.matrix( rv217.nflg.1m.results.matches.gold.is.multiple[ , 3:ncol( rv217.nflg.1m.results.matches.gold.is.multiple ) ] );
suppressWarnings( mode( rv217.nflg.1m.results.prep ) <- "numeric" );

rv217.nflg.1m.cols.to.keep <-
    grep( "FULL_SEQUENCE", grep( "\\.sum\\.ignoringGaps", colnames( rv217.nflg.1m.results.prep ), value = TRUE ), invert = TRUE, value = TRUE );
rv217.nflg.1m.cols.to.keep.prefixes <- gsub( "^(.+)\\.perspective.+", "\\1", rv217.nflg.1m.cols.to.keep );
## Ok so now we need to recompute per-ptid averages; first we sum like columns for a ptid.
rv217.nflg.1m.results.prep2 <- t(
    sapply( rv217.nflg.1m.results.in.ptids, function( .ptid ) {
        .ptid.sums <- 
            apply( rv217.nflg.1m.results.prep[ rv217.nflg.1m.results.matches.gold.is.multiple.truths.ptid == .ptid, rv217.nflg.1m.cols.to.keep, drop = FALSE ], 2, sum, na.rm = T );
        .rv <- sapply( unique( rv217.nflg.1m.cols.to.keep.prefixes ), function( .prefix ) {
            .ptid.sums[ paste( .prefix, "perspective.HD.sum.ignoringGaps", sep = "." ) ] /
            .ptid.sums[ paste( .prefix, "perspective.denominator.sum.ignoringGaps", sep = "." ) ]
        } );
        names( .rv ) <- unique( rv217.nflg.1m.cols.to.keep.prefixes );
        return( .rv );
    } ) );
## And finally replace it all with averages (we average the truths-perspective and the estimates-perspective).
rv217.nflg.1m.results.prep2.prefixes <- gsub( "^(.+)\\.(truths|estimates)", "\\1", colnames( rv217.nflg.1m.results.prep2 ) );

rv217.nflg.1m.results <- sapply( unique( rv217.nflg.1m.results.prep2.prefixes ), function( .prefix ) {
    apply( rv217.nflg.1m.results.prep2[ , paste( .prefix, c( "truths", "estimates" ), sep = "." ), drop = FALSE ], 1, mean, na.rm = T )
} );

#### ERE I AM. THAT SEEMS TO HAVE WORKED OUT OK, BUT MAYBE IF an is.multiple match isn't available, have a fallback option?  Anyway, this works enough for now.

rv217.nflg.6m.results.in <- read.delim( "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw/nflg/6m/evaluateFounders.tbl", sep = "\t" );
rv217.nflg.6m.results <- as.matrix( rv217.nflg.6m.results.in );
suppressWarnings( mode( rv217.nflg.6m.results ) <- "numeric" );

rv217.nflg.1m6m.results.in <- read.delim( "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw/nflg/1m6m/evaluateFounders.tbl", sep = "\t" );
rv217.nflg.1m6m.results <- as.matrix( rv217.nflg.1m6m.results.in );
suppressWarnings( mode( rv217.nflg.1m6m.results ) <- "numeric" );

## Read in the rv217 nflg evaluateFounders results for the identify-founders method.
caprisa002.nflg.1m.results.in <- read.delim( "/fh/fast/edlefsen_p/bakeoff_analysis_results/raw/v3/1m/evaluateFounders.tbl", sep = "\t" );
caprisa002.nflg.1m.results <- as.matrix( caprisa002.nflg.1m.results.in );
suppressWarnings( mode( caprisa002.nflg.1m.results ) <- "numeric" );

rv217.nflg.1m.results.mean.by.stat <- 
  apply( rv217.nflg.1m.results, 2, function( .column ) {
      .column[ !is.finite( .column ) ] <- NA;
      return( mean( .column, na.rm = T ) );
  } );
rv217.nflg.6m.results.mean.by.stat <- 
  apply( rv217.nflg.6m.results, 2, function( .column ) {
      .column[ !is.finite( .column ) ] <- NA;
      return( mean( .column, na.rm = T ) );
  } );
rv217.nflg.1m6m.results.mean.by.stat <- 
  apply( rv217.nflg.1m6m.results, 2, function( .column ) {
      .column[ !is.finite( .column ) ] <- NA;
      return( mean( .column, na.rm = T ) );
  } );


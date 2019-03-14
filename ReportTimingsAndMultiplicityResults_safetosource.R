library( "xtable" )
library( "gplots" ) # for heatmap.2
library( "ggplot2" )
library( "forestplot" )
library( "reshape2" )

.owd <- getwd();
setwd( "~/src/from-git/hiv-founder-id/" );
source( "getFilteredResultsTables_safetosource.R" );
source( "evaluateTimings_safetosource.R" );
source( "evaluateIsMultiple_safetosource.R" )
source( "daysFromLambda_safetosource.R" )
setwd( .owd );

## CONFIGURATION and DEFAULTS will usually be set by the caller.
## See ReportTimingsAndMultiplicityResults.Rnw where these are set. Eg:
#############
# CONFIGURATION 
############
# timings.include.methods <- c( "glm" ); # c( "glm", "step", "lasso" )
# timings.include.training.codes <- c( "RT", "Rt", "rT", "rt" ); 
# 
# ismultiple.show.heatmap = FALSE;
# if( ismultiple.show.heatmap ) {
#     ismultiple.include.training.codes <- c( "RT", "Rt", "rT", "rt" );
#     ismultiple.only.is.lasso <- FALSE;
#     ismultiple.exclude.continuous.predictors <- FALSE;
# } else {
#     ismultiple.include.training.codes <-  "RT";
#     ismultiple.only.is.lasso <- TRUE;
#     ismultiple.exclude.continuous.predictors <- TRUE;
# }
# 
# ###########
# ## DEFAULTS
# ###########
# 
# ## Configure which results we are displaying.
# BakeOff.RESULTS.DIR <- Sys.getenv( "BakeOff_RESULTS_DIR" );
# if( BakeOff.RESULTS.DIR != "" ) {
#     RESULTS.DIR <- BakeOff.RESULTS.DIR;
# } else {
#     #RESULTS.DIR <- "/fast/bakeoff_merged_analysis_sequences_unfiltered/results/";
#     RESULTS.DIR <- "/fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/";
#     #RESULTS.DIR <- "/fast/bakeoff_merged_analysis_sequences_filtered2019/results/";
# }
# RESULTS.DIRNAME <- "raw_fixed";
# 
# ## NOTE even setting this to TRUE will not force all to be recomputed. See below where force.recomputation = FALSE.
# FORCE.RECOMPUTATION <- FALSE;#TRUE;



#############
# FUNCTIONS
############

# For RV217 we are just interested in mean depth over NFLG-equivs. For now this means that we count half-genomes as 0.5 and we don't count env seqs.
remove.env.from.nflg.subtable <- function( .subtable, .ptid ) {
    # Special: for this table we are excluding the 21 env seqs available for ptid "20263" or the 20 env seqs for ptid "40100", since they aren't really NFLG-equivs and it's too complicated to explain this sublety in the table.
    if( length( .which.are.env <- grep( "_env\\.fasta$", .subtable[ , "file" ] ) ) > 0 ) {
        stopifnot( ( .ptid == "20263" ) || ( .ptid == "40100" ) );
        if( .ptid == "20263" ) {
            warning( "We are intentionally not counting the 21 env seqs from RV217 ptid \"20263\" in the table of input data" );
        } else {
            warning( "We are intentionally not counting the 20 env seqs from RV217 ptid \"40100\" in the table of input data" );
        }
        .subtable <- .subtable[ setdiff( 1:nrow( .subtable ), .which.are.env ), , drop = FALSE ];
    }
    return( .subtable );
} # remove.env.from.nflg.subtable (..)
get.nflg.depth.from.ptid.subtable <- function( .subtable ) {
    # Ok, RH and LH seqs each count for a half of an nflg-equiv. Note that this means oddly that if there's only RH (as there are for two of the ptids: 10066 and 20507), the mean depth is half the number of sequences -- true but maybe unclear.
    .which.are.RH <- grep( "_RH.*\\.fasta$", .subtable[ , "file" ] );
    .which.are.LH <- grep( "_LH.*\\.fasta$", .subtable[ , "file" ] );
    .which.are.NFLG <- grep( "_NFLG.*\\.fasta$", .subtable[ , "file" ] );
    stopifnot( sort( c( .which.are.RH, .which.are.LH, .which.are.NFLG ) ) == 1:nrow( .subtable ) );
    stopifnot( length( .which.are.RH ) <= 1 );
    stopifnot( length( .which.are.LH ) <= 1 );
    stopifnot( length( .which.are.NFLG ) <= 1 );
    .nflg.depth.for.ptid <- 0;
    if( length( .which.are.RH ) == 1 ) {
        .nflg.depth.for.ptid <- .nflg.depth.for.ptid + 0.5*.subtable[ .which.are.RH, "num.seqs" ];
    }
    if( length( .which.are.LH ) == 1 ) {
        .nflg.depth.for.ptid <- .nflg.depth.for.ptid + 0.5*.subtable[ .which.are.LH, "num.seqs" ];
    }
    if( length( .which.are.NFLG ) == 1 ) {
        .nflg.depth.for.ptid <- .nflg.depth.for.ptid + 1.0*.subtable[ .which.are.NFLG, "num.seqs" ];
    }
    return( .nflg.depth.for.ptid );
} # get.nflg.depth.from.ptid.subtable (..)
# For RV217 we are interested in three widths, one for nflg, one for lh, one for rh.
get.nflg.widths.from.ptid.subtable <- function( .subtable ) {
    # NOTE: there's only RH for two of the ptids at 1m: 10066 and 20507.
    .which.are.RH <- grep( "_RH.*\\.fasta$", .subtable[ , "file" ] );
    .which.are.LH <- grep( "_LH.*\\.fasta$", .subtable[ , "file" ] );
    .which.are.NFLG <- grep( "_NFLG.*\\.fasta$", .subtable[ , "file" ] );
    stopifnot( sort( c( .which.are.RH, .which.are.LH, .which.are.NFLG ) ) == 1:nrow( .subtable ) );
    stopifnot( length( .which.are.RH ) <= 1 );
    stopifnot( length( .which.are.LH ) <= 1 );
    stopifnot( length( .which.are.NFLG ) <= 1 );
    .nflg <- NA;
    .lh <- NA;
    .rh <- NA;
    if( length( .which.are.NFLG ) == 1 ) {
        .nflg <- .subtable[ .which.are.NFLG, "PFitter.nbases" ];
    }
    if( length( .which.are.RH ) == 1 ) {
        .rh <- .subtable[ .which.are.RH, "PFitter.nbases" ];
    }
    if( length( .which.are.LH ) == 1 ) {
        .lh <- .subtable[ .which.are.LH, "PFitter.nbases" ];
    }
    return( c(
        "nflg" = .nflg,
        "lh" = .lh,
        "rh" = .rh
        ) );
} # get.nflg.widths.from.ptid.subtable (..)

    missing.row.safe.cbind <- function ( matA, matB ) {
        all.rownames <- union( rownames( matA ), rownames( matB ) );
        .rv <- matrix( NA, ncol = ncol( matA ) + ncol( matB ), nrow = length( all.rownames ) );
        rownames( .rv ) <- all.rownames;
        colnames( .rv ) <-
            c( colnames( matA ), colnames( matB ) );
        .rv[ rownames( matA ), 1:ncol( matA ) ] <-
            as.matrix( matA );
        .rv[ rownames( matB ), ncol( matA ) + 1:ncol( matB ) ] <-
            as.matrix( matB );
        
        return( .rv );
    } # missing.row.safe.cbind (..)


clean.timings.mat <- function ( timings.mat ) {
    # merge bounded Infer rows
    rownames( timings.mat ) <- gsub( "Infer sampledwidth uniform .*$", "Infer (w/ prior bounds)", rownames( timings.mat ) );
    #rownames( timings.mat ) <- gsub( "Infer uniform .*$", "Infer (w/ prior bounds)", rownames( timings.mat ) );
    if( length( grep( "Infer \\(w/ prior bounds\\)", rownames( timings.mat ) ) ) > 1 ) {
        .rows <- grep( "Infer \\(w/ prior bounds\\)", rownames( timings.mat ) );
        .submat <- timings.mat[ .rows, , drop = FALSE ];
        .new.row <- apply( .submat, 2, function( .col ) { .non.na <- .col[ !is.na( .col ) ]; if( length( .non.na ) == 0 ) { return( NA ); }; stopifnot( length( .non.na ) == 1 ); return( .non.na ); } );
        # Replace the first of the rows with the merged one.
        timings.mat[ .rows[ 1 ], ] <- .new.row;
        # Get rid of the rest.
        timings.mat <- timings.mat[ -.rows[ -1 ], , drop = FALSE ];
    }
    # remove unused bounded Infer rows
    #timings.mat <- timings.mat[ grep( "Infer uniform .*$", rownames( timings.mat ), invert = TRUE ), , drop = FALSE ];
    
    # merge bounded COB rows, but first remove the non-sampledwidth ones.
    # remove unused bounded COB rows
    timings.mat <- timings.mat[ grep( "COB uniform .*$", rownames( timings.mat ), invert = TRUE ), , drop = FALSE ];
    timings.mat <- timings.mat[ grep( "COB exponentialwidth .*$", rownames( timings.mat ), invert = TRUE ), , drop = FALSE ];
    rownames( timings.mat ) <- gsub( "COB sampledwidth uniform .*$", "COB (w/ prior bounds)", rownames( timings.mat ) );
    if( length( grep( "COB \\(w/ prior bounds\\)", rownames( timings.mat ) ) ) > 1 ) {
        .rows <- grep( "COB \\(w/ prior bounds\\)", rownames( timings.mat ) );
        .submat <- timings.mat[ .rows, , drop = FALSE ];
        .new.row <- apply( .submat, 2, function( .col ) { .non.na <- .col[ !is.na( .col ) ]; if( length( .non.na ) == 0 ) { return( NA ); }; stopifnot( length( .non.na ) == 1 ); return( .non.na ); } );
        # Replace the first of the rows with the merged one.
        timings.mat[ .rows[ 1 ], ] <- .new.row;
        # Get rid of the rest.
        timings.mat <- timings.mat[ -.rows[ -1 ], , drop = FALSE ];
    }
    
    # Rename rows
    rownames( timings.mat ) <- gsub( "^COB .*$", "center of prior bounds", rownames( timings.mat ) );
    rownames( timings.mat ) <- gsub( "Infer", "PREAST", rownames( timings.mat ) );
    # Reorder the rows, so "none" is first, then "PREAST" alone, then "PREAST (w/ prior bounds)"
    .new.row.order <- c( "none", "center of prior bounds", grep( "InSites", rownames( timings.mat ), value = TRUE ), grep( "PREAST", rownames( timings.mat ), value = TRUE ), grep( "PFitter", rownames( timings.mat ), value = TRUE ) );
    .new.row.order <- .new.row.order[ .new.row.order %in% rownames( timings.mat ) ];
    stopifnot( length( setdiff( rownames( timings.mat ), .new.row.order ) ) == 0 );
    timings.mat <- timings.mat[ .new.row.order, , drop = FALSE ];
    #timings.mat <- timings.mat[ , 1:4 ];
    
    return( timings.mat );
} # clean.timings.mat (..)

prepare.timings.mat <- function ( timings.mat, include.intercept, include.methods = timings.include.methods, include.training.codes = timings.include.training.codes ) {
    # Also exclude nointercept no-bounds columns.
    #timings.mat <- timings.mat[ , grep( "(?:glm nointercept validation results|lasso)", colnames( timings.mat ), invert = TRUE ), drop = FALSE ];
    
    rows.that.are.1m <- grep( " mtn003", rownames( timings.mat ) );
    rows.that.are.6m <- grep( " hvtn502", rownames( timings.mat ) );
    #rows.that.are.1m.6m <- grep( " 1mmtn003 6mhvtn502", rownames( timings.mat ) );
    #cols.that.are.1m <- grep( "\\.1m\\.[^6]", colnames( timings.mat ) );
    cols.that.are.1m.6m <- grep( "\\.1m\\.6m", colnames( timings.mat ) );
    #cols.that.are.6m <- setdiff( cols.that.are.1m.6m, grep( "\\.6m", colnames( timings.mat ) ) );
    
    # Ensure that only time-matching cells are non-na
    timings.mat[ rows.that.are.1m, cols.that.are.1m.6m ] <- NA;
    timings.mat[ rows.that.are.6m, cols.that.are.1m.6m ] <- NA;

    is.glm <- grep( "glm", colnames( timings.mat ) );
    is.step <- grep( "step", colnames( timings.mat ) );
    is.lasso <- grep( "lasso", colnames( timings.mat ) );

    # Filter out if we're not using all the methods.
    .exclude.cols <- c();
    if( !( "glm" %in% include.methods ) ) {
        .exclude.cols <- is.glm;
    }
    if( !( "step" %in% include.methods ) ) {
        .exclude.cols <- c( .exclude.cols, is.step );
    }
    if( !( "lasso" %in% include.methods ) ) {
        .exclude.cols <- c( .exclude.cols, is.lasso );
    }
    if( length( .exclude.cols ) > 0 ) {
        timings.mat <- timings.mat[ , -.exclude.cols, drop = FALSE ];
        
        is.glm <- grep( "glm", colnames( timings.mat ) );
        is.step <- grep( "step", colnames( timings.mat ) );
        is.lasso <- grep( "lasso", colnames( timings.mat ) );
    }
    has.nflg <- grep( "nflg\\.", colnames( timings.mat ) );
    has.v3 <- grep( "v3\\.", colnames( timings.mat ) );
    has.1m <- grep( "mtn003|1m\\.", colnames( timings.mat ) );
    has.6m <- grep( "hvtn502|\\.6m", colnames( timings.mat ) );
    with.bounds <- grep( "withbounds", colnames( timings.mat ) );
    
    # Rename cols:
    colnames( timings.mat ) <- gsub( "rmse", "RMSE", colnames( timings.mat ) );
    
    new.colnames <- colnames( timings.mat );
    for( col.i in 1:length( new.colnames ) ) {
        ## TODO: Put this information onto the plot by creating new rows.
        if( col.i %in% has.nflg ) {
            if( col.i %in% has.v3 ) {
                region.part <- "R";#"both regions";
            } else {
                region.part <- "r"; #"same region";
            }
        } else {
            stopifnot( col.i %in% has.v3 );
            region.part <- "r";#"same region";
        }

        if( col.i %in% has.1m ) {
            if( col.i %in% has.6m ) {
                time.part <- "T";#"both times";
            } else {
                time.part <- "t";#"same time";
            }
        } else {
            stopifnot( col.i %in% has.6m );
            time.part <- "t";#"same time";
        }
        if( col.i %in% c( is.glm, is.step, is.lasso ) ) {
            if( col.i %in% with.bounds ) {
                covars.part <- "[bound]+";
            } else {
                covars.part <- "";
            }
            if( col.i %in% is.step ) {
              if( col.i %in% with.bounds ) {
                covars.part <- paste( covars.part, "[step]+", sep = "" );
              } else { # if no bounds, don't need to mention other covars, as all are included.
                covars.part <- "[step]+";
              }
            }
            if( col.i %in% is.lasso ) {
              if( col.i %in% with.bounds ) {
                covars.part <- paste( covars.part, "[lasso]+", sep = "" );
              } else { # if no bounds, don't need to mention other covars, as all are included.
                covars.part <- "[lasso]+";
              }
            }
            if( include.intercept ) {
                intercept.part <- "1+";
            } else {
                intercept.part <- "0+";
            }
            type.part <- paste( "days~", intercept.part, covars.part, sep = "" );
        } else {
            type.part <- "uncalibrated";
        }
        new.colname <- paste( paste( "[", region.part, time.part, "]", sep = "" ), type.part );
        #new.colname <- type.part;
        #print( new.colname );
        new.colnames[ col.i ] <- new.colname;
    }
    #colnames( timings.mat ) <- gsub( "(nflg\\.)?(v3\\.)?(1m\\.)?(6m\\.)?", "", colnames( timings.mat ) );
    colnames( timings.mat ) <- new.colnames;

    ## The "uncalibrated" cols should all be identical; keep only the first.
    .cols <- grep( "uncalibrated", colnames( timings.mat ) );
    if( length( .cols ) > 1 ) {
        timings.mat.uncalibrated <- timings.mat[ , .cols, drop = FALSE ];
        stopifnot( all( apply( timings.mat.uncalibrated, 1, function( .row ) { all( .row == .row[1], na.rm = TRUE ); } ) ) );
        colnames( timings.mat )[ .cols[ 1 ] ] <- "uncalibrated";
        timings.mat <- timings.mat[ , -.cols[ -1 ], drop = FALSE ]; # keep only the first.
    }

    cols.with.included.training.code <- unlist( lapply( include.training.codes, function( .training.code ) { grep( paste( "\\[", .training.code, "\\]", sep = "" ), colnames( timings.mat ) ) } ) );
    
    # Filter down to the requested training codes. Always include the first column, which is uncalibrated.
    timings.mat <- timings.mat[ , c( 1, cols.with.included.training.code ), drop = FALSE ];
    timings.mat <- clean.timings.mat( timings.mat );
    
    return( timings.mat );
} # prepare.timings.mat ( .. )

prepare.usage.mat <- function ( usage.mat ) {
    # Rename cols:
    colnames( usage.mat )[ grep( "upper", colnames( usage.mat ) ) ] <- "[bound]";
    
    usage.mat <- clean.timings.mat( usage.mat );
    
    return( usage.mat );
} # prepare.usage.mat ( .. )

getTimingsResults <- function (
  include.intercept = FALSE,
  mutation.rate.calibration = FALSE,
  include.all.vars.in.lasso = TRUE,
  helpful.additional.cols = c( "lPVL" ),
  helpful.additional.cols.with.interactions = c(),#c( "v3_not_nflg", "X6m.not.1m" );
  use.gold.is.multiple = FALSE,
  force.recomputation = FORCE.RECOMPUTATION
) {
    # Ensure results exist.
    .results.by.region.and.time.Rda.filename <- evaluateTimings(
         include.intercept = include.intercept,
         include.all.vars.in.lasso = include.all.vars.in.lasso,
         helpful.additional.cols = helpful.additional.cols,
         helpful.additional.cols.with.interactions = helpful.additional.cols.with.interactions,
         use.gold.is.multiple = use.gold.is.multiple,
         force.recomputation = force.recomputation,
        mutation.rate.calibration = mutation.rate.calibration,
         RESULTS.DIR = RESULTS.DIR,
         results.dirname = RESULTS.DIRNAME
    );
    
    # evaluateTimings.compute.config.string(..) is defined in evaluateTimings.R.
    config.string <- evaluateTimings.compute.config.string(
            include.intercept = include.intercept,
            include.all.vars.in.lasso = include.all.vars.in.lasso,
            helpful.additional.cols = helpful.additional.cols,
            helpful.additional.cols.with.interactions = helpful.additional.cols.with.interactions,
            use.gold.is.multiple = use.gold.is.multiple
    );
    
    ### evaluateTimings
        if( config.string == "" ) {
            evaluateTimings.tab.file.suffix <- "_evaluateTimings.tab";
        } else {
            evaluateTimings.tab.file.suffix <- paste( "_evaluateTimings_", config.string, ".tab", sep = "" );
        }
        
    
    ###
    ## sampledwidth_uniform_1mmtn003_6mhvtn502 1m.v3
    evaluateTimings.nflg.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 <- getFilteredResultsTables( evaluateTimings.tab.file.suffix, the.region = "nflg", to.region = "v3", the.time = "1m.6m", sort.column = NULL, the.bounds.type = "sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3", column.pattern = "rmse", RESULTS.DIR = RESULTS.DIR, results.dirname = RESULTS.DIRNAME );
    colnames( evaluateTimings.nflg.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 ) <- paste( "nflg.v3.1m.6m", colnames( evaluateTimings.nflg.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 ), sep = "." );
    evaluateTimings.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 <- getFilteredResultsTables( evaluateTimings.tab.file.suffix, the.region = "v3", the.time = "1m.6m", sort.column = NULL, the.bounds.type = "sampledwidth_uniform_1mmtn003_6mhvtn502.1m", column.pattern = "rmse", RESULTS.DIR = RESULTS.DIR, results.dirname = RESULTS.DIRNAME );
    colnames( evaluateTimings.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 ) <- paste( "v3.1m.6m", colnames( evaluateTimings.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 ), sep = "." );
    evaluateTimings.nflg.v3.1m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 <- getFilteredResultsTables( evaluateTimings.tab.file.suffix, the.region = "nflg", to.region = "v3", the.time = "1m", sort.column = NULL, the.bounds.type = "sampledwidth_uniform_mtn003.v3", column.pattern = "rmse", RESULTS.DIR = RESULTS.DIR, results.dirname = RESULTS.DIRNAME );
    colnames( evaluateTimings.nflg.v3.1m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 ) <- paste( "nflg.v3.1m", colnames( evaluateTimings.nflg.v3.1m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 ), sep = "." );
    evaluateTimings.v3.1m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 <- getFilteredResultsTables( evaluateTimings.tab.file.suffix, the.region = "v3", the.time = "1m", sort.column = NULL, the.bounds.type = "sampledwidth_uniform_mtn003", column.pattern = "rmse", RESULTS.DIR = RESULTS.DIR, results.dirname = RESULTS.DIRNAME );
    colnames( evaluateTimings.v3.1m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 ) <- paste( "v3.1m", colnames( evaluateTimings.v3.1m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 ), sep = "." );
    
    evaluateTimings.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 <- missing.row.safe.cbind( missing.row.safe.cbind( missing.row.safe.cbind( evaluateTimings.nflg.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3, evaluateTimings.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 ), evaluateTimings.nflg.v3.1m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 ), evaluateTimings.v3.1m.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 );
    ###
    
    ###
    evaluateTimings.nflg.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 <- getFilteredResultsTables( evaluateTimings.tab.file.suffix, the.region = "nflg", to.region = "v3", the.time = "1m.6m", sort.column = NULL, the.bounds.type = "sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3", column.pattern = "rmse", RESULTS.DIR = RESULTS.DIR, results.dirname = RESULTS.DIRNAME );
    colnames( evaluateTimings.nflg.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 ) <- paste( "nflg.v3.1m.6m", colnames( evaluateTimings.nflg.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 ), sep = "." );
    evaluateTimings.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 <- getFilteredResultsTables( evaluateTimings.tab.file.suffix, the.region = "v3", the.time = "1m.6m", sort.column = NULL, the.bounds.type = "sampledwidth_uniform_1mmtn003_6mhvtn502.6m", column.pattern = "rmse", RESULTS.DIR = RESULTS.DIR, results.dirname = RESULTS.DIRNAME );
    colnames( evaluateTimings.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 ) <- paste( "v3.1m.6m", colnames( evaluateTimings.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 ), sep = "." );
    evaluateTimings.nflg.v3.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 <- getFilteredResultsTables( evaluateTimings.tab.file.suffix, the.region = "nflg", to.region = "v3", the.time = "6m", sort.column = NULL, the.bounds.type = "sampledwidth_uniform_hvtn502.v3", column.pattern = "rmse", RESULTS.DIR = RESULTS.DIR, results.dirname = RESULTS.DIRNAME );
    colnames( evaluateTimings.nflg.v3.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 ) <- paste( "nflg.v3.6m", colnames( evaluateTimings.nflg.v3.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 ), sep = "." );
    evaluateTimings.v3.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 <- getFilteredResultsTables( evaluateTimings.tab.file.suffix, the.region = "v3", the.time = "6m", sort.column = NULL, the.bounds.type = "sampledwidth_uniform_hvtn502", column.pattern = "rmse", RESULTS.DIR = RESULTS.DIR, results.dirname = RESULTS.DIRNAME );
    colnames( evaluateTimings.v3.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 ) <- paste( "v3.6m", colnames( evaluateTimings.v3.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 ), sep = "." );
    
    evaluateTimings.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 <- missing.row.safe.cbind( missing.row.safe.cbind( missing.row.safe.cbind( evaluateTimings.nflg.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3, evaluateTimings.v3.1m.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 ), evaluateTimings.nflg.v3.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 ), evaluateTimings.v3.6m.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 );
    ###
    
    timings.mat.evaluateTimings.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3 <-
        prepare.timings.mat( evaluateTimings.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3, include.intercept = include.intercept );
    timings.mat.evaluateTimings.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 <-
        prepare.timings.mat( evaluateTimings.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3, include.intercept = include.intercept );
    timings.mat.evaluateTimings.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.nflg <-
        prepare.timings.mat( evaluateTimings.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.nflg, include.intercept = include.intercept );
    timings.mat.evaluateTimings.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.nflg <-
        prepare.timings.mat( evaluateTimings.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.nflg, include.intercept = include.intercept );
    
    return( list( "v3" =
                      list( "1m" = timings.mat.evaluateTimings.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.v3,
                           "6m" = timings.mat.evaluateTimings.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.v3 ),
                 "nflg" =
                      list( "1m" = timings.mat.evaluateTimings.sampledwidth_uniform_1mmtn003_6mhvtn502.1m.nflg,
                           "6m" = timings.mat.evaluateTimings.sampledwidth_uniform_1mmtn003_6mhvtn502.6m.nflg )
                 )
           );

} # getTimingsResults (..)


prepare.summary.timings.mat <- function ( include.training.codes, the.region, the.time, include.V3 = FALSE, include.6m = FALSE, include.V3.6m = FALSE ) {
    timings.mat.nointercept.nolPVL.noV3.no6m.noGIM <-
        timings.mats.nointercept.nolPVL.noV3.no6m.noGIM[[ the.region ]][[ the.time ]];
    cols.with.included.training.code <- unlist( lapply( include.training.codes, function( .training.code ) { grep( paste( "\\[", .training.code, "\\]", sep = "" ), colnames( timings.mat.nointercept.nolPVL.noV3.no6m.noGIM ) ) } ) );
    # Filter down to the requested training codes.
    # The first column, which is uncalibrated, is the same in all of these matrices.
    timings.mat.uncalibrated.column <- 
        timings.mat.nointercept.nolPVL.noV3.no6m.noGIM[ , 1, drop = FALSE ];
    timings.mat.nointercept.nolPVL.noV3.no6m.noGIM <-
        timings.mat.nointercept.nolPVL.noV3.no6m.noGIM[ , cols.with.included.training.code, drop = FALSE ];
    # Remove the training code.
    colnames( timings.mat.nointercept.nolPVL.noV3.no6m.noGIM ) <-
        gsub( paste( "\\[", paste( include.training.codes, collapse = "|" ), "\\] ", sep = "" ), "", colnames( timings.mat.nointercept.nolPVL.noV3.no6m.noGIM ) );

    # Note that the cell without an intercept for "none" is actually the same as the one with an intercept, by convention. Here we replace it with NA.
    timings.mat.nointercept.nolPVL.noV3.no6m.noGIM[ "none", 1 ] <- NA;
    
    timings.mat.nointercept.yeslPVL.noV3.no6m.noGIM <-
        timings.mats.nointercept.yeslPVL.noV3.no6m.noGIM[[ the.region ]][[ the.time ]];
    cols.with.included.training.code <- unlist( lapply( include.training.codes, function( .training.code ) { grep( paste( "\\[", .training.code, "\\]", sep = "" ), colnames( timings.mat.nointercept.yeslPVL.noV3.no6m.noGIM ) ) } ) );
    timings.mat.nointercept.yeslPVL.noV3.no6m.noGIM <-
        timings.mat.nointercept.yeslPVL.noV3.no6m.noGIM[ , cols.with.included.training.code, drop = FALSE ];
    # Remove the training code.
    colnames( timings.mat.nointercept.yeslPVL.noV3.no6m.noGIM ) <-
        gsub( paste( "\\[", paste( include.training.codes, collapse = "|" ), "\\] ", sep = "" ), "", colnames( timings.mat.nointercept.yeslPVL.noV3.no6m.noGIM ) );
    # Add "lPVL+"
    colnames( timings.mat.nointercept.yeslPVL.noV3.no6m.noGIM ) <-
        gsub( "\\+$", "+lPVL+", colnames( timings.mat.nointercept.yeslPVL.noV3.no6m.noGIM ) );

    timings.mat.yesintercept.nolPVL.noV3.no6m.noGIM <-
        timings.mats.yesintercept.nolPVL.noV3.no6m.noGIM[[ the.region ]][[ the.time ]];
    cols.with.included.training.code <- unlist( lapply( include.training.codes, function( .training.code ) { grep( paste( "\\[", .training.code, "\\]", sep = "" ), colnames( timings.mat.yesintercept.nolPVL.noV3.no6m.noGIM ) ) } ) );
    timings.mat.yesintercept.nolPVL.noV3.no6m.noGIM <-
        timings.mat.yesintercept.nolPVL.noV3.no6m.noGIM[ , cols.with.included.training.code, drop = FALSE ];
    # Remove the training code.
    colnames( timings.mat.yesintercept.nolPVL.noV3.no6m.noGIM ) <-
        gsub( paste( "\\[", paste( include.training.codes, collapse = "|" ), "\\] ", sep = "" ), "", colnames( timings.mat.yesintercept.nolPVL.noV3.no6m.noGIM ) );
    
    timings.mat.yesintercept.yeslPVL.noV3.no6m.noGIM <-
        timings.mats.yesintercept.yeslPVL.noV3.no6m.noGIM[[ the.region ]][[ the.time ]];
    cols.with.included.training.code <- unlist( lapply( include.training.codes, function( .training.code ) { grep( paste( "\\[", .training.code, "\\]", sep = "" ), colnames( timings.mat.yesintercept.yeslPVL.noV3.no6m.noGIM ) ) } ) );
    timings.mat.yesintercept.yeslPVL.noV3.no6m.noGIM <-
        timings.mat.yesintercept.yeslPVL.noV3.no6m.noGIM[ , cols.with.included.training.code, drop = FALSE ];
    # Remove the training code.
    colnames( timings.mat.yesintercept.yeslPVL.noV3.no6m.noGIM ) <-
        gsub( paste( "\\[", paste( include.training.codes, collapse = "|" ), "\\] ", sep = "" ), "", colnames( timings.mat.yesintercept.yeslPVL.noV3.no6m.noGIM ) );
    # Add "lPVL+"
    colnames( timings.mat.yesintercept.yeslPVL.noV3.no6m.noGIM ) <-
        gsub( "\\+$", "+lPVL+", colnames( timings.mat.yesintercept.yeslPVL.noV3.no6m.noGIM ) );
    
    timings.mat <-
        cbind( timings.mat.uncalibrated.column,
              timings.mat.yesintercept.nolPVL.noV3.no6m.noGIM,
              timings.mat.yesintercept.yeslPVL.noV3.no6m.noGIM,
              timings.mat.nointercept.nolPVL.noV3.no6m.noGIM,
              timings.mat.nointercept.yeslPVL.noV3.no6m.noGIM
        );
    
    if( include.V3 ) {
        timings.mat.nointercept.nolPVL.yesV3.no6m.noGIM <-
            timings.mats.nointercept.nolPVL.yesV3.no6m.noGIM[[ the.region ]][[ the.time ]];
        cols.with.included.training.code <- unlist( lapply( include.training.codes, function( .training.code ) { grep( paste( "\\[", .training.code, "\\]", sep = "" ), colnames( timings.mat.nointercept.nolPVL.yesV3.no6m.noGIM ) ) } ) );
        # Filter down to the requested training codes.
        # The first column, which is uncalibrated, is the same in all of these matrices.
        timings.mat.uncalibrated.column <- 
            timings.mat.nointercept.nolPVL.yesV3.no6m.noGIM[ , 1, drop = FALSE ];
        timings.mat.nointercept.nolPVL.yesV3.no6m.noGIM <-
            timings.mat.nointercept.nolPVL.yesV3.no6m.noGIM[ , cols.with.included.training.code, drop = FALSE ];
        # Remove the training code.
        colnames( timings.mat.nointercept.nolPVL.yesV3.no6m.noGIM ) <-
            gsub( paste( "\\[", paste( include.training.codes, collapse = "|" ), "\\] ", sep = "" ), "", colnames( timings.mat.nointercept.nolPVL.yesV3.no6m.noGIM ) );
        # Add "isV3+"
        colnames( timings.mat.nointercept.nolPVL.yesV3.no6m.noGIM ) <-
            gsub( "\\+$", "+isV3+", colnames( timings.mat.nointercept.nolPVL.yesV3.no6m.noGIM ) );
        
        timings.mat.nointercept.yeslPVL.yesV3.no6m.noGIM <-
            timings.mats.nointercept.yeslPVL.yesV3.no6m.noGIM[[ the.region ]][[ the.time ]];
        cols.with.included.training.code <- unlist( lapply( include.training.codes, function( .training.code ) { grep( paste( "\\[", .training.code, "\\]", sep = "" ), colnames( timings.mat.nointercept.yeslPVL.yesV3.no6m.noGIM ) ) } ) );
        timings.mat.nointercept.yeslPVL.yesV3.no6m.noGIM <-
            timings.mat.nointercept.yeslPVL.yesV3.no6m.noGIM[ , cols.with.included.training.code, drop = FALSE ];
        # Remove the training code.
        colnames( timings.mat.nointercept.yeslPVL.yesV3.no6m.noGIM ) <-
            gsub( paste( "\\[", paste( include.training.codes, collapse = "|" ), "\\] ", sep = "" ), "", colnames( timings.mat.nointercept.yeslPVL.yesV3.no6m.noGIM ) );
        # Add "isV3+"
        colnames( timings.mat.nointercept.yeslPVL.yesV3.no6m.noGIM ) <-
            gsub( "\\+$", "+isV3+", colnames( timings.mat.nointercept.yeslPVL.yesV3.no6m.noGIM ) );
        # Add "lPVL+"
        colnames( timings.mat.nointercept.yeslPVL.yesV3.no6m.noGIM ) <-
            gsub( "\\+$", "+lPVL+", colnames( timings.mat.nointercept.yeslPVL.yesV3.no6m.noGIM ) );
        timings.mat <- cbind( timings.mat,
                             timings.mat.nointercept.nolPVL.yesV3.no6m.noGIM,
                             timings.mat.nointercept.yeslPVL.yesV3.no6m.noGIM );
    } # End if include.V3    
    if( include.6m ) {
        timings.mat.nointercept.nolPVL.noV3.yes6m.noGIM <-
            timings.mats.nointercept.nolPVL.noV3.yes6m.noGIM[[ the.region ]][[ the.time ]];
        cols.with.included.training.code <- unlist( lapply( include.training.codes, function( .training.code ) { grep( paste( "\\[", .training.code, "\\]", sep = "" ), colnames( timings.mat.nointercept.nolPVL.noV3.yes6m.noGIM ) ) } ) );
        # Filter down to the requested training codes.
        # The first column, which is uncalibrated, is the same in all of these matrices.
        timings.mat.uncalibrated.column <- 
            timings.mat.nointercept.nolPVL.noV3.yes6m.noGIM[ , 1, drop = FALSE ];
        timings.mat.nointercept.nolPVL.noV3.yes6m.noGIM <-
            timings.mat.nointercept.nolPVL.noV3.yes6m.noGIM[ , cols.with.included.training.code, drop = FALSE ];
        # Remove the training code.
        colnames( timings.mat.nointercept.nolPVL.noV3.yes6m.noGIM ) <-
            gsub( paste( "\\[", paste( include.training.codes, collapse = "|" ), "\\] ", sep = "" ), "", colnames( timings.mat.nointercept.nolPVL.noV3.yes6m.noGIM ) );
        # Add "is6m+"
        colnames( timings.mat.nointercept.nolPVL.noV3.yes6m.noGIM ) <-
            gsub( "\\+$", "+is6m+", colnames( timings.mat.nointercept.nolPVL.noV3.yes6m.noGIM ) );
        
        timings.mat.nointercept.yeslPVL.noV3.yes6m.noGIM <-
            timings.mats.nointercept.yeslPVL.noV3.yes6m.noGIM[[ the.region ]][[ the.time ]];
        cols.with.included.training.code <- unlist( lapply( include.training.codes, function( .training.code ) { grep( paste( "\\[", .training.code, "\\]", sep = "" ), colnames( timings.mat.nointercept.yeslPVL.noV3.yes6m.noGIM ) ) } ) );
        timings.mat.nointercept.yeslPVL.noV3.yes6m.noGIM <-
            timings.mat.nointercept.yeslPVL.noV3.yes6m.noGIM[ , cols.with.included.training.code, drop = FALSE ];
        # Remove the training code.
        colnames( timings.mat.nointercept.yeslPVL.noV3.yes6m.noGIM ) <-
            gsub( paste( "\\[", paste( include.training.codes, collapse = "|" ), "\\] ", sep = "" ), "", colnames( timings.mat.nointercept.yeslPVL.noV3.yes6m.noGIM ) );
        # Add "is6m+"
        colnames( timings.mat.nointercept.yeslPVL.noV3.yes6m.noGIM ) <-
            gsub( "\\+$", "+is6m+", colnames( timings.mat.nointercept.yeslPVL.noV3.yes6m.noGIM ) );
        # Add "lPVL+"
        colnames( timings.mat.nointercept.yeslPVL.noV3.yes6m.noGIM ) <-
            gsub( "\\+$", "+lPVL+", colnames( timings.mat.nointercept.yeslPVL.noV3.yes6m.noGIM ) );
        timings.mat <- cbind( timings.mat,
                             timings.mat.nointercept.nolPVL.noV3.yes6m.noGIM,
                             timings.mat.nointercept.yeslPVL.noV3.yes6m.noGIM );
    } # End if include.6m    
    if( include.V3.6m ) {
        timings.mat.nointercept.nolPVL.yesV3.yes6m.noGIM <-
            timings.mats.nointercept.nolPVL.yesV3.yes6m.noGIM[[ the.region ]][[ the.time ]];
        cols.with.included.training.code <- unlist( lapply( include.training.codes, function( .training.code ) { grep( paste( "\\[", .training.code, "\\]", sep = "" ), colnames( timings.mat.nointercept.nolPVL.yesV3.yes6m.noGIM ) ) } ) );
        # Filter down to the requested training codes.
        # The first column, which is uncalibrated, is the same in all of these matrices.
        timings.mat.uncalibrated.column <- 
            timings.mat.nointercept.nolPVL.yesV3.yes6m.noGIM[ , 1, drop = FALSE ];
        timings.mat.nointercept.nolPVL.yesV3.yes6m.noGIM <-
            timings.mat.nointercept.nolPVL.yesV3.yes6m.noGIM[ , cols.with.included.training.code, drop = FALSE ];
        # Remove the training code.
        colnames( timings.mat.nointercept.nolPVL.yesV3.yes6m.noGIM ) <-
            gsub( paste( "\\[", paste( include.training.codes, collapse = "|" ), "\\] ", sep = "" ), "", colnames( timings.mat.nointercept.nolPVL.yesV3.yes6m.noGIM ) );
        # Add "isV3+is6m+"
        colnames( timings.mat.nointercept.nolPVL.yesV3.yes6m.noGIM ) <-
            gsub( "\\+$", "+isV3+is6m+", colnames( timings.mat.nointercept.nolPVL.yesV3.yes6m.noGIM ) );
        
        timings.mat.nointercept.yeslPVL.yesV3.yes6m.noGIM <-
            timings.mats.nointercept.yeslPVL.yesV3.yes6m.noGIM[[ the.region ]][[ the.time ]];
        cols.with.included.training.code <- unlist( lapply( include.training.codes, function( .training.code ) { grep( paste( "\\[", .training.code, "\\]", sep = "" ), colnames( timings.mat.nointercept.yeslPVL.yesV3.yes6m.noGIM ) ) } ) );
        timings.mat.nointercept.yeslPVL.yesV3.yes6m.noGIM <-
            timings.mat.nointercept.yeslPVL.yesV3.yes6m.noGIM[ , cols.with.included.training.code, drop = FALSE ];
        # Remove the training code.
        colnames( timings.mat.nointercept.yeslPVL.yesV3.yes6m.noGIM ) <-
            gsub( paste( "\\[", paste( include.training.codes, collapse = "|" ), "\\] ", sep = "" ), "", colnames( timings.mat.nointercept.yeslPVL.yesV3.yes6m.noGIM ) );
        # Add "isV3+is6m+"
        colnames( timings.mat.nointercept.yeslPVL.yesV3.yes6m.noGIM ) <-
            gsub( "\\+$", "+isV3+is6m+", colnames( timings.mat.nointercept.yeslPVL.yesV3.yes6m.noGIM ) );
        # Add "lPVL+"
        colnames( timings.mat.nointercept.yeslPVL.yesV3.yes6m.noGIM ) <-
            gsub( "\\+$", "+lPVL+", colnames( timings.mat.nointercept.yeslPVL.yesV3.yes6m.noGIM ) );
        timings.mat <- cbind( timings.mat,
                             timings.mat.nointercept.nolPVL.yesV3.yes6m.noGIM,
                             timings.mat.nointercept.yeslPVL.yesV3.yes6m.noGIM );
    } # End if include.V3.6m    
    
    return( timings.mat );     
} # prepare.summary.timings.mat ( .. )

get.days.since.infection <-
          function ( results.by.region.and.time, regions = c( "nflg", "v3" ), times = c( "1m", "6m" ) ) {
            if( length( times ) == 2 ) {
                the.time <- "1m.6m";
            } else {
                the.time <- times;
            }
            if( length( regions ) == 2 ) {
                .results.for.region <- results.by.region.and.time[[3]][[1]][[1]];
            } else {
                .results.for.region <- results.by.region.and.time[[ regions ]];
            }
            results.covars.per.person.df <-
                data.frame( .results.for.region[[ the.time ]][[ "results.covars.per.person.with.extra.cols" ]] );

            return( .results.for.region[[ the.time ]][["days.since.infection" ]][ rownames( results.covars.per.person.df ) ] );
} # get.days.since.infection (..)

getTimingsResultFormula <- function (
  varname = "PFitter.mut.rate.coef",
  withbounds = FALSE,
  regions = c( "nflg", "v3" ),
  times = c( "1m", "6m" ),
  include.intercept = FALSE,
  mutation.rate.calibration = FALSE,
  include.all.vars.in.lasso = TRUE,
  helpful.additional.cols = c( "lPVL" ),
  helpful.additional.cols.with.interactions = c(),#c( "v3_not_nflg", "X6m.not.1m" );
  use.gold.is.multiple = FALSE,
  force.recomputation = FORCE.RECOMPUTATION,#FALSE,
  evaluate.regions = regions,
  evaluate.times = times
) {

# Ensure results exist.
.results.by.region.and.time.Rda.filename <- evaluateTimings(
     include.intercept = include.intercept,
     include.all.vars.in.lasso = include.all.vars.in.lasso,
     helpful.additional.cols = helpful.additional.cols,
     helpful.additional.cols.with.interactions = helpful.additional.cols.with.interactions,
     use.gold.is.multiple = use.gold.is.multiple,
     force.recomputation = force.recomputation,
     mutation.rate.calibration = mutation.rate.calibration,
     RESULTS.DIR = RESULTS.DIR,
     results.dirname = RESULTS.DIRNAME
);

    if( length( times ) == 2 ) { 
        the.bounds.type <- "sampledwidth_uniform_1mmtn003_6mhvtn502";
    } else if( times == "1m" ) {
        the.bounds.type <- "sampledwidth_uniform_mtn003";
    } else {
        stopifnot( times == "6m" );
        the.bounds.type <- "sampledwidth_uniform_hvtn502";
    }
    
    
#print( .results.by.region.and.time.Rda.filename );
local({
        load( .results.by.region.and.time.Rda.filename ); # adds results.by.region.and.time to environment.
        .formulas.table <- get.formulas( results.by.region.and.time, .varname = varname, model.type = "glm", withbounds = withbounds, regions = regions, times = times, the.bounds.type = the.bounds.type );
        .evaluated.formula <-
            evaluate.specific.timings.model.formula( results.by.region.and.time, names( which.max( .formulas.table ) ), step = FALSE, regions = evaluate.regions, times = evaluate.times );
        .evaluated.formula.table <- cbind( coef( summary( .evaluated.formula ) )[ , "Estimate", drop = FALSE ], confint( .evaluated.formula ), coef( summary( .evaluated.formula ) )[ , "Pr(>|t|)", drop = FALSE ] );
        return( list( formulas = .formulas.table, evaluated.formula = summary( .evaluated.formula ), evaluated.formula.table = .evaluated.formula.table ) );
     });
} # getTimingsResultFormula (..)

createMutationRateScalarPlot <- function (
  varname = "PFitter.mut.rate.coef",
  varname.pretty = "PFitter",
  force.recomputation = FORCE.RECOMPUTATION,
  evaluate.regions = regions,
  evaluate.times = times
) {

    timings.result.formula.v3.6m <- 
        getTimingsResultFormula( varname, include.intercept = FALSE, helpful.additional.cols = c(), force.recomputation = FORCE.RECOMPUTATION, regions = "v3", times = "6m" );
    
    timings.result.formula.nflg.6m <- 
        getTimingsResultFormula( varname, include.intercept = FALSE, helpful.additional.cols = c(), force.recomputation = FORCE.RECOMPUTATION, regions = "nflg", times = "6m" );
    timings.result.formula.nflgv3.6m <- 
        getTimingsResultFormula( varname, include.intercept = FALSE, helpful.additional.cols = c(), force.recomputation = FORCE.RECOMPUTATION, regions = c( "nflg", "v3" ), times = "6m" );
    
    timings.result.formula.v3.1m6m <- 
        getTimingsResultFormula( varname, include.intercept = FALSE, helpful.additional.cols = c(), force.recomputation = FORCE.RECOMPUTATION, regions = c( "v3" ), times = c( "1m", "6m" ) );
    
    timings.result.formula.nflgv3.1m6m <- 
        getTimingsResultFormula( varname, include.intercept = FALSE, helpful.additional.cols = c(), force.recomputation = FORCE.RECOMPUTATION, regions = c( "nflg", "v3" ), times = c( "1m", "6m" ) );
    
    ## 1m
    timings.result.formula.v3.1m <- 
        getTimingsResultFormula( varname, include.intercept = FALSE, helpful.additional.cols = c(), force.recomputation = FORCE.RECOMPUTATION, regions = "v3", times = "1m" );
    
    timings.result.formula.nflg.1m <- 
        getTimingsResultFormula( varname, include.intercept = FALSE, helpful.additional.cols = c(), force.recomputation = FORCE.RECOMPUTATION, regions = "nflg", times = "1m" );
    timings.result.formula.nflg.1m6m <- 
        getTimingsResultFormula( varname, include.intercept = FALSE, helpful.additional.cols = c(), force.recomputation = FORCE.RECOMPUTATION, regions = "nflg", times = c( "1m", "6m" ) );
    
    timings.result.formula.nflgv3.1m <- 
        getTimingsResultFormula( varname, include.intercept = FALSE, helpful.additional.cols = c(), force.recomputation = FORCE.RECOMPUTATION, regions = c( "nflg", "v3" ), times = "1m" );
    
    if( length( grep( "mut.rate.coef", varname ) ) > 0 ) {
        maybe.invert.function <- function ( x ) { 1 / x }
    } else {
        maybe.invert.function <- function ( x ) { x }
        stop( "this is designed to show PFitter's results. TODO: Make a version that uses a days estimate instead of a mut.rate.coef." );
    }
    
    ## Forestplot the updated PFitter mutation rates.
    .mat <- rbind(
        c( region = "v3", time = "1m", signif( maybe.invert.function( timings.result.formula.v3.1m$evaluated.formula.table[1,c(1,3,2)]), 3 ) ),
        c( region = "nflg", time = "1m", signif( maybe.invert.function( timings.result.formula.nflg.1m$evaluated.formula.table[1,c(1,3,2)]), 3 ) ),
        c( region = "nflg, v3", time = "1m", signif( maybe.invert.function( timings.result.formula.nflgv3.1m$evaluated.formula.table[1,c(1,3,2)]), 3 ) ),
        c( region = "v3", time = "6m", signif( maybe.invert.function( timings.result.formula.v3.6m$evaluated.formula.table[1,c(1,3,2)]), 3 ) ),
        c( region = "nflg", time = "6m", signif( maybe.invert.function( timings.result.formula.nflg.6m$evaluated.formula.table[1,c(1,3,2)]), 3 ) ),
        c( region = "nflg, v3", time = "6m", signif( maybe.invert.function( timings.result.formula.nflgv3.6m$evaluated.formula.table[1,c(1,3,2)]), 3 ) ),
        c( region = "v3", time = "1m, 6m", signif( maybe.invert.function( timings.result.formula.v3.1m6m$evaluated.formula.table[1,c(1,3,2)]), 3 ) ),
        c( region = "nflg", time = "1m, 6m", signif( maybe.invert.function( timings.result.formula.nflg.1m6m$evaluated.formula.table[1,c(1,3,2)]), 3 ) ),
        c( region = "nflg, v3", time = "1m, 6m", signif( maybe.invert.function( timings.result.formula.nflgv3.1m6m$evaluated.formula.table[1,c(1,3,2)]), 3 ) )
    );
    .mat.for.table <- .mat[ , 1:2, drop = FALSE ];
    .mat.for.table <- rbind( c( "Regions", "Samples" ), .mat.for.table ); # Add a row for the column headings
        
    .mat.for.CIs <- .mat[ , c( "Estimate", "2.5 %", "97.5 %" ) ];
    mode( .mat.for.CIs ) <- "numeric";
    .mat.for.CIs <- rbind( rep( NA, ncol( .mat.for.CIs ) ), .mat.for.CIs ); # Add a row for the column headings

    pdf( paste( varname.pretty, " Calibrated Mutation Rates forestplot.pdf", sep = FALSE ) );
    forestplot( .mat.for.table, (.mat.for.CIs*10^5)/1.19, is.summary = c( TRUE, rep( FALSE, nrow( .mat ) ) ),
               new_page = FALSE, 
               txt_gp = fpTxtGp(label = list(gpar(fontfamily = "",
                                                  col = "#660000"),
                                             gpar(fontfamily = "",
                                                  col = "#660000")),
                                ticks = gpar(fontfamily = "", cex=1),
                                xlab  = gpar(fontfamily = "HersheySerif", cex = 1.5)),
               xlab="Mutation rate (x 10^-5 substitutions per day)",
                   xlog=TRUE, 
                   col=fpColors(box="royalblue",line="darkblue", summary="royalblue"), 
                   vertices = TRUE, boxsize = 0.1, line.margin = 0.05, zero = 1.19 ); # note zero is the PFitter default.
    dev.off();

    pdf( paste( varname, " Coefficient forestplot.pdf", sep = "" ) );
    forestplot( .mat.for.table, (.mat.for.CIs*10^5), is.summary = c( TRUE, rep( FALSE, nrow( .mat ) ) ),
           new_page = FALSE, 
           txt_gp = fpTxtGp(label = list(gpar(fontfamily = "",
                                              col = "#660000"),
                                         gpar(fontfamily = "",
                                              col = "#660000")),
                            ticks = gpar(fontfamily = "", cex=1),
                            xlab  = gpar(fontfamily = "HersheySerif", cex = 1.5)),
           xlab=paste( "Multiple of ", varname.pretty, " Days Estimate", sep = "" ),
               xlog=TRUE, 
               col=fpColors(box="royalblue",line="darkblue", summary="royalblue"), 
               vertices = TRUE, boxsize = 0.1, line.margin = 0.05, zero = 1.19/1.19 ); # note zero is the PFitter default.
    dev.off();
} # createMutationRateScalarPlot (..)

get.predicted.values <-
    function ( results.by.region.and.time, varname, withbounds = FALSE, method = "glm", regions = c( "nflg", "v3" ), times = c( "1m", "6m" )  ) {
    if( method == "glm" ) {
        if( withbounds ) {
            the.prefix <- "glm.withbounds.validation.results.";
        } else {
            the.prefix <- "glm.validation.results.";
        }
    } else {
        stopifnot( method == "uncalibrated" );
        stopifnot( withbounds == FALSE );
        the.prefix <- "";
    }
        
    if( length( times ) == 2 ) { 
        the.bounds.type <- "sampledwidth_uniform_1mmtn003_6mhvtn502";
    } else if( times == "1m" ) {
        the.bounds.type <- "sampledwidth_uniform_mtn003";
    } else {
        stopifnot( times == "6m" );
        the.bounds.type <- "sampledwidth_uniform_hvtn502";
    }
    

            if( length( times ) == 2 ) {
                the.time <- "1m.6m";
            } else {
                the.time <- times;
            }
            if( length( regions ) == 2 ) {
                .results.for.region <- results.by.region.and.time[[3]][[1]][[1]];
            } else {
                .results.for.region <- results.by.region.and.time[[ regions ]];
            }
            results.covars.per.person.df <-
                data.frame( .results.for.region[[ the.time ]][[ "results.covars.per.person.with.extra.cols" ]] );

    if( !( paste( the.prefix, varname, sep = "" ) %in% colnames( .results.for.region[[ the.time ]][["evaluated.results"]][[the.bounds.type]][["results.per.person"]] ) ) ) {
        warning( paste( paste( the.prefix, varname, sep = "" ), "is not a valid column name; returning a column of NAs" ) );
        .rv <- .results.for.region[[ the.time ]][["evaluated.results"]][[the.bounds.type]][["results.per.person"]][rownames(results.covars.per.person.df), 1 ];
        .rv[ 1:length( .rv ) ] <- NA;
        return( .rv );
    }
    return( .results.for.region[[ the.time ]][["evaluated.results"]][[the.bounds.type]][["results.per.person"]][rownames(results.covars.per.person.df), paste( the.prefix, varname, sep = "" ) ] );
    } # get.predicted.values (..)
getV3Timings <- function (
                          ...
                          ) {
    getTimings( "v3", ... );
} # getV3Timings (..)
getNFLGTimings <- function (
                          ...
                          ) {
    getTimings( "nflg", ... );
} # getNFLGTimings (..)
getTimings <- function (
  the.region = "both", # region of data to get.
  regions = c( "nflg", "v3" ),
  times = c( "1m", "6m" ),
  include.intercept = FALSE,
  mutation.rate.calibration = FALSE,
  include.all.vars.in.lasso = TRUE,
  helpful.additional.cols = c(),
  helpful.additional.cols.with.interactions = c(),
  use.gold.is.multiple = FALSE,
  force.recomputation = FALSE
) {

    # Ensure results exist.
    .results.by.region.and.time.Rda.filename <- evaluateTimings(
         include.intercept = include.intercept,
         include.all.vars.in.lasso = include.all.vars.in.lasso,
         helpful.additional.cols = helpful.additional.cols,
         helpful.additional.cols.with.interactions = helpful.additional.cols.with.interactions,
         use.gold.is.multiple = use.gold.is.multiple,
         force.recomputation = force.recomputation,
         mutation.rate.calibration = mutation.rate.calibration,
         RESULTS.DIR = RESULTS.DIR,
         results.dirname = RESULTS.DIRNAME
    );
    
    #print( .results.by.region.and.time.Rda.filename );
    local({
            load( .results.by.region.and.time.Rda.filename ); # adds results.by.region.and.time to environment.
            days.since.infection <- get.days.since.infection( results.by.region.and.time, regions = regions, times = times );
            if( the.region == "both" ) {
                return( days.since.infection );
            }
            if( the.region == "v3" ) {
                days.since.infection.v3 <- days.since.infection[ grep( "nflg$", names( days.since.infection ), invert = TRUE ) ];
                names( days.since.infection.v3 ) <- gsub( "\\.v3$", "", names( days.since.infection.v3 ) );
                return( days.since.infection.v3 );
            }
            stopifnot( the.region == "nflg" );
            days.since.infection.nflg <- days.since.infection[ grep( "v3", names( days.since.infection ), invert = TRUE ) ];
            names( days.since.infection.nflg ) <- gsub( "\\.nflg$", "", names( days.since.infection.nflg ) );
            return( days.since.infection.nflg );
    } );
} # getTimings ( .. )

getV3Predictions <- function ( ... ) {
    getPredictions( "v3", ... )
} # getV3Predictions (...)
getNFLGPredictions <- function ( ... ) {
    getPredictions( "nflg", ... )
} # getNFLGPredictions (...)
getPredictions <- function (
  the.region = "both",
  varname = "none",
  withbounds = FALSE,
  method = "glm", # or "uncalibrated"
  regions = c( "nflg", "v3" ),
  times = c( "1m", "6m" ),
  include.intercept = FALSE,
  mutation.rate.calibration = FALSE,
  include.all.vars.in.lasso = TRUE,
  helpful.additional.cols = c(),
  helpful.additional.cols.with.interactions = c(),
  use.gold.is.multiple = FALSE,
  force.recomputation = FALSE
) {

    # Ensure results exist.
    .results.by.region.and.time.Rda.filename <- evaluateTimings(
         include.intercept = include.intercept,
         include.all.vars.in.lasso = include.all.vars.in.lasso,
         helpful.additional.cols = helpful.additional.cols,
         helpful.additional.cols.with.interactions = helpful.additional.cols.with.interactions,
         use.gold.is.multiple = use.gold.is.multiple,
         force.recomputation = force.recomputation,
        mutation.rate.calibration = mutation.rate.calibration,
         RESULTS.DIR = RESULTS.DIR,
         results.dirname = RESULTS.DIRNAME
    );
    
    #print( .results.by.region.and.time.Rda.filename );
    local({
            load( .results.by.region.and.time.Rda.filename ); # adds results.by.region.and.time to environment.
            predicted.values <-
                get.predicted.values( results.by.region.and.time, varname = varname, withbounds = withbounds, method = method, regions = regions, times = times );
            if( the.region == "both" ) {
                return( predicted.values );
            }
            if( the.region == "v3" ) {
                predicted.values.v3 <- predicted.values[ grep( "nflg$", names( predicted.values ), invert = TRUE ) ];
                names( predicted.values.v3 ) <- gsub( "\\.v3$", "", names( predicted.values.v3 ) );
                return( predicted.values.v3 );
            }
            stopifnot( the.region == "nflg" );
            predicted.values.nflg <- predicted.values[ grep( "v3$", names( predicted.values ), invert = TRUE ) ];
            names( predicted.values.nflg ) <- gsub( "\\.nflg$", "", names( predicted.values.nflg ) );
            return( predicted.values.nflg );
    } );
} # getPredictions ( .. )

makeGGPFromDataFrame <- function ( .the.data.frame, .Biases = NULL, .RMSEs = NULL, .Biases.and.RMSEs.yloc = NULL ) {
    .ggp <-
        ggplot( .the.data.frame, aes(x=Estimator,y=Days)) + geom_boxplot(aes(fill=Estimator) ) + guides(fill=FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      axis.text.x =
          element_text(
              angle = 90, vjust = 0.95, hjust=0.95,
              colour="black",size=14,face="bold"
           ),
      axis.text.y =
          element_text(
              colour="black",size=18,face="bold"
           ),
      axis.title.y =
          element_text(
              colour="black",size=18,face="bold"
           )
        );
    if( !is.null( .Biases ) || !is.null( .RMSEs ) ) {
        if( is.null( .Biases.and.RMSEs.yloc ) ) {
          if( "layout" %in% names( ggplot_build( .ggp ) ) ) {
              .yrange <- ggplot_build(.ggp)$layout$panel_ranges[[1]]$y.range;
          } else {
              .yrange <- ggplot_build(.ggp)$panel$ranges[[1]]$y.range;
          }
          # This is a magic #! Same as in addYLims: 1.2.
          .lim <- 1.2*max( abs( .yrange ) );
          .maxError <- max(.the.data.frame$Days);
          .yloc <- mean( c( .maxError, .lim ) );
        } else { 
            .yloc <- .Biases.and.RMSEs.yloc;
        }
        if( !is.null( .Biases ) && !is.null( .RMSEs ) ) {
            stopifnot( length( .Biases ) == length( .RMSEs ) );
            .Biases.and.RMSEs <- sapply( 1:length( .Biases ), function( .i ) {
                sprintf( "%s (%s)", .Biases[ .i ], .RMSEs[ .i ] );
             } );
            .ggp <- .ggp +
                annotate("text", x = 1:length(.Biases), y=.yloc, label = .Biases.and.RMSEs )
        } else if( !is.null( .Biases ) ) {
            .ggp <- .ggp +
                annotate("text", x = 1:length(.Biases), y=.yloc, label = .Biases)
        } else if( !is.null( .RMSEs ) ) {
            .ggp <- .ggp +
                annotate("text", x = 1:length(.RMSEs), y=.yloc, label = .RMSEs)
        }
    }
  #if( the.var.scale.is.discrete ) {
  #    .ggp <- .ggp + scale_y_discrete();
  #}
  .ggp <- .ggp + ylab( "Days" ) + xlab( "" );
return( .ggp );
    } # makeGGPFromDataFrame (..)

addYLims <- function ( .ggp, max.or.ggp.source.for.lims, y.low = NULL ) {
    if( class( max.or.ggp.source.for.lims ) == "numeric" ) {
        .lim <- max.or.ggp.source.for.lims;
    } else {
        .ggp.source.for.lims <- max.or.ggp.source.for.lims;
        if( "layout" %in% names( ggplot_build( .ggp.source.for.lims ) ) ) {
            .yrange <- ggplot_build(.ggp.source.for.lims)$layout$panel_ranges[[1]]$y.range;
        } else {
            .yrange <- ggplot_build(.ggp.source.for.lims)$panel$ranges[[1]]$y.range;
        }
        .lim <- 1.2*max( abs( .yrange ) );
    }
    if( !is.null( y.low ) ) {
        return( .ggp + ylim( y.low, .lim ) );
    }
    return( .ggp + ylim( -.lim, .lim ) );
} # addYLims (..)

calculateRMSEsFromPredictionErrors <- 
    function( .col ) { sprintf( "%0.1f", sqrt( mean( .col*.col, na.rm = T ) ) ) };    

calculateBiasesFromPredictionErrors <- 
    function( .col ) { sprintf( "%0.1f", mean( .col, na.rm = T ) ) };    

## This creates pdf and csv file outputs for the timings results
createResultsPlotAndSpreadsheet <- function ( include.intercept, withbounds, train.regions, region, time, mutation.rate.calibration = FALSE, ylims.1m = 300, ylims.6m = 300 ) {
    if( time == "1w" ) {
        ylims = ylims.1m;## TODO: Add a 1w version in the argument list?
        .Biases.and.RMSEs.yloc = -ylims.1m;
    } else if( time == "1m" ) {
        ylims = ylims.1m;
        .Biases.and.RMSEs.yloc = -ylims.1m;
    } else if( time == "6m" ) {
        ylims = ylims.6m;
        .Biases.and.RMSEs.yloc = -ylims.6m;
    } else {
        stop( "time should be 1m or 6m" );
    }

    gold.standard <- getTimings( region, include.intercept = FALSE, helpful.additional.cols = c(), regions = region, times = time );

    .config.string <- paste( region, " Days at ", time, sep = "" );
    if( withbounds ) {
        .config.string <- paste( .config.string, " with bound", sep = "" );
    }
    if( include.intercept ) {
        .config.string <- paste( .config.string, " including intercept", sep = "" );
    }
    if( !( ( length( the.region ) == length( train.regions ) ) && ( all( sort( the.region ) == sort( train.regions ) ) ) ) ) {
        .config.string <- paste( .config.string, " trained using regions ", paste( train.regions, collapse = ", " ), sep = "" );
    }

    bounded.COB.var <- ifelse( time == "1m", "COB.sampledwidth.uniform.mtn003.time.est", "COB.sampledwidth.uniform.hvtn502.time.est" );
    bounded.PrankenBeast.var <- ifelse( time == "1m", "Infer.sampledwidth.uniform.mtn003.time.est", "Infer.sampledwidth.uniform.hvtn502.time.est" );

    predicted.values.calibrated.none <- getPredictions( region, "none", include.intercept = include.intercept, mutation.rate.calibration = mutation.rate.calibration, withbounds = withbounds, helpful.additional.cols = c( "lPVL" ), regions = train.regions, times = time, method = "glm" );
#    formula.calibrated.withintercept.withlPVL.withBounds.none <- getTimingsResultFormula( "none", withbounds = TRUE, include.intercept = TRUE, helpful.additional.cols = c("lPVL"), regions = region, times = time );
#    formula.calibrated.nointercept.withlPVL.withBounds.none <- getTimingsResultFormula( "none", withbounds = TRUE, include.intercept = FALSE, helpful.additional.cols = c("lPVL"), regions = region, times = time );
#    formula.calibrated.nointercept.withlPVL.noBounds.none <- getTimingsResultFormula( "none", withbounds = FALSE, include.intercept = FALSE, helpful.additional.cols = c("lPVL"), regions = region, times = time );
    # This one is good for showing that for 1m the relationship to lPVL is positive and for 6m it is negative.
#    formula.calibrated.withintercept.withlPVL.noBounds.none <- getTimingsResultFormula( "none", withbounds = FALSE, include.intercept = TRUE, helpful.additional.cols = c("lPVL"), regions = region, times = time );
    predicted.values.calibrated.COB <- getPredictions( region, bounded.COB.var, include.intercept = include.intercept, mutation.rate.calibration = mutation.rate.calibration, withbounds = withbounds, helpful.additional.cols = c( "lPVL" ), regions = train.regions, times = time, method = "glm" );
#    formula.calibrated.withintercept.withlPVL.withBounds.COB <- getTimingsResultFormula( bounded.COB.var, withbounds = TRUE, include.intercept = TRUE, helpful.additional.cols = c("lPVL"), regions = region, times = time );
#    formula.calibrated.withintercept.nolPVL.withBounds.COB <- getTimingsResultFormula( bounded.COB.var, withbounds = TRUE, include.intercept = TRUE, helpful.additional.cols = c(), regions = region, times = time );
    ## This one best demonstrates that the relationship with lPVL is negative for ~6m, positive for 1-2m:
 #   formula.calibrated.withintercept.withlPVL.noBounds.COB <- getTimingsResultFormula( bounded.COB.var, withbounds = FALSE, include.intercept = TRUE, helpful.additional.cols = c("lPVL"), regions = region, times = time );
    ## TODO: REMOVE
#    print( formula.calibrated.withintercept.withlPVL.noBounds.COB );
#    formula.calibrated.nointercept.withlPVL.withBounds.COB <- getTimingsResultFormula( bounded.COB.var, withbounds = TRUE, include.intercept = FALSE, helpful.additional.cols = c("lPVL"), regions = region, times = time );
    # For 1m results, we need to see these: 
#    formula.calibrated.nointercept.withlPVL.noBounds.COB <- getTimingsResultFormula( bounded.COB.var, withbounds = FALSE, include.intercept = FALSE, helpful.additional.cols = c("lPVL"), regions = region, times = time );
    predicted.values.calibrated.Prankenbeast <- getPredictions( region, bounded.PrankenBeast.var, include.intercept = include.intercept, mutation.rate.calibration = mutation.rate.calibration, withbounds = withbounds, helpful.additional.cols = c( "lPVL" ), regions = train.regions, times = time, method = "glm" );
#    formula.calibrated.withintercept.withlPVL.withBounds.Prankenbeast <- getTimingsResultFormula( bounded.PrankenBeast.var, withbounds = TRUE, include.intercept = TRUE, helpful.additional.cols = c("lPVL"), regions = region, times = time );
#    formula.calibrated.nointercept.withlPVL.withBounds.Prankenbeast <- getTimingsResultFormula( bounded.PrankenBeast.var, withbounds = TRUE, include.intercept = FALSE, helpful.additional.cols = c("lPVL"), regions = region, times = time );
    predicted.values.calibrated.PFitter <- getPredictions( region, "PFitter.mut.rate.coef", include.intercept = include.intercept, mutation.rate.calibration = mutation.rate.calibration, withbounds = withbounds, helpful.additional.cols = c( "lPVL" ), regions = train.regions, times = time, method = "glm" );

#     formula.calibrated.withintercept.withlPVL.withBounds.PFitter <- getTimingsResultFormula( "PFitter.mut.rate.coef", withbounds = TRUE, include.intercept = TRUE, helpful.additional.cols = c("lPVL"), regions = region, times = time );
#     formula.calibrated.nointercept.withlPVL.withBounds.PFitter <- getTimingsResultFormula( "PFitter.mut.rate.coef", withbounds = TRUE, include.intercept = FALSE, helpful.additional.cols = c("lPVL"), regions = region, times = time );
#     formula.calibrated.nointercept.withlPVL.noBounds.PFitter <- getTimingsResultFormula( "PFitter.mut.rate.coef", withbounds = FALSE, include.intercept = FALSE, helpful.additional.cols = c("lPVL"), regions = region, times = time );
#     ## For the concept of calibrating the mutation rate.
#     formula.calibrated.nointercept.nolPVL.noBounds.PFitter <- getTimingsResultFormula( "PFitter.mut.rate.coef", withbounds = FALSE, include.intercept = FALSE, helpful.additional.cols = c(), regions = region, times = time );
#     ## For the concept of calibrating the mutation rate with lPVL. 
#     formula.calibrated.mutation.rate.calibration.withlPVL.noBounds.PFitter <- getTimingsResultFormula( "PFitter.mut.rate.coef", withbounds = FALSE, mutation.rate.calibration = TRUE, helpful.additional.cols = c("lPVL"), regions = region, times = time );
#     formula.calibrated.mutation.rate.calibration.nolPVL.noBounds.PFitter <- getTimingsResultFormula( "PFitter.mut.rate.coef", withbounds = FALSE, mutation.rate.calibration = TRUE, helpful.additional.cols = c("lPVL"), regions = region, times = time );

    predicted.values.calibrated.Synonymous.PFitter <- getPredictions( region, "Synonymous.PFitter.mut.rate.coef", include.intercept = include.intercept, mutation.rate.calibration = mutation.rate.calibration, withbounds = withbounds, helpful.additional.cols = c( "lPVL" ), regions = train.regions, times = time, method = "glm" );
    predicted.values.calibrated.multifounder.PFitter <- getPredictions( region, "multifounder.PFitter.mut.rate.coef", include.intercept = include.intercept, mutation.rate.calibration = mutation.rate.calibration, withbounds = withbounds, helpful.additional.cols = c( "lPVL" ), regions = train.regions, times = time, method = "glm" );
    predicted.values.calibrated.multifounder.Synonymous.PFitter <- getPredictions( region, "multifounder.Synonymous.PFitter.mut.rate.coef", include.intercept = include.intercept, mutation.rate.calibration = mutation.rate.calibration, withbounds = withbounds, helpful.additional.cols = c( "lPVL" ), regions = train.regions, times = time, method = "glm" );
    
    .mat <- cbind( "Gold Standard" = gold.standard, "Center of Bounds" = predicted.values.calibrated.COB, "PrankenBeast" = predicted.values.calibrated.Prankenbeast, "PFitter" = predicted.values.calibrated.PFitter, "(syn) PFitter" = predicted.values.calibrated.Synonymous.PFitter, "(w/in c) PFitter" = predicted.values.calibrated.multifounder.PFitter, "(w/in c+syn) PFitter" = predicted.values.calibrated.multifounder.Synonymous.PFitter );
    .mat.fewer <- cbind( "Gold Standard" = gold.standard, "Center of Bounds" = predicted.values.calibrated.COB, "PrankenBeast" = predicted.values.calibrated.Prankenbeast, "PFitter" = predicted.values.calibrated.PFitter );
    
    ## Do this twice, once with days, once with prediction errors.
    .the.data.frame <- melt(.mat);
    .the.data.frame.fewer <- melt(.mat.fewer);
    colnames( .the.data.frame ) <- c( "Person", "Estimator", "Days" );
    colnames( .the.data.frame.fewer ) <- c( "Person", "Estimator", "Days" );
    .ggp.days.calibrated <- makeGGPFromDataFrame( .the.data.frame );
    .ggp.days.calibrated.fewer <- makeGGPFromDataFrame( .the.data.frame.fewer );

    ## TODO: ? These are always ylim 0 to 300 for now.
    pdf( paste( "Calibrated Estimators of ", .config.string, ".pdf" ) );
    print( addYLims( .ggp.days.calibrated, 300, 0 ) )
    dev.off()

    ## TODO: ? These are always ylim 0 to 300 for now.
    pdf( paste( "Calibrated Estimators of", .config.string, "Fewer.pdf" ) );
    print( addYLims( .ggp.days.calibrated.fewer, 300, 0 ) )
    dev.off()
    
    ## second time, using errors.
    .mat.prediction.errors <- .mat - .mat[ , 1 ];
    .mat.prediction.errors.fewer <- .mat.fewer - .mat.fewer[ , 1 ];
    
    .the.data.frame <- melt(.mat.prediction.errors);
    .the.data.frame.fewer <- melt(.mat.prediction.errors.fewer);
    colnames( .the.data.frame ) <- c( "Person", "Estimator", "Days" );
    colnames( .the.data.frame.fewer ) <- c( "Person", "Estimator", "Days" );
    
    .Biases <-
        apply( .mat.prediction.errors, 2, calculateBiasesFromPredictionErrors );
    .Biases.fewer <-
        apply( .mat.prediction.errors.fewer, 2, calculateBiasesFromPredictionErrors );

    .RMSEs <-
        apply( .mat.prediction.errors, 2, calculateRMSEsFromPredictionErrors );
    .RMSEs.fewer <-
        apply( .mat.prediction.errors.fewer, 2, calculateRMSEsFromPredictionErrors );
    
    .ggp.prediction.errors.calibrated <- makeGGPFromDataFrame( .the.data.frame, .Biases, .RMSEs, .Biases.and.RMSEs.yloc = .Biases.and.RMSEs.yloc );
    .ggp.prediction.errors.calibrated.fewer <- makeGGPFromDataFrame( .the.data.frame.fewer, .Biases.fewer, .RMSEs.fewer, .Biases.and.RMSEs.yloc = .Biases.and.RMSEs.yloc );
    
    pdf( paste( "Prediction Errors of Calibrated Estimators of ", .config.string, ".pdf" ) );
    print( addYLims( .ggp.prediction.errors.calibrated, ylims ) )
    dev.off()
    
    pdf( paste( "Prediction Errors of Calibrated Estimators of", .config.string, "Fewer.pdf" ) );
    print( addYLims( .ggp.prediction.errors.calibrated.fewer, ylims ) )
    dev.off()
    
    colnames( .mat ) <- gsub( "\\s", " ", colnames( .mat ) ); # Remove newlines before saving it.
    write.csv( .mat, paste( "Calibrated Estimators of", .config.string, "data.csv" ) );
    colnames( .mat.prediction.errors ) <- gsub( "\\s", " ", colnames( .mat.prediction.errors ) ); # Remove newlines before saving it.
    write.csv( .mat.prediction.errors, paste( "Prediction Errors of Calibrated Estimators of", .config.string, "data.csv" ) );
    
    ## Now do the same thing but for the uncalibrated estimates.
    predicted.values.calibrated.none <- getPredictions( region, "none", include.intercept = TRUE, withbounds = FALSE, helpful.additional.cols = c(), regions = train.regions, times = time, method = "glm" );
    predicted.values.uncalibrated.COB <- getPredictions( region, bounded.COB.var, include.intercept = FALSE, withbounds = FALSE, helpful.additional.cols = c(), regions = train.regions, times = time, method = "uncalibrated" );
    predicted.values.uncalibrated.Prankenbeast <- getPredictions( region, bounded.PrankenBeast.var, include.intercept = FALSE, withbounds = FALSE, helpful.additional.cols = c(), regions = train.regions, times = time, method = "uncalibrated" );
    predicted.values.uncalibrated.PFitter <- getPredictions( region, "PFitter.time.est", include.intercept = FALSE, withbounds = FALSE, helpful.additional.cols = c(), regions = train.regions, times = time, method = "uncalibrated" );
    predicted.values.uncalibrated.Synonymous.PFitter <- getPredictions( region, "Synonymous.PFitter.time.est", include.intercept = FALSE, withbounds = FALSE, helpful.additional.cols = c(), regions = train.regions, times = time, method = "uncalibrated" );
    predicted.values.uncalibrated.multifounder.PFitter <- getPredictions( region, "multifounder.PFitter.time.est", include.intercept = FALSE, withbounds = FALSE, helpful.additional.cols = c(), regions = train.regions, times = time, method = "uncalibrated" );
    predicted.values.uncalibrated.multifounder.Synonymous.PFitter <- getPredictions( region, "multifounder.Synonymous.PFitter.time.est", include.intercept = FALSE, withbounds = FALSE, helpful.additional.cols = c(), regions = train.regions, times = time, method = "uncalibrated" );
    .mat <- cbind( "Gold Standard" = gold.standard, "Center of Bounds" = predicted.values.uncalibrated.COB, "PrankenBeast" = predicted.values.uncalibrated.Prankenbeast, "PFitter" = predicted.values.uncalibrated.PFitter, "(syn) PFitter" = predicted.values.uncalibrated.Synonymous.PFitter, "(w/in c) PFitter" = predicted.values.uncalibrated.multifounder.PFitter, "(w/in c+syn) PFitter" = predicted.values.uncalibrated.multifounder.Synonymous.PFitter );
    .mat.fewer <- cbind( "Gold Standard" = gold.standard, "Center of Bounds" = predicted.values.uncalibrated.COB, "PrankenBeast" = predicted.values.uncalibrated.Prankenbeast, "PFitter" = predicted.values.uncalibrated.PFitter );
    
    ## Do this twice, once with days, once with prediction errors.
    .the.data.frame <- melt(.mat);
    .the.data.frame.fewer <- melt(.mat.fewer);
    colnames( .the.data.frame ) <- c( "Person", "Estimator", "Days" );
    colnames( .the.data.frame.fewer ) <- c( "Person", "Estimator", "Days" );
    .ggp.days.uncalibrated <- makeGGPFromDataFrame( .the.data.frame );
    .ggp.days.uncalibrated.fewer <- makeGGPFromDataFrame( .the.data.frame.fewer );
    
    
    ## TODO: ? These are always ylim 0 to 300 for now.
    pdf( paste( "Uncalibrated Estimators of ", .config.string, ".pdf" ) );
    print( addYLims( .ggp.days.uncalibrated, 300, 0 ) )
    dev.off()

    ## TODO: ? These are always ylim 0 to 300 for now.
    pdf( paste( "Uncalibrated Estimators of", .config.string, "Fewer.pdf" ) );
    print( addYLims( .ggp.days.uncalibrated.fewer, 300, 0 ) )
    dev.off()
    
    ## second time, using errors.
    .mat.prediction.errors <- .mat - .mat[ , 1 ];
    .mat.prediction.errors.fewer <- .mat.fewer - .mat.fewer[ , 1 ];
    
    .the.data.frame <- melt(.mat.prediction.errors);
    .the.data.frame.fewer <- melt(.mat.prediction.errors.fewer);
    colnames( .the.data.frame ) <- c( "Person", "Estimator", "Days" );
    colnames( .the.data.frame.fewer ) <- c( "Person", "Estimator", "Days" );
    
    .Biases <-
        apply( .mat.prediction.errors, 2, calculateBiasesFromPredictionErrors );
    .Biases.fewer <-
        apply( .mat.prediction.errors.fewer, 2, calculateBiasesFromPredictionErrors );
    
    .RMSEs <-
        apply( .mat.prediction.errors, 2, calculateRMSEsFromPredictionErrors );
    .RMSEs.fewer <-
        apply( .mat.prediction.errors.fewer, 2, calculateRMSEsFromPredictionErrors );
    
    .ggp.prediction.errors.uncalibrated <- makeGGPFromDataFrame( .the.data.frame, .Biases, .RMSEs, .Biases.and.RMSEs.yloc = .Biases.and.RMSEs.yloc );
    .ggp.prediction.errors.uncalibrated.fewer <- makeGGPFromDataFrame( .the.data.frame.fewer, .Biases.fewer, .RMSEs.fewer, .Biases.and.RMSEs.yloc = .Biases.and.RMSEs.yloc );
    
    pdf( paste( "Prediction Errors of Uncalibrated Estimators of ", .config.string, ".pdf" ) );
    print( addYLims( .ggp.prediction.errors.uncalibrated, ylims ) )
    dev.off()
    
    pdf( paste( "Prediction Errors of Uncalibrated Estimators of", .config.string, "Fewer.pdf" ) );
    print( addYLims( .ggp.prediction.errors.uncalibrated.fewer, ylims ) )
    dev.off()
    
    colnames( .mat ) <- gsub( "\\s", " ", colnames( .mat ) ); # Remove newlines before saving it.
    write.csv( .mat, paste( "Uncalibrated Estimators of", .config.string, "data.csv" ) );
    colnames( .mat.prediction.errors ) <- gsub( "\\s", " ", colnames( .mat.prediction.errors ) ); # Remove newlines before saving it.
    write.csv( .mat.prediction.errors, paste( "Prediction Errors of Uncalibrated Estimators of", .config.string, "data.csv" ) );

# Next, plot calibrated vs uncalibrated estimates. These are for PFitter. Below are for COB.
    predicted.values.calibrated <- getPredictions( region, "PFitter.mut.rate.coef", include.intercept = FALSE, helpful.additional.cols = c(), regions = train.regions, times = time, method = "glm" );
    predicted.values.calibrated.withlPVL <- getPredictions( region, "PFitter.mut.rate.coef", include.intercept = FALSE, helpful.additional.cols = c("lPVL"), regions = train.regions, times = time, method = "glm" );
    predicted.values.calibrated.withBounds <- getPredictions( region, "PFitter.mut.rate.coef", withbounds = TRUE, include.intercept = FALSE, helpful.additional.cols = c(), regions = train.regions, times = time, method = "glm" );
    predicted.values.calibrated.withlPVL.withBounds <- getPredictions( region, "PFitter.mut.rate.coef", withbounds = TRUE, include.intercept = FALSE, helpful.additional.cols = c("lPVL"), regions = train.regions, times = time, method = "glm" );
    .mat <- cbind( "Gold Standard" = gold.standard, PFitter = predicted.values.uncalibrated.PFitter, "Calibrated\n(Constant)" = predicted.values.calibrated, "Calibrated\nwith PVL" = predicted.values.calibrated.withlPVL, "Calibrated\nwith Bounds" = predicted.values.calibrated.withBounds, "Calibrated\nwith Bounds and PVL" = predicted.values.calibrated.withlPVL.withBounds );
    .mat.fewer <- cbind( "Gold Standard" = gold.standard, PFitter = predicted.values.uncalibrated.PFitter, "Calibrated\n(Constant)" = predicted.values.calibrated, "Calibrated\n(Person Specific)" = predicted.values.calibrated.withlPVL );
    
    .the.data.frame <- melt(.mat);
    .the.data.frame.fewer <- melt(.mat.fewer);
    colnames( .the.data.frame ) <- c( "Person", "Estimator", "Days" );
    colnames( .the.data.frame.fewer ) <- c( "Person", "Estimator", "Days" );
    .ggp.PFitter.days.calibrated <- makeGGPFromDataFrame( .the.data.frame );
    .ggp.PFitter.days.calibrated.fewer <- makeGGPFromDataFrame( .the.data.frame.fewer );
    
    colnames( .mat ) <- gsub( "\\s", " ", colnames( .mat ) ); # Remove newlines before saving it.
    write.csv( .mat, paste( "Calibrated PFitter Estimates of", .config.string, "boxplot data.csv" ) );
    
    ### Now same thing but for COB
    # Next, plot calibrated vs uncalibrated estimates. These are for COB.
    predicted.values.calibrated <- getPredictions( region, bounded.COB.var, include.intercept = FALSE, helpful.additional.cols = c(), regions = train.regions, times = time, method = "glm" );
    predicted.values.calibrated.withlPVL <- getPredictions( region, bounded.COB.var, include.intercept = FALSE, helpful.additional.cols = c("lPVL"), regions = train.regions, times = time, method = "glm" );
    predicted.values.calibrated.withBounds <- getPredictions( region, bounded.COB.var, withbounds = TRUE, include.intercept = FALSE, helpful.additional.cols = c(), regions = train.regions, times = time, method = "glm" );
    predicted.values.calibrated.withlPVL.withBounds <- getPredictions( region, bounded.COB.var, withbounds = TRUE, include.intercept = FALSE, helpful.additional.cols = c("lPVL"), regions = train.regions, times = time, method = "glm" );
    
    .mat <- cbind( "Gold Standard" = gold.standard, "Center of Bounds" = predicted.values.uncalibrated.COB, "Calibrated\n(Constant)" = predicted.values.calibrated, "Calibrated\nwith PVL" = predicted.values.calibrated.withlPVL, "Calibrated\nwith Bounds" = predicted.values.calibrated.withBounds, "Calibrated\nwith Bounds and PVL" = predicted.values.calibrated.withlPVL.withBounds );
    .mat.fewer <- cbind( "Gold Standard" = gold.standard, "Center of Bounds" = predicted.values.uncalibrated.COB, "Calibrated\n(Constant)" = predicted.values.calibrated, "Calibrated\n(Person Specific)" = predicted.values.calibrated.withlPVL );
    
    .the.data.frame <- melt(.mat);
    .the.data.frame.fewer <- melt(.mat.fewer);
    colnames( .the.data.frame ) <- c( "Person", "Estimator", "Days" );
    colnames( .the.data.frame.fewer ) <- c( "Person", "Estimator", "Days" );
    .ggp.COB.days.calibrated <- makeGGPFromDataFrame( .the.data.frame );
    .ggp.COB.days.calibrated.fewer <- makeGGPFromDataFrame( .the.data.frame.fewer );
    colnames( .mat ) <- gsub( "\\s", " ", colnames( .mat ) ); # Remove newlines before saving it.
    write.csv( .mat, paste( "Calibrated Center of Bounds Estimates of", .config.string, "boxplot data.csv" ) );
        
    pdf( paste( "Calibrated PFitter Estimates of", .config.string, "boxplot.pdf" ) )
    print( .ggp.PFitter.days.calibrated )
    dev.off()
    
    pdf( paste( "Calibrated Center of Bounds Estimates of", .config.string, "boxplot.pdf" ) )
    print( .ggp.COB.days.calibrated )
    dev.off()
    
    return( NULL );
} # createResultsPlotAndSpreadsheet (..)

# The default option of setting train.regions to NULL will be the same as the.region except for the "6m" time point, which will use both regions. To override this, set train.regions to NA to use the.region always, or to any other value to use a specified train.regions set for all results.
createAllResultsPlotsAndSpreadsheets <- function ( regions = c( "v3", "nflg" ), times = c( "1m", "6m" ), train.regions = NULL, withbounds = c( FALSE, TRUE ), include.intercept = c( FALSE, TRUE ) ) {
    #### NOTE 1w and 1m results ARE ESTIMATES FROM the same region only; for the 6m results below we are using estimates from both models.
    for( the.region in regions ) {
        for( the.time in times ) {
            for( the.withbounds in withbounds ) {
                for( the.include.intercept in include.intercept ) {
                    if( the.time == "6m" ) {
                        train.regions <- c( "nflg", "v3" );
                    } else {
                        train.regions <- the.region;
                    }
                    createResultsPlotAndSpreadsheet(
                                                    include.intercept = the.include.intercept,
                                                    withbounds = the.withbounds,
                                                    train.regions = train.regions,
                                                    region = the.region,
                                                    time = the.time
                                                   );
                } # End foreach include.intercept 
            } # End foreach withbounds
        } # End foreach the.time
    } # End foreach the.region
} # createAllResultsPlotsAndSpreadsheets (..)

prepare.ismultiple.mat <- function ( ismultiple.mat, only.is.lasso = ismultiple.only.is.lasso, exclude.continuous.predictors = ismultiple.exclude.continuous.predictors, include.training.codes = ismultiple.include.training.codes ) {
    has.nflg <- grep( "nflg\\.", colnames( ismultiple.mat ) );
    has.v3 <- grep( "v3\\.", colnames( ismultiple.mat ) );
    has.1m <- grep( "mtn003|1m\\.", colnames( ismultiple.mat ) );
    has.6m <- grep( "hvtn502|\\.6m", colnames( ismultiple.mat ) );
    is.glm <- grep( "glm", colnames( ismultiple.mat ) );
    is.lasso <- grep( "lasso", colnames( ismultiple.mat ) );
    
    new.colnames <- colnames( ismultiple.mat );
    for( col.i in 1:length( new.colnames ) ) {
        ## TODO: Put this information onto the plot by creating new rows.
        if( col.i %in% has.nflg ) {
            if( col.i %in% has.v3 ) {
                region.part <- "R";#"both regions";
            } else {
                region.part <- "r"; #"same region";
            }
        } else {
            stopifnot( col.i %in% has.v3 );
            region.part <- "r";#"same region";
        }

        if( col.i %in% has.1m ) {
            if( col.i %in% has.6m ) {
                time.part <- "T";#"both times";
            } else {
                time.part <- "t";#"same time";
            }
        } else {
            stopifnot( col.i %in% has.6m );
            time.part <- "t";#"same time";
        }
        covars.part <- "[bound]+";
        if( col.i %in% is.glm ) {
            intercept.part <- "1+";
            type.part <- paste( "single-founder~", intercept.part, covars.part, sep = "" );
        } else if( col.i %in% is.lasso ) {
            intercept.part <- "1+";
            type.part <- paste( "single-founder~", intercept.part, covars.part, "[lasso]+", sep = "" );
        } else {
            type.part <- "uncalibrated";
        }
        new.colname <- paste( paste( "[", region.part, time.part, "]", sep = "" ), type.part );
        #new.colname <- type.part;
        #print( new.colname );
        new.colnames[ col.i ] <- new.colname;
    }
    #colnames( ismultiple.mat ) <- gsub( "(nflg\\.)?(v3\\.)?(1m\\.)?(6m\\.)?", "", colnames( ismultiple.mat ) );
    colnames( ismultiple.mat ) <- new.colnames;

    ## The "uncalibrated" cols should all be identical; keep only the first.
    .cols <- grep( "uncalibrated", colnames( ismultiple.mat ) );
    if( length( .cols ) > 1 ) {
        ismultiple.mat.uncalibrated <- ismultiple.mat[ , .cols, drop = FALSE ];
        stopifnot( all( apply( ismultiple.mat.uncalibrated, 1, function( .row ) { all( .row == .row[1], na.rm = TRUE ); } ) ) );
        colnames( ismultiple.mat )[ .cols[ 1 ] ] <- "uncalibrated";
        ismultiple.mat <- ismultiple.mat[ , -.cols[ -1 ], drop = FALSE ]; # keep only the first.
    }
    
    cols.with.included.training.code <- unlist( lapply( include.training.codes, function( .training.code ) { grep( paste( "\\[", .training.code, "\\]", sep = "" ), colnames( ismultiple.mat ) ) } ) );

    if( only.is.lasso ) {
        cols.with.included.training.code <-
            intersect( cols.with.included.training.code, is.lasso );
    }
    
    # Filter down to the requested training codes. Always include the first column, which is uncalibrated.
    ismultiple.mat <- ismultiple.mat[ , c( 1, cols.with.included.training.code ), drop = FALSE ];
    
    if( exclude.continuous.predictors ) {
        # the continous predictors have no uncalibrated AUC value.
        ismultiple.mat <-
            ismultiple.mat[ !is.na( ismultiple.mat[ , "uncalibrated" ] ), , drop = FALSE ];
    }
    
    # Change the first one from "InSites single-founder" to "Rolland HVTN"
    rownames( ismultiple.mat )[ rownames( ismultiple.mat ) == "InSites single-founder" ] <-
        "Rolland HVTN";
    # Change the other ones, "(w/in clusts) (syn)" to "(w/in c+syn)" and "(w/in clusts)" to "(w/in c)"
    rownames( ismultiple.mat ) <- gsub( "\\(w/in clusts\\) \\(syn\\)", "(w/in c+syn)", rownames( ismultiple.mat ) );
    rownames( ismultiple.mat ) <- gsub( "\\(w/in clusts\\)", "(w/in c)", rownames( ismultiple.mat ) );
    return( ismultiple.mat );
} # prepare.ismultiple.mat ( .. )

plotIsMultipleMat <- function ( mat.data, do.barplots = !ismultiple.show.heatmap, start.barplot.from.half = FALSE ) {
    if( do.barplots ) {
        ### An alternative depiction, using side-by-side barplots.
        # Ok, use the AUCs as heights of bars.
        colnames( mat.data ) <- c( "Uncalibrated", "Calibrated" );
        mat.data.melted <- melt( mat.data );

        colnames( mat.data.melted ) <- c( "Predictor", "Calibration", "AUC" );
        mat.data.melted$Predictor <- factor( mat.data.melted$Predictor, levels = rownames( mat.data ) );
        mat.data.melted$Calibration <- relevel( mat.data.melted$Calibration, "Uncalibrated" );

        ggplot( mat.data.melted, aes( x = Predictor, y = AUC, fill = Calibration ) ) +
            geom_bar(stat="identity", color="black", position=position_dodge() ) +
            guides(fill=FALSE) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
      axis.text.x =
          element_text(
              angle = 90, vjust = 0.95, hjust=0.95,
              colour="black",size=16,face="bold"
           ),
      axis.text.y =
          element_text(
              colour="black",size=16,face="bold"
           ),
      axis.title.y =
          element_text(
              colour="black",size=16,face="bold"
           )
        ) + ylab( "Area Under the ROC Curve" ) + xlab( "" ) + scale_y_continuous( limits = c( 0, 1 ) ) + coord_cartesian( ylim = c( ifelse( start.barplot.from.half, 0.5, 0 ), 1 ) );
    } else {
        # do heatmaps.
        mat.text <- mat.data;

        my.palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                                       
        col.breaks = c(seq(0,0.5,length=100),  # for red
                       seq(0.5,0.9,length=101)[-1],              # for yellow
                       seq(.9,1,length=101)[-1])              # for green

        heatmap.2(
            mat.data,
            col = my.palette, colsep=c(1:62),rowsep=(1:62), sepwidth=c(0.05,0.05), sepcolor="white",
            trace="none", Rowv=F,Colv=F, scale="none", dendrogram="none",key=F,
            cellnote = mat.text, notecol = "black", notecex = 0.9,
            srtCol=45, adjCol = c(1,1), cexCol = 0.7,
            cexRow = 0.7,
            margins = c( 11.5, 9.5 ), breaks = col.breaks
        )
    }
} # plotIsMultipleMat ( mat.data, do.barplots )

getIsMultipleResultFormula <- function (
  varname = "none",
  regions = c( "nflg", "v3" ),
  times = c( "1m", "6m" ),
  helpful.additional.cols = c( "diversity", "priv.sites", "DSStarphyTest.R", "inf.sites.clusters", "lPVL" ),
  force.recomputation = FALSE,
  evaluate.regions = regions,
  evaluate.times = times,
  evaluate.varname = varname
) {

    stopifnot( varname %in% rownames( uses.by.evaluator.1m.6m.raw ) );
    
# Ensure results exist.
.results.by.region.and.time.Rda.filename <- evaluateIsMultiple(
     helpful.additional.cols = helpful.additional.cols,
     force.recomputation = force.recomputation,
     RESULTS.DIR = RESULTS.DIR,
     results.dirname = RESULTS.DIRNAME
);

    if( length( times ) == 2 ) { 
        the.bounds.type <- "sampledwidth_uniform_1mmtn003_6mhvtn502";
    } else if( times == "1m" ) {
        the.bounds.type <- "sampledwidth_uniform_mtn003";
    } else {
        stopifnot( times == "6m" );
        the.bounds.type <- "sampledwidth_uniform_hvtn502";
    }
    
#print( .results.by.region.and.time.Rda.filename );
local({
        load( .results.by.region.and.time.Rda.filename ); # adds results.by.region.and.time to environment.
        model.vars <- 
            setdiff( names( which( uses.by.evaluator.1m.6m.raw[ varname, ] >= 0.5 ) )[ -1 ], varname ); # remove (intercept), we next put it back as 1+
        .formula.txt <- paste( "is.one.founder ~ 1 + ", paste( model.vars, collapse = "+" ) );
        if( evaluate.varname != "none" ) {
            .formula.txt <- paste( .formula.txt, evaluate.varname, sep = "+" );
        }
        .formula <- 
            as.formula( .formula.txt );
        .evaluated.formula <-
            evaluate.specific.isMultiple.model.formula( results.by.region.and.time, .formula, step = FALSE, regions = evaluate.regions, times = evaluate.times );
        .evaluated.formula.table <- cbind( coef( summary( .evaluated.formula ) )[ , "Estimate", drop = FALSE ], confint( .evaluated.formula ), coef( summary( .evaluated.formula ) )[ , "Pr(>|z|)", drop = FALSE ] );
        return( list( formula = .formula, evaluated.formula = summary( .evaluated.formula ), evaluated.formula.table = .evaluated.formula.table ) );
     });
} # getIsMultipleResultFormula (..)


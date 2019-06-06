source( "evaluateIsMultiple_safetosource.R" )
source( "getResultsByRegionAndTime_safetosource.R" )
source( "evaluateIsMultiple_safetosource.R" )
source( "getFilteredResultsTables_safetosource.R" )

FORCE.RECOMPUTATION <- FALSE;

## Configure which results we are displaying.
BakeOff.RESULTS.DIR <- Sys.getenv( "BakeOff_RESULTS_DIR" );
if( BakeOff.RESULTS.DIR != "" ) {
    RESULTS.DIR <- BakeOff.RESULTS.DIR;
} else {
    #RESULTS.DIR <- "/fast/bakeoff_merged_analysis_sequences_unfiltered/results/";
    RESULTS.DIR <- "~/bakeoff_merged_analysis_sequences_filteredPre2017/results/";"/fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/";
    #RESULTS.DIR <- "/fast/bakeoff_merged_analysis_sequences_filtered2019/results/";
}
results.dirname <- "raw_fixed";

is.multiple.results.by.region.and.time.Rda.filename <-
    paste( RESULTS.DIR, results.dirname, "/isMultiple.results.by.region.and.time.Rda", sep = "" );

# load the results (must have been saved previously; see evaluateIsMultiple_safetosource.R)
        # loads results.by.region.and.time
load( file = is.multiple.results.by.region.and.time.Rda.filename );

calculateRSquaredValuesForIsMultiple <- function( results.by.region.and.time, step = FALSE, regions = c( "nflg", "v3" ), times = c( "1m", "6m", unbounded = FALSE ) ) {
    if( length( times == 2 ) ) {
        # good.
    } else {
        stop( "need to implement a little more to get the individual time-trained results here" );
    }
    ## TODO: Fix this, broken because the varnames are renamed:
#     if( length( regions == 2 ) ) {
#         uses.by.evaluator.1m.6m <- getFilteredLassoUsageTables( "isMultiple", the.region = "nflg", to.region = "v3", the.time = "1m.6m", RESULTS.DIR = RESULTS.DIR, results.dirname = RESULTS.DIRNAME );
#     } else {
#         stop( "need to implement a little more to get the individual region-trained results here" );
                                        #     }
                                        #    model.vars <- setdiff( colnames( uses.by.evaluator.1m.6m ), c( "Intercept", "(Intercept)", "the.intercept" ) );
    model.vars <- c( "lPVL", "bound", "time since diagnosis", "diversity", "inf.sites.clusters", "priv.sites", "DS.Starphy.R", "multifounder.Synonymous.DS.StarPhy.R" );
       
    the.bound <- "";
    if( unbounded ) {
        the.bound <- "unbounded";
    }
    if( length( times ) == 2 ) {
        time <- "1m.6m";
        if( the.bound != "unbounded" ) {
            the.bound <- "sampledwidth_uniform_1mmtn003_6mhvtn502";
        }
    } else {
        time <- times;
                if( time == "6m" ) {
                    if( the.bound != "unbounded" ) {
                        the.bound <- "sampledwidth_uniform_hvtn502";
                    }
                    stopifnot( length( times ) == 1 );
                    stopifnot( times == "6m" );
                } else {
                    if( the.bound != "unbounded" ) {
                        the.bound <- "sampledwidth_uniform_mtn003";
                    }
                    stopifnot( length( times ) == 1 );
                    stopifnot( times == "1m" );
                }
            }

.R.values <- 
    sapply( model.vars, function ( .var ) {
        print( .var )
        # replace "bound" with the appropriate var
        if( .var == "bound" ) {
            .var <- paste( the.bound, ".upper", sep = "" );
        }
        if( .var == "time since diagnosis" ) {
            .var <- paste( the.bound, ".lower", sep = "" );
        }
        if( .var == "unbounded" ) {
            return( NA );
        }
        if( length( grep( "DS.Starphy.R", .var ) ) > 0 ) {
            .var <- gsub( "DS.Starphy.R", "DSStarPhyTest.R", .var );
        }
        if( length( grep( "DS.StarPhy.R", .var ) ) > 0 ) {
            .var <- gsub( "DS.StarPhy.R", "DSStarPhyTest.R", .var );
        }
        .formula <- as.formula( paste( "is.one.founder ~ ", .var ) );
        .lm <- evaluate.specific.isMultiple.model.formula( results.by.region.and.time, .formula, step = step, regions = regions, times = times )
        compute.pearson.R.of.specific.isMultiple.model.formula.residuals.with.gold.standard( .lm )
    } )
    .R.squareds <- ( .R.values * .R.values );
    percent.of.variance.explained <- sprintf( "%0.2f%%", 100 * .R.squareds );
    names( percent.of.variance.explained ) <- names( .R.values );
    return( percent.of.variance.explained );
} # calculateRSquaredValuesForIsMultiple (..)

calculateRSquaredValuesForIsMultiple( results.by.region.and.time )

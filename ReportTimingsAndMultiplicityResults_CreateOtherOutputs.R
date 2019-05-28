source( "ReportTimingsAndMultiplicityResults_safetosource.R" );

#############
# CONFIGURATION 
############
timings.include.methods <- c( "glm" ); # c( "glm", "step", "lasso" )
timings.include.training.codes <- c( "RT", "Rt", "rT", "rt" ); 

ismultiple.show.heatmap = FALSE;
if( ismultiple.show.heatmap ) {
    ismultiple.include.training.codes <- c( "RT", "Rt", "rT", "rt" );
    ismultiple.only.is.lasso <- FALSE;
    ismultiple.exclude.continuous.predictors <- FALSE;
} else {
    ismultiple.include.training.codes <-  "RT";
    ismultiple.only.is.lasso <- TRUE;
    ismultiple.exclude.continuous.predictors <- TRUE;
}

###########
## DEFAULTS
###########

## Configure which results we are displaying.
BakeOff.RESULTS.DIR <- Sys.getenv( "BakeOff_RESULTS_DIR" );
if( BakeOff.RESULTS.DIR != "" ) {
    RESULTS.DIR <- BakeOff.RESULTS.DIR;
} else {
    #RESULTS.DIR <- "/fast/bakeoff_merged_analysis_sequences_unfiltered/results/";
    RESULTS.DIR <- "/fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/";
    #RESULTS.DIR <- "/fast/bakeoff_merged_analysis_sequences_filtered2019/results/";
}
RESULTS.DIRNAME <- "raw_fixed";

## NOTE even setting this to TRUE will not force all to be recomputed. See below where force.recomputation = FALSE.
FORCE.RECOMPUTATION <- FALSE;#TRUE;

## Ensure results exist and create the all of the output files.
createAllResultsPlotsAndSpreadsheets();


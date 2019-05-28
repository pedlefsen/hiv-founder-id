## First do all the stuff in README.preprocessing.txt.

source( "getDaysSinceInfection_safetosource.R" );

## Configure which results we are displaying.
BakeOff.RESULTS.DIR <- Sys.getenv( "BakeOff_RESULTS_DIR" );
if( BakeOff.RESULTS.DIR != "" ) {
     RESULTS.DIR <- BakeOff.RESULTS.DIR;
} else {
    #RESULTS.DIR <- "/fast/bakeoff_merged_analysis_sequences_unfiltered/results/";
    RESULTS.DIR <- "/fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/";
    #RESULTS.DIR <- "/fast/bakeoff_merged_analysis_sequences_filtered2019/results/";
}

results.dirname <- "raw_fixed";
#results.dirname <- "raw";

bounds.subdirname <- "bounds";

# force.recomputation <- FALSE;

rv217.gold.standard.infection.dates.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/rv217_gold_standard_timings.csv" );
rv217.gold.standard.infection.dates <- as.Date( as.character( rv217.gold.standard.infection.dates.in[,2] ), "%m/%d/%y" );
names( rv217.gold.standard.infection.dates ) <- as.character( rv217.gold.standard.infection.dates.in[,1] );

caprisa002.gold.standard.infection.dates.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/caprisa_002_gold_standard_timings.csv" );
caprisa002.gold.standard.infection.dates <- as.Date( as.character( caprisa002.gold.standard.infection.dates.in[,2] ), "%Y/%m/%d" );
names( caprisa002.gold.standard.infection.dates ) <- as.character( caprisa002.gold.standard.infection.dates.in[,1] );

# These are from HVTN 502 / HVTN 504 (STEP), which had a scheduled visit every 6 months (I believe) during the period after the vaccination phase, but these values I believe include everyone (those infected during vaccination, during the main trial, and even during unblinded followup (HVTN 504). I excluded one very high outlier value (1496).  The range of the remaining are from 18 to 495.
hvtn502.timing.windows.of.infecteds <- read.table( file = paste( RESULTS.DIR, results.dirname, "/", bounds.subdirname, "/infectionWindowInDays_v502.csv", sep = "" ), header = TRUE )[[1]];

# These are from HVTN 503 (Phambili), which had the same visit schedule as HVTN 502, but was halted early. I've excluded the two zero-valued days_negminpos and the two high outliers: 742 and 1008.  The rest are in the reasonable 1-month to 1-year range.
hvtn503.timing.windows.of.infecteds <- read.table( file = paste( RESULTS.DIR, results.dirname, "/", bounds.subdirname, "/infectionWindowInDays_v503.csv", sep = "" ), header = TRUE )[[1]];

# These are from HVTN 505, which during this period (blinded-phase of study, not unblinded followup) had a scheduled visit every 3 months, though it may have included more frequent visits during the vaccination phase (I do not know).
hvtn505.timing.windows.of.infecteds <- read.table( file = paste( RESULTS.DIR, results.dirname, "/", bounds.subdirname, "/infectionWindowInDays_v505.csv", sep = "" ), header = TRUE )[[1]];

# These are from MTN 003, which had a 4-week visit window but this have a long large tail.
mtn003.timing.windows.of.infecteds <- read.table( file = paste( RESULTS.DIR, results.dirname, "/", bounds.subdirname, "/infectionWindowInDays_m003.csv", sep = "" ), header = TRUE )[[1]];

## Some exploration suggests that the 503 data are quite different from the 502 data.
## Also, the distribution of the 502 and 505 data is pretty similar if you scale by the ratio of means (which is nearly 2 = 6 months / 3 months, the ratio of the window widths).
# qqplot( hvtn502.timing.windows.of.infecteds, ( mean(hvtn502.timing.windows.of.infecteds)/mean(hvtn505.timing.windows.of.infecteds) ) *hvtn505.timing.windows.of.infecteds )
# boxplot( hvtn502.timing.windows.of.infecteds, ( mean(hvtn502.timing.windows.of.infecteds)/mean(hvtn505.timing.windows.of.infecteds) ) *hvtn505.timing.windows.of.infecteds )
## There is a heavy upper tail. The data follow something like a mixture of a normal with range 0 to 1.5*windowsize and a beta over that range with mode at the target time (90 days or 180 days) and a third component describing the heavy right tail. -- but the beta part is not necessary see below qqnorm plot (pooling the data after scaling).
# qqplot( hvtn502.timing.windows.of.infecteds[ hvtn502.timing.windows.of.infecteds < 270 ], ( mean(hvtn502.timing.windows.of.infecteds)/mean(hvtn505.timing.windows.of.infecteds) ) *hvtn505.timing.windows.of.infecteds[ hvtn502.timing.windows.of.infecteds < 130 ] )
# qqnorm( c( hvtn502.timing.windows.of.infecteds[ hvtn502.timing.windows.of.infecteds < 270 ], ( mean(hvtn502.timing.windows.of.infecteds)/mean(hvtn505.timing.windows.of.infecteds) ) *hvtn505.timing.windows.of.infecteds[ hvtn502.timing.windows.of.infecteds < 130 ] ) )
.scaled.pooled.widths.excludingrighttail <- c( hvtn502.timing.windows.of.infecteds[ hvtn502.timing.windows.of.infecteds < 270 ], ( mean(hvtn502.timing.windows.of.infecteds)/mean(hvtn505.timing.windows.of.infecteds) ) *hvtn505.timing.windows.of.infecteds[ hvtn502.timing.windows.of.infecteds < 130 ] );
.scaled.pooled.widths <- c( hvtn502.timing.windows.of.infecteds, ( mean(hvtn502.timing.windows.of.infecteds)/mean(hvtn505.timing.windows.of.infecteds) ) *hvtn505.timing.windows.of.infecteds );

# What would the mean be if these were 30-day windows instead of 90 or 180?  What is the relationship bn target and mean?
mean(hvtn505.timing.windows.of.infecteds) / 90
# [1] 0.9445988
mean(hvtn502.timing.windows.of.infecteds) / 180
# [1] 0.859767
180 - mean(hvtn502.timing.windows.of.infecteds)
# [1] 25.24194
90 - mean(hvtn505.timing.windows.of.infecteds)
# [1] 4.986111

# Basically if you double the window length you make the difference between the mean and the target greater by 5-fold.
x = c( 90, 180 ); y = c( 90-mean(hvtn505.timing.windows.of.infecteds), 180-mean(hvtn502.timing.windows.of.infecteds) )
y[2]/y[1]
# [1] 5.062449
x[2]/x[1]
# [1] 2
the.slope.logpace <- ( log( y[2]/y[1] ) ) / ( log( x[2]/x[1] ) );
new.y.from.new.x <- function ( new.x ) {
    new.run <- log(new.x/x[1]);
    new.rise <- new.run * the.slope.logpace;
    exp( new.rise ) * y[1];
}
stopifnot( abs( new.y.from.new.x( x[1] ) - y[1] ) == 0 );
stopifnot( abs( new.y.from.new.x( x[2] ) - y[2] ) == 0 );

## So for the 1-month data:
x.onemonth <- 30;
y.onemonth <- new.y.from.new.x( x.onemonth );
# y.onemonth is the diff between the mean of the distribution and the target, so here it is.
mean.onemonth <- x.onemonth - y.onemonth;
sd.onemonth <- ( mean.onemonth / mean( .scaled.pooled.widths.excludingrighttail ) ) * sd( .scaled.pooled.widths.excludingrighttail );
# Create the fake 1-month data by sampling from the actual pooled widths, rescaled.
onemonth.rescaled.pooled.widths <- .scaled.pooled.widths * ( mean.onemonth / mean( .scaled.pooled.widths ) );
## Now we sample from these for the 1-month fake boundaries.
## For the 6-month fake boundaries we use this, which basically just pools in the scaled 505 (3-month) windows into the 502 6-month windows.
x.sixmonths <- 180;
y.sixmonths <- new.y.from.new.x( x.sixmonths );
# y.sixmonths is the diff between the mean of the distribution and the target, so here it is.
mean.sixmonths <- x.sixmonths - y.sixmonths;
sd.sixmonths <- ( mean.sixmonths / mean( .scaled.pooled.widths.excludingrighttail ) ) * sd( .scaled.pooled.widths.excludingrighttail );
# Create the fake 1-month data by sampling from the actual pooled widths, rescaled.
sixmonths.rescaled.pooled.widths <- .scaled.pooled.widths * ( mean.sixmonths / mean( .scaled.pooled.widths ) );

regions <- c( "nflg", "v3" );

# Eg interval.center = 0, interval.width = 35, interval.center.percentile = 0.5
create.deterministic.interval.generation.fn <- function( interval.center.percentile ) {
    #print( interval.center.percentile );
    stopifnot( is.numeric( interval.center.percentile ) && ( 0 <= interval.center.percentile ) && ( 1 >= interval.center.percentile ) );

    function ( interval.center, interval.width ) {
        #print( interval.center );
        #print( interval.width );
      interval.width.below.center <- round( interval.width * interval.center.percentile );
      interval.width.above.center <- round( interval.width - interval.width.below.center );
        #print( interval.width.below.center );
        #print( interval.width.above.center );
      return( c( interval.center - interval.width.below.center, interval.center + interval.width.above.center ) );
    }
} # create.deterministic.interval.generation.fn ( .. )

# This generates an interval of a fixed width but randomly places it; simply calls create.deterministic.interval.generation.fn with a uniformly-chosen percentile.
uniform.interval.generation.fn <- function ( interval.center, interval.width ) {
    return( create.deterministic.interval.generation.fn( runif( 1 ) )( interval.center, interval.width ) );
} # uniform.interval.generation.fn ( .. )

# This generates an interval of a width sampled from the given set, and then randomly places it; simply calls uniform.interval.generation.fn with a randomly-chosen width.
sampledwidth.uniform.interval.generation.fn <- function ( interval.center, interval.widths.to.sample.from ) {
    return( create.deterministic.interval.generation.fn( runif( 1 ) )( interval.center, sample( interval.widths.to.sample.from, 1 ) ) );
} # sampledwidth.uniform.interval.generation.fn ( .. )

# This generates an interval of random (exponentially-distributed) width and then randomly places it; simply calls create.deterministic.interval.generation.fn with a uniformly-chosen percentile and an exponentially-distributed width.
exponential.uniform.interval.generation.fn <- function ( interval.center, mean.interval.width ) {
    return( create.deterministic.interval.generation.fn( runif( 1 ) )( interval.center, rexp( 1, 1 / mean.interval.width ) ) );
} # exponential.uniform.interval.generation.fn ( .. )

# This generates an interval of random (gamma-distributed) width and then randomly places it; simply calls create.deterministic.interval.generation.fn with a uniformly-chosen percentile and a gamma-distributed width.
gamma.uniform.interval.generation.fn <- function ( interval.center, mean.interval.width, sd.interval.width = (1/5)*mean.interval.width ) {
    # mean of a gamma is .shape * .scale.
    .shape <- 1 / ( sd.interval.width / mean.interval.width )^2;
    .scale <- mean.interval.width/.shape;
    
    return( create.deterministic.interval.generation.fn( runif( 1 ) )( interval.center, rgamma( 1, shape = .shape, rate = 1/.scale ) ) );
} # gamma.uniform.interval.generation.fn ( .. )

test.uniform.interval.generation.fn <- function () {
    stopifnot( all.equal( mean( replicate( 10000, uniform.interval.generation.fn( 0, 100 )[ 1 ] ) ), -50, tol = 1E-2 ) );
    stopifnot( all.equal( mean( replicate( 100000, uniform.interval.generation.fn( 100, 1000 )[ 2 ] ) ), 600, tol = 1 ) );
} # test.uniform.interval.generation.fn ()

test.exponential.uniform.interval.generation.fn <- function () {
    stopifnot( all.equal( mean( apply( replicate( 10000, exponential.uniform.interval.generation.fn( 100, 5*7 ) ), 2, diff ) ), ( 5 * 7 ), tol = 1 ) )
} # test.exponential.uniform.interval.generation.fn ()

test.gamma.uniform.interval.generation.fn <- function () {
    stopifnot( all.equal( mean( apply( replicate( 10000, gamma.uniform.interval.generation.fn( 100, 5*7, 2*7 ) ), 2, diff ) ), ( 5 * 7 ), tol = 1 ) )
    stopifnot( all.equal( sd( apply( replicate( 10000, gamma.uniform.interval.generation.fn( 100, 5*7, 2*7 ) ), 2, diff ) ), ( 2 * 7 ), tol = 1 ) )
} # test.gamma.uniform.interval.generation.fn ()

## For now "interval.center.percentile" should be a fraction, so eg 0.5 means that the interval will be centered on the true infection date, 0.1 means that the 10th percentile of the interval will be placed at the true value, and 1.0 means that the interval will end at the true infection date.
## By default it'll be deterministic, but if you change interval.generation.fn to a function, it should accept the following parameters: [TODO]
## Note that it is up to the CALLER to set the random seed if the fn is random.
createArtificialBoundsOnInfectionDate <-
    function ( interval.width.in.days = ( 5 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.5 ), output.file.suffix = "deterministic_5weeks_centered.tab", times = c( "1m", "6m", "1m6m" ) )
{
    .result.ignored <-
        lapply( regions, function( the.region ) {
            ## TODO: REMOVE
            cat( the.region, fill = T );
            ..result.ignored <- 
                lapply( times, function( the.time ) {
                    ## TODO: REMOVE
                    cat( the.time, fill = T );

                    .days.since.infection.filename <-
                        paste( RESULTS.DIR, results.dirname, "/", the.region, "/", the.time, "/sampleDates.tbl", sep = "" );
                    if( the.region == "v3" ) {
                        days.since.infection <-
                            getDaysSinceInfection(
                                .days.since.infection.filename,
                                caprisa002.gold.standard.infection.dates
                            );
                    } else {
                        stopifnot( ( the.region == "nflg" ) || ( length( grep( "rv217", the.region ) ) > 0 ) );
                        days.since.infection <-
                            getDaysSinceInfection(
                                .days.since.infection.filename,
                                rv217.gold.standard.infection.dates
                            );
                    }

                    the.interval.by.ppt <-
                        t( sapply( days.since.infection, interval.generation.fn, interval.width.in.days ) );
                    colnames( the.interval.by.ppt ) <- c( "lower", "upper" );
                    
                    ## Write it out as we go.
                    .artificial.bounds.dirname <-
                        paste( RESULTS.DIR, results.dirname, "/", bounds.subdirname, "/", the.region, "/", the.time, "/", sep = "" );
                    dir.create( .artificial.bounds.dirname, recursive = TRUE, showWarnings = FALSE );
                    .artificial.bounds.filename <-
                        paste( .artificial.bounds.dirname, "artificialBounds_", output.file.suffix, sep = "" );

                    write.table( the.interval.by.ppt, .artificial.bounds.filename, quote = FALSE, sep = "\t" );
                    
                    return( NULL );
                } ); # End foreach the.time

            return( NULL );
        } ); # End foreach the.region
    
    return( NULL );
} # createArtificialBoundsOnInfectionDate (..)

## DO it.

## Uniform-center, sampled from onemonth.rescaled.pooled.widths.
# set.seed( 98103 );
# createArtificialBoundsOnInfectionDate( interval.width.in.days = onemonth.rescaled.pooled.widths, interval.generation.fn = sampledwidth.uniform.interval.generation.fn, output.file.suffix = "sampledwidth_uniform_onemonth.tab", times = c( "1m", "1m6m" ) );
# 1w was done separately:
set.seed( 98103 );
createArtificialBoundsOnInfectionDate( interval.width.in.days = onemonth.rescaled.pooled.widths, interval.generation.fn = sampledwidth.uniform.interval.generation.fn, output.file.suffix = "sampledwidth_uniform_onemonth.tab", times = "1w" )

## Uniform-center, sampled from mtn003.timing.windows.of.infecteds (one-monthly, nominally)
set.seed( 98103 );
createArtificialBoundsOnInfectionDate( interval.width.in.days = mtn003.timing.windows.of.infecteds, interval.generation.fn = sampledwidth.uniform.interval.generation.fn, output.file.suffix = "sampledwidth_uniform_mtn003.tab", times = c( "1m", "1m6m" ) );
# 1w was done separately:
set.seed( 98103 );
createArtificialBoundsOnInfectionDate( interval.width.in.days = mtn003.timing.windows.of.infecteds, interval.generation.fn = sampledwidth.uniform.interval.generation.fn, output.file.suffix = "sampledwidth_uniform_mtn003.tab", times = "1w" );

## Uniform-center, sampled from sixmonth.rescaled.pooled.widths.
# set.seed( 98103 );
# createArtificialBoundsOnInfectionDate( interval.width.in.days = sixmonths.rescaled.pooled.widths, interval.generation.fn = sampledwidth.uniform.interval.generation.fn, output.file.suffix = "sampledwidth_uniform_sixmonths.tab", times = c( "6m" ) );

## Uniform-center, sampled from hvtn502.timing.windows.of.infecteds (sixmonthly, nominally)
set.seed( 98103 );
createArtificialBoundsOnInfectionDate( interval.width.in.days = hvtn502.timing.windows.of.infecteds, interval.generation.fn = sampledwidth.uniform.interval.generation.fn, output.file.suffix = "sampledwidth_uniform_hvtn502.tab", times = c( "6m" ) );

## ARCHIVE:
# ## Deterministic, 5 weeks, centered.
# createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 5 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.5 ), output.file.suffix = "deterministic_5weeks_centered.tab" );
# ## Deterministic, 5 weeks, 10th percentile.
# createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 5 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.1 ), output.file.suffix = "deterministic_5weeks_percentile10.tab" );
# ## Deterministic, 5 weeks, 90th percentile.
# createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 5 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.9 ), output.file.suffix = "deterministic_5weeks_percentile90.tab" );

# ## Exponential-width Uniform-center, 5 weeks.
# set.seed( 98103 );
# createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 5 * 7 ), interval.generation.fn = exponential.uniform.interval.generation.fn, output.file.suffix = "exponentialwidth_uniform_5weeks.tab", times = c( "1m", "1m6m" ) );

# ## Gamma-width Uniform-center, mean 5 weeks, SD 1 week.
# set.seed( 98103 );
# createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 5 * 7 ), interval.generation.fn = gamma.uniform.interval.generation.fn, output.file.suffix = "gammawidth_uniform_5weeks.tab", times = c( "1m", "1m6m" ) );

# ## Uniform-center, 5 weeks.
# set.seed( 98103 );
# createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 5 * 7 ), interval.generation.fn = exponential.uniform.interval.generation.fn, output.file.suffix = "uniform_5weeks.tab", times = c( "6m" ) );

# ## Deterministic, 30 weeks, centered.
# createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 30 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.50 ), output.file.suffix = "deterministic_30weeks_centered.tab", times = c( "6m" ) );
# ## Deterministic, 30 weeks, 10th percentile.
# createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 30 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.1 ), output.file.suffix = "deterministic_30weeks_percentile10.tab", times = c( "6m" ) );
# ## Deterministic, 30 weeks, 90th percentile.
# createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 30 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.9 ), output.file.suffix = "deterministic_30weeks_percentile90.tab", times = c( "6m" ) );

# ## Exponential-width Uniform-center, 30 weeks.
# set.seed( 98103 );
# createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 30 * 7 ), interval.generation.fn = exponential.uniform.interval.generation.fn, output.file.suffix = "exponentialwidth_uniform_30weeks.tab", times = c( "6m" ) );

## Gamma-width Uniform-center, mean 30 weeks, SD 6 weeks.
# set.seed( 98103 );
# createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 30 * 7 ), interval.generation.fn = gamma.uniform.interval.generation.fn, output.file.suffix = "gammawidth_uniform_30weeks.tab", times = c( "6m" ) );

## Uniform-center, 30 weeks.
# set.seed( 98103 );
# createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 30 * 7 ), interval.generation.fn = uniform.interval.generation.fn, output.file.suffix = "uniform_30weeks.tab", times = c( "6m" ) );



## First do all the stuff in README.preprocessing.txt.

source( "getDaysSinceInfection_safetosource.R" );

results.dirname <- "raw_edited_20160216";
#results.dirname <- "raw";
bounds.subdirname <- "bounds";

# force.recomputation <- FALSE;

rv217.gold.standard.infection.dates.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/rv217_gold_standard_timings.csv" );
rv217.gold.standard.infection.dates <- as.Date( as.character( rv217.gold.standard.infection.dates.in[,2] ), "%m/%d/%y" );
names( rv217.gold.standard.infection.dates ) <- as.character( rv217.gold.standard.infection.dates.in[,1] );

caprisa002.gold.standard.infection.dates.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/caprisa_002_gold_standard_timings.csv" );
caprisa002.gold.standard.infection.dates <- as.Date( as.character( caprisa002.gold.standard.infection.dates.in[,2] ), "%Y/%m/%d" );
names( caprisa002.gold.standard.infection.dates ) <- as.character( caprisa002.gold.standard.infection.dates.in[,1] );

regions <- c( "nflg", "v3", "rv217_v3" );
times <- c( "1m", "6m", "1m6m" );

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

# This generates an interval of random (exponentially-distributed) width and then randomly places it; simply calls create.deterministic.interval.generation.fn with a uniformly-chosen percentile and an exponentially-distributed width.
exponential.uniform.interval.generation.fn <- function ( interval.center, mean.interval.width ) {
    return( create.deterministic.interval.generation.fn( runif( 1 ) )( interval.center, rexp( 1, 1 / mean.interval.width ) ) );
} # exponential.uniform.interval.generation.fn ( .. )

test.uniform.interval.generation.fn <- function () {
    stopifnot( all.equal( mean( replicate( 10000, uniform.interval.generation.fn( 0, 100 )[ 1 ] ) ), -50, tol = 1E-2 ) );
    stopifnot( all.equal( mean( replicate( 100000, uniform.interval.generation.fn( 100, 1000 )[ 2 ] ) ), 600, tol = 1 ) );
} # test.uniform.interval.generation.fn ()

test.exponential.uniform.interval.generation.fn <- function () {
    stopifnot( all.equal( mean( apply( replicate( 10000, exponential.uniform.interval.generation.fn( 100, 5*7 ) ), 2, diff ) ), ( 5 * 7 ), tol = 1 ) )
} # test.exponential.uniform.interval.generation.fn ()

## For now "interval.center.percentile" should be a fraction, so eg 0.5 means that the interval will be centered on the true infection date, 0.1 means that the 10th percentile of the interval will be placed at the true value, and 1.0 means that the interval will end at the true infection date.
## By default it'll be deterministic, but if you change interval.generation.fn to a function, it should accept the following parameters: [TODO]
## Note that it is up to the CALLER to set the random seed if the fn is random.
createArtificialBoundsOnInfectionDate <-
    function ( interval.width.in.days = ( 5 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.5 ), output.file.suffix = "deterministic_5weeks_centered.tab" )
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
                        paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/sampleDates.tbl", sep = "" );
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
                        paste( "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/", results.dirname, "/", bounds.subdirname, "/", the.region, "/", the.time, "/", sep = "" );
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

## Deterministic, 5 weeks, centered.
createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 5 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.5 ), output.file.suffix = "deterministic_5weeks_centered.tab" );
## Deterministic, 5 weeks, 10th percentile.
createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 5 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.1 ), output.file.suffix = "deterministic_5weeks_percentile10.tab" );
## Deterministic, 5 weeks, 90th percentile.
createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 5 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.9 ), output.file.suffix = "deterministic_5weeks_percentile90.tab" );

## Exponential-width Uniform-center, 5 weeks.
set.seed( 98103 );
createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 5 * 7 ), interval.generation.fn = exponential.uniform.interval.generation.fn, output.file.suffix = "exponentialwidth_uniform_5weeks.tab" );

## Uniform-center, 5 weeks.
set.seed( 98103 );
createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 5 * 7 ), interval.generation.fn = exponential.uniform.interval.generation.fn, output.file.suffix = "uniform_5weeks.tab" );

## Deterministic, 20 weeks, centered.
createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 20 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.50 ), output.file.suffix = "deterministic_20weeks_centered.tab" );
## Deterministic, 20 weeks, 10th percentile.
createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 20 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.1 ), output.file.suffix = "deterministic_20weeks_percentile10.tab" );
## Deterministic, 20 weeks, 90th percentile.
createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 20 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.9 ), output.file.suffix = "deterministic_20weeks_percentile90.tab" );

## Exponential-width Uniform-center, 20 weeks.
set.seed( 98103 );
createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 20 * 7 ), interval.generation.fn = exponential.uniform.interval.generation.fn, output.file.suffix = "exponentialwidth_uniform_20weeks.tab" );

## Uniform-center, 20 weeks.
set.seed( 98103 );
createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 20 * 7 ), interval.generation.fn = uniform.interval.generation.fn, output.file.suffix = "uniform_20weeks.tab" );

## Deterministic, 30 weeks, centered.
createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 30 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.50 ), output.file.suffix = "deterministic_30weeks_centered.tab" );
## Deterministic, 30 weeks, 10th percentile.
createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 30 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.1 ), output.file.suffix = "deterministic_30weeks_percentile10.tab" );
## Deterministic, 30 weeks, 90th percentile.
createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 30 * 7 ), interval.generation.fn = create.deterministic.interval.generation.fn( 0.9 ), output.file.suffix = "deterministic_30weeks_percentile90.tab" );

## Exponential-width Uniform-center, 30 weeks.
set.seed( 98103 );
createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 30 * 7 ), interval.generation.fn = exponential.uniform.interval.generation.fn, output.file.suffix = "exponentialwidth_uniform_30weeks.tab" );

## Uniform-center, 30 weeks.
set.seed( 98103 );
createArtificialBoundsOnInfectionDate( interval.width.in.days = ( 30 * 7 ), interval.generation.fn = uniform.interval.generation.fn, output.file.suffix = "uniform_30weeks.tab" );



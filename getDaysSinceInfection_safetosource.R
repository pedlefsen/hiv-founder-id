getDaysSinceInfection <- function ( sample.dates.tbl.filename, gold.standard.infection.dates, return.sample.dates.in = FALSE ) {
    sample.dates.in <- read.delim( sample.dates.tbl.filename, sep = " ", header = F, fill = T );
    colnames( sample.dates.in ) <- c( "ptid", "date" );
    ## Remove anything that's not really a ptid/date combo
    sample.dates.in <-
        sample.dates.in[ grep( "^\\d+$", as.character( sample.dates.in[ , 1 ] ) ), , drop = FALSE ];
    # remove anything with a missing date.
    sample.dates.in <-
        sample.dates.in[ sample.dates.in[ , 2 ] != "", , drop = FALSE ];
    if( return.sample.dates.in ) {
        return( sample.dates.in );
    }
    
    days.since.infection <-
        sapply( 1:nrow( sample.dates.in ), function( .i ) {
            as.numeric( as.Date( as.character( sample.dates.in[ .i, 2 ] ) ) - gold.standard.infection.dates[ as.character( sample.dates.in[ .i, 1 ] ) ] )
        } );
    names( days.since.infection ) <- sample.dates.in[ , "ptid" ];

    return( days.since.infection );
} # getDaysSinceInfection ( .. )

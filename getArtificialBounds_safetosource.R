getArtificialBounds <- function (
    the.region,
    the.time,
    results.dirname,
    bounds.subdirname = "bounds"
) {
    ## Ideally we'd use bounds
    ## on the actual infection time computed from the
    ## dates and results of the HIV positivity tests
    ## (antibody or PCR).  The typical approach is to
    ## perform antibody testing every X days
    ## (historically this is 6 months in most HIV
    ## vaccine trials, except during the vaccination
    ## phase there are more frequent visits and on
    ## every visit HIV testing is conducted).  The
    ## (fake) bounds used here are calculated in the
    ## createArtificialBoundsOnInfectionDate.R file.
    ## The actual bounds would be too tight, since the
    ## participants were detected HIV+ earlier in these
    ## people than what we expect to see in a trial in
    ## which testing is conducted every X days.  For
    ## the center of bounds approach we load the bounds
    ## files in subdirs of the "bounds" subdirectory eg
    ## at
    ## /fh/fast/edlefsen_p/bakeoff/analysis_sequences/bounds/nflg/1m/.
    ## These files have names beginning with
    ## "artificialBounds_" and ending with ".tab".
    .artificial.bounds.dirname <-
        paste( "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/", results.dirname, "/", bounds.subdirname, "/", the.region, "/", the.time, "/", sep = "" );
    artificial.bounds.filenames <-
        dir( .artificial.bounds.dirname, pattern = "artificialBounds_.*.tab", recursive = FALSE, full.names = TRUE );
    names( artificial.bounds.filenames ) <- gsub( "^.*artificialBounds_(.*).tab$", "\\1", artificial.bounds.filenames );
    
    the.artificial.bounds <- lapply( names( artificial.bounds.filenames ), function ( .artificial.bounds.name ) {
        .tbl <- read.table( artificial.bounds.filenames[[ .artificial.bounds.name ]], header = TRUE, sep = "\t" );
     # Special: for v3, only use caprisa seqs (not rv217, for now).
        if( the.region == "v3" ) {
            .tbl <-
                .tbl[ grep( "^100\\d\\d\\d", rownames( .tbl ) ), , drop = FALSE ];
        } else if( the.region == "rv217_v3" ) {
            .tbl <-
                .tbl[ grep( "^100\\d\\d\\d", rownames( .tbl ), invert = TRUE ), , drop = FALSE ];
        }
        
        return( .tbl );
    } );
    names( the.artificial.bounds ) <- names( artificial.bounds.filenames );

    ## Special: exclude the deterministic bounds, for now.
    the.artificial.bounds <-
        the.artificial.bounds[ grep( "deterministic", names( the.artificial.bounds ), invert = TRUE ) ];
    
    return( the.artificial.bounds );
} # get.artificial.bounds (..)


### R code from vignette source '/Users/Paul/src/from-git/hiv-founder-id/DSStarPhyTest.Rnw'
### Encoding: ASCII

###################################################
### code chunk number 1: DSStarPhyTest.Rnw:18-30
###################################################
## R packages needed
#library( "xtable" )
#library( "coin" )
#library( "ggplot2" )
#library( "binom" )

library(tools)
# For running this in the same style as PFitter.R
args = commandArgs(trailingOnly=TRUE)

# Setup for prettier Sweave output.
old.continue.option <- options( continue = " " )


###################################################
### code chunk number 2: DSStarPhyTest.Rnw:78-81
###################################################
    # MAGIC N: DS.NDRAWS
    DS.NDRAWS <- 1000;



###################################################
### code chunk number 3: DSStarPhyTest.Rnw:84-405
###################################################
calculateSumOfDistancesAccountingForSequenceMultiplicity <- function ( intersequence.dlist ) {
    .mult.1 <- suppressWarnings( as.numeric( gsub( "^.*_(\\d+)$", "\\1", intersequence.dlist[ , 1 ] ) ) );
    .mult.1[ is.na( .mult.1 ) | ( .mult.1 == 0 ) ] <- 1; # if it is missing the _nnn, just call it 1.
    .mult.2 <- suppressWarnings( as.numeric( gsub( "^.*_(\\d+)$", "\\1", intersequence.dlist[ , 2 ] ) ) );

    .mult.2[ is.na( .mult.2 ) | ( .mult.2 == 0 ) ] <- 1; # if it is missing the _nnn, just call it 1.
    return( list( sum = sum( intersequence.dlist[ , 3 ] * .mult.1 * .mult.2, na.rm = T ), count = sum( .mult.1 * .mult.2, na.rm = T ) ) ); 
} # calculateSumOfDistancesAccountingForSequenceMultiplicity (..)

replicateDistancesForSequenceMultiplicity <- function ( any.dlist, include.intersequence.distances.of.duplicates = FALSE, missing.seqnames = NULL ) {
    # Note that this isn't the entire thing because we also have to add the "0" distance distances among the duplicated seqs; see below.
    #print( missing.seqnames );
    if( is.null( any.dlist ) || is.null( dim( any.dlist ) ) || ( nrow( any.dlist ) == 0 ) ) {
        maybe.longer.dlist.list <- NULL;
    } else {
      maybe.longer.dlist.list <-
      lapply( 1:nrow( any.dlist ), function( .row.i ) {
        .row <- any.dlist[ .row.i, ];
        .mult.1 <- NA;
        if( length( grep( "^.*_(\\d+)$", .row[ 1 ], perl = TRUE ) ) > 0 ) {
            .mult.1 <- as.numeric( gsub( "^.*_(\\d+)$", "\\1", .row[ 1 ], perl = TRUE ) );
        }
        if( is.na( .mult.1 ) ) {
            .mult.1 <- 1; # if it is missing the _nnn, just call it 1.
        }
        .mult.2 <- NA;
        if( length( grep( "^.*_(\\d+)$", .row[ 2 ], perl = TRUE ) ) > 0 ) {
            .mult.2 <- as.numeric( gsub( "^.*_(\\d+)$", "\\1", .row[ 2 ], perl = TRUE ) );
        }
        if( is.na( .mult.2 ) ) {
            .mult.2 <- 1; # if it is missing the _nnn, just call it 1.
        }
        if( ( .mult.1 * .mult.2 ) > 1 ) {
          .rv <- matrix( "", ncol = 3, nrow = ( .mult.1 * .mult.2 ) );
          if( .mult.1 > 1 ) {
              .rv[ , 1 ] <- paste( .row[[ 1 ]], "_replicate", 1:.mult.1, "_1", sep = "" );
          } else {
              # If they end in _###, replace that number with 1.
              .rv[ , 1 ] <- gsub( "(^.+)_(\\d+)$", "\\1_1", .row[[ 1 ]] );
          }
          if( .mult.2 > 1 ) {
              .rv[ , 2 ] <- paste( .rv[ , 2 ], "_replicate", 1:.mult.2, "_1", sep = "" );
          } else {
              # If they end in _###, replace that number with 1.
              .rv[ , 2 ] <- gsub( "(^.+)_(\\d+)$", "\\1_1", .row[[ 2 ]] );
          }
          .rv[ , 2 ] <- paste( .row[[ 1 ]], "_replicate", 1:.mult.2, "_1", sep = "" )
          .rv[ , 3 ] <- .row[[ 3 ]];
          ## TODO: REMOVE
          # print( c( .mult.1, .mult.2, nrow( .rv ) ) );
          return( .rv );
        } else {
          .rv <- matrix( .row, nrow = 1 );
          mode( .rv ) <- "character";
          return( .rv );
        }
      } );
    }
    if( !is.null( maybe.longer.dlist.list ) && ( length( maybe.longer.dlist.list ) > 0 ) ) {
      maybe.longer.dlist <- do.call( rbind, maybe.longer.dlist.list );
      colnames( maybe.longer.dlist ) <- colnames( any.dlist );
    } else {
      maybe.longer.dlist <- NULL;
    }
    if( include.intersequence.distances.of.duplicates ) {
      duplicated.entries.dlist.list <-
        lapply( unique( c( any.dlist[ , 1 ], any.dlist[ , 2 ], missing.seqnames ) ), function( .seq.name ) {
            .mult <- NA;
            if( length( grep( "^.*_(\\d+)$", .seq.name, perl = TRUE ) ) > 0 ) {
                .mult <- as.numeric( gsub( "^.*_(\\d+)$", "\\1", .seq.name, perl = TRUE ) );
            }
            if( is.na( .mult ) ) {
                .mult <- 1; # if it is missing the _nnn, just call it 1.
            }
            if( .mult == 1 ) {
                return();
            }
            .rv <- matrix( "", ncol = 3, nrow = choose( .mult, 2 ) );
            .rv[ , 1 ] <- .seq.name;
            .rv[ , 2 ] <- .seq.name;
            .rv[ , 3 ] <- 0;
            return( .rv );
      } );
      if( !all( unlist( lapply( duplicated.entries.dlist.list, is.null ) ) ) ) {
        duplicated.entries.dlist <- do.call( rbind, duplicated.entries.dlist.list );
        if( is.null( dim( duplicated.entries.dlist ) ) ) {
            duplicated.entries.dlist <- matrix( duplicated.entries.dlist, nrow = 1, ncol = 3 );
            mode( duplicated.entries.dlist ) <- "character";
        }
        colnames( duplicated.entries.dlist ) <- colnames( any.dlist );
        if( is.null( maybe.longer.dlist ) ) {
          return( duplicated.entries.dlist );
        } else {
          # Put zero entries first to maintian sort order jic.
          return( rbind( duplicated.entries.dlist, maybe.longer.dlist ) );
        }
      } else {
        return( maybe.longer.dlist );
      }
    } else {
      return( maybe.longer.dlist );
    } # If include.intersequence.distances.of.duplicates ..
} # replicateDistancesForSequenceMultiplicity (..)

makeDList <- function ( consensus.distances, intersequence.distances ) {
    dlist <- matrix( NA, nrow = length( consensus.distances ) + length( intersequence.distances ), ncol = 3 );
    colnames( dlist ) <- c( "from.seq", "to.seq", "dist" );
    dlist[ 1:length( consensus.distances ), 1 ] <- "Consensus";
    dlist[ 1:length( consensus.distances ), 2 ] <- names( consensus.distances );
    dlist[ 1:length( consensus.distances ), 3 ] <- consensus.distances;
    fromnames <- gsub( "^(.+) to .+$", "\\1", names( intersequence.distances ) );
    tonames <- gsub( "^.+ to (.+)$", "\\1", names( intersequence.distances ) );
    dlist[ length( consensus.distances ) + 1:length( intersequence.distances ), 1 ] <- fromnames;
    dlist[ length( consensus.distances ) + 1:length( intersequence.distances ), 2 ] <- tonames;
    dlist[ length( consensus.distances ) + 1:length( intersequence.distances ), 3 ] <- intersequence.distances;
    return( as.data.frame( dlist ) );
} # makeDList (..)

# Calculate and return the approximate one-sided (upper) p-value using the given null statistics (usually computed through permutation).
calculateUpperSidedPValue <- function ( the.stat, null.stats, add.one.for.observed.test.statistic = TRUE, use.ci = TRUE ) {
    if( is.na( the.stat ) ) {
        return( NA );
    } else {
        .nulls.more.extreme.than.the.stat <-
            sum( null.stats >= the.stat, na.rm=T )
            + as.numeric( add.one.for.observed.test.statistic );
        .total.nulls <-
            sum( !is.na( null.stats ) ) +
                as.numeric( add.one.for.observed.test.statistic );
        if( use.ci ) {
            if( !require( "binom" ) ) {
                install.packages( "binom", dependencies = TRUE );
                if( !require( "binom" ) ) {
                    stop( "ERROR LOADING \"binom\" package" );
                }
            }
            return( binom.confint( .nulls.more.extreme.than.the.stat, .total.nulls, conf.level = 0.95, methods = "wilson" )$upper );
        } else {
            return( .nulls.more.extreme.than.the.stat / .total.nulls );
        }
    }
} # calculateUpperSidedPValue (..)

DSStarPhyTest <- function (
  infile = args[1],
  dlist = read.table( file=infile, sep="\t", stringsAsFactors=F ),
  epsilon = c(as.numeric(args[2])), # 2.16e-05
  nbases = c(as.numeric(args[3])),
  be.verbose = TRUE,
  MAXIMUM.DS.SAMPLE.SIZE = 2500
) {
    ### Consensus distances.
    ## Sort first for efficiency.  Then expand.
    .con.mat <- dlist[ which(dlist[,1]==dlist[1,1]), , drop = FALSE ];
    .con.mat <- .con.mat[ order( .con.mat[ , 3 ] ), , drop = FALSE ]; # Sort it by column 3, distance.
    rownames( .con.mat ) <- .con.mat[ , 2 ];
    sorted.consensus.distances <-
        as.numeric( replicateDistancesForSequenceMultiplicity( .con.mat, include.intersequence.distances.of.duplicates = FALSE )[ , 3 ] );
    consensus.distances <-
        sorted.consensus.distances;
    
    ### Intersequence distances.
    ## Sort first for efficiency.  Then expand.
    .mat <- dlist[ -which(dlist[,1]==dlist[1,1]), , drop = FALSE ];
    .mat <- .mat[ order( .mat[ , 3 ] ), , drop = FALSE ]; # Sort it by column 3, distance.
    rownames( .mat ) <- apply( .mat[ , 1:2 ], 1, paste, collapse = " to " );
    sorted.intersequence.distances <-
        as.numeric( replicateDistancesForSequenceMultiplicity( .mat, include.intersequence.distances.of.duplicates = TRUE, missing.seqnames = .con.mat[ , 2 ] )[ , 3 ] );
    intersequence.distances <-
        sorted.intersequence.distances;

    # Do this before reducing the sample size..
    DS.lambda <- mean( sorted.intersequence.distances, na.rm = TRUE );

    ## Reduce the sample sizes to make it computable.
    .num.consensus.distances <- length( sorted.consensus.distances );
    if( .num.consensus.distances > MAXIMUM.DS.SAMPLE.SIZE ) {
        new.sorted.consensus.distances.table <- round( MAXIMUM.DS.SAMPLE.SIZE * table( sorted.consensus.distances ) / .num.consensus.distances )
        sorted.consensus.distances <-
            unlist( sapply( 1:length( new.sorted.consensus.distances.table ), function( .i ) { rep( .i - 1, new.sorted.consensus.distances.table[ .i ]  ) } ) );
    }
    .num.intersequence.distances <- length( sorted.intersequence.distances );
    if( .num.intersequence.distances > MAXIMUM.DS.SAMPLE.SIZE ) {
        new.sorted.intersequence.distances.table <- round( MAXIMUM.DS.SAMPLE.SIZE * table( sorted.intersequence.distances ) / .num.intersequence.distances )
        sorted.intersequence.distances <-
            unlist( sapply( 1:length( new.sorted.intersequence.distances.table ), function( .i ) { rep( .i - 1, new.sorted.intersequence.distances.table[ .i ]  ) } ) );
    }
    
    #print( DS.lambda );
    
    draw.one.nonconflicted.pepr.sample <- function ( sorted.observed.data, conflict.max = 5000, method = "fast" ) {
        ## The naive way, just keep drawing until one works.
          found.it <- FALSE;
          conflict <- 0;
          while( !found.it && conflict < conflict.max ) {
            # Directly draw from the Poisson DSM
            lambda.low <-
                rgamma( n = length( sorted.observed.data ), shape = sorted.observed.data );
              lambda.width <-
                  rgamma( n = length( sorted.observed.data ), shape = 1 );
              if( any( max( lambda.low ) > ( lambda.low + lambda.width ) ) ) {
                  if( method == "naive" ) {
                    ## CONFLICT.
                    conflict <- conflict + 1;
                    cat( "." );
                    next;
                  } else {
                    # All of them need to overlap, so every value of lambda.width must be at least as long as max( lambda.low ) - lambda.low.  They are iid, so the condition that any one of them is at least that value has no effect on the rest of them.
                    # But also, the distribution of an exponential is memoryless, so we can just shift everything. Note we do need to draw them again so they are not truncated (because we've conditioned on them being a minimum, but they're memoryless..)
                      ## ERE I AM.  This should just be done outright; the point is right, easily explained using memorylessness.
                      lambda.width <-
                          ( max( lambda.low ) - lambda.low ) +
                          rgamma( n = length( sorted.observed.data ), shape = 1 );
                  } # End if method == "naive" .. else ..
              } # End if the draws aren't long enough.
              stopifnot( !any( max( lambda.low ) > ( lambda.low + lambda.width ) ) );
            found.it <- 1;
          } # end while( !found.it )
          if( !found.it ) { stop( "TOO MUCH CONFLICT!" ) }
          
        return( list( lambda.low = lambda.low, lambda.width = lambda.width ) );
    } # draw.one.nonconflicted.pepr.sample 
    
    draw.once.from.each.and.evaluate.assertion <- function( assertion.is.that.1.is.x.times.2.low = 1.5, assertion.is.that.1.is.x.times.2.high = 2.5, sorted.observed.data.1 = sorted.consensus.distances, sorted.observed.data.2 = sorted.intersequence.distances, be.verbose = FALSE, weakening.type = NA ) {
        stopifnot( length( assertion.is.that.1.is.x.times.2.low ) == length( assertion.is.that.1.is.x.times.2.high ) );
        stopifnot( all( assertion.is.that.1.is.x.times.2.low <= assertion.is.that.1.is.x.times.2.high ) );
        .sample1 <- draw.one.nonconflicted.pepr.sample( sorted.observed.data.1 );
        if( is.na( weakening.type ) ) {
            .weakening.dependent.part.high <- min( .sample1$lambda.low + .sample1$lambda.width );
        } else if( weakening.type == "rotate" ) {
            .all.width.rotations <- sapply( 1:length( .sample1$lambda.width ), function( .start.i ) { c( .sample1$lambda.width, .sample1$lambda.width )[ .start.i + 0:( length( .sample1$lambda.width ) - 1 ) ] } );
            .weakening.dependent.part.high <- max( apply( .all.width.rotations, 1, function( .widths ) { min( .sample1$lambda.low + .widths ) } ) );
        } else if( weakening.type == "complete" ) {
            .weakening.dependent.part.high <- min( sort( .sample1$lambda.low ) + sort( .sample1$lambda.width, decreasing = TRUE ) );
        } else {
            stop( paste( "unknown weakening type:", weakening.type ) );
        }
        .sample1.max = lapply( 1:length( assertion.is.that.1.is.x.times.2.low ), function( .i ) {
            assertion.is.that.1.is.x.times.2.high[ .i ] * .weakening.dependent.part.high;
        } );
        .weakening.dependent.part.low <- max( .sample1$lambda.low ); # Not actually dependent on how we weaken.
        .sample1.min = lapply( 1:length( assertion.is.that.1.is.x.times.2.low ), function( .i ) {
            assertion.is.that.1.is.x.times.2.low[ .i ] * .weakening.dependent.part.low;
        } );
        stopifnot( all( unlist( lapply( 1:length( assertion.is.that.1.is.x.times.2.low ), function( .i ) { .sample1.min[[ .i ]] <= .sample1.max[[ .i ]] } ) ) ) );
        if( FALSE && be.verbose ) {
            cat( paste( sapply( 1:length( .sample1.min ), function( .i ) { paste( c( ".sample1.min", .sample1.min[[ .i ]],  ".sample1.max", .sample1.max[[ .i ]] ), collapse = " " ) } ), collapse = "\n" ), fill = T );
        }
        .sample2 <- draw.one.nonconflicted.pepr.sample( sorted.observed.data.2 );
        .sample2.max = min( .sample2$lambda.low + .sample2$lambda.width );
        .sample2.min = max( .sample2$lambda.low );
        stopifnot( .sample2.min <= .sample2.max );
        if( FALSE && be.verbose ) {
            cat( paste( c( ".sample2.min", .sample2.min,  ".sample2.max", .sample2.max ), collapse = " " ), fill = T );
        }
        .rv <-
            sapply( 1:length( .sample1.min ), function( .i ) {
        if( ( .sample1.max[[ .i ]] < .sample2.min ) ||
                ( .sample1.min[[ .i ]] > .sample2.max ) ) {
            # Evidence against.
            if( FALSE && be.verbose ) {
                cat( "Q", fill = T );
            }
            return( list( P = 0, Q = 1, R = 0 ) );
        }
        if( ( .sample1.min[[ .i ]] >= .sample2.min ) && ( .sample1.max[[ .i ]] <= .sample2.max ) ) {
            # Evidence in support.
            if( FALSE && be.verbose ) {
                cat( "P", fill = T );
            }
            return( list( P = 1, Q = 0, R = 0 ) );
        }
        # If it's not evidence for or against, then we are in a "don't know" situation.
        if( FALSE && be.verbose ) {
            cat( "R", fill = T );
        }
        return( list( P = 0, Q = 0, R = 1 ) );
            } );
        if( is.null( dim( .rv ) ) ) {
            .rv <- matrix( .rv, nrow = 3 );
        }
        rownames( .rv ) <- c( "P", "Q", "R" );
        return( .rv );
    } # draw.once.from.each.and.evaluate.assertion (..)

    # GET A BUNCH OF DRAWS
    .result <- sapply( 1:DS.NDRAWS, function( .i ) { draw.once.from.each.and.evaluate.assertion( be.verbose = be.verbose ) } );
    rownames( .result ) <- c( "P", "Q", "R" );
    DS.P <- mean( unlist( .result[ "P", ] ) );
    DS.Q <- mean( unlist( .result[ "Q", ] ) );
    DS.R <- mean( unlist( .result[ "R", ] ) );
    if( DS.P == 0 ) {
        DS.Ptext <- paste( "<=", ( 1 / DS.NDRAWS ) );
    } else {
        DS.Ptext <- sprintf( paste( "= %0.", ceiling( log10( DS.NDRAWS ) ), "f", sep = "" ), DS.P )
    }
    if( DS.Q == 0 ) {
        DS.Qtext <- paste( "<=", ( 1 / DS.NDRAWS ) );
    } else {
        DS.Qtext <- sprintf( paste( "= %0.", ceiling( log10( DS.NDRAWS ) ), "f", sep = "" ), DS.Q )
    }
    if( DS.R == 0 ) {
        DS.Rtext <- paste( "<=", ( 1 / DS.NDRAWS ) );
    } else {
        DS.Rtext <- sprintf( paste( "= %0.", ceiling( log10( DS.NDRAWS ) ), "f", sep = "" ), DS.R )
    }
    #print( paste( "The DS evidence against the assertion that the Poisson rate between sequences is twice the rate of sequences to the consensus is Q ", DS.Qtext, ". The remaining evidence (", sprintf( paste( "%0.", ceiling( log10( DS.NDRAWS ) ), "f", sep = "" ), DS.R ), ") neither supports nor contradicts the assertion.", sep = "" ) );
    if( be.verbose ) {
      if( DS.Q >= 0.95 ) {
          cat( "DSStarPhyTest that intersequence rate = 2 x seq-consensus rate: BAD", fill = TRUE );
          cat( paste( "\tThere is evidence against the assertion that the Poisson rate between sequences is between 1.5 and 2.5 times the rate of sequences to the consensus (P ", DS.Ptext, ", Q ", DS.Qtext, ", R ", DS.Rtext, ").", sep = "" ), fill = TRUE );
      } else {
          cat( "DSStarPhyTest that intersequence rate = 2 x seq-consensus rate: OK", fill = TRUE );
          cat( paste( "\tThere is not sufficient evidence against the assertion that the Poisson rate between sequences is between 1.5 and 2.5 times the rate of sequences to the consensus (P ", DS.Ptext, ", Q ", DS.Qtext, ", R ", DS.Rtext, ").", sep = "" ), fill = TRUE );
      }
    }
    dspfitter.results <- list( twox = list( "P" = DS.P, "Q" = DS.Q, "R" = DS.R ) );

    return( dspfitter.results );
} # DSStarPhyTest (..)



###################################################
### code chunk number 4: DSStarPhyTest.Rnw:417-418
###################################################
.result.ignored <- DSStarPhyTest( be.verbose = TRUE );


###################################################
### code chunk number 5: DSStarPhyTest.Rnw:425-427
###################################################
# (un)Setup for prettier Sweave output.
options( continue = old.continue.option$continue )



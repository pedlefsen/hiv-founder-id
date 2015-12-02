library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment", "consensus"

# for maskSynonymousCodonsInAlignedFasta(..)
source( "maskSynonymousCodonsInAlignedFasta_safetosource.R" );

# for removeDuplicateSequencesFromAlignedFasta(..)
source( "removeDuplicateSequencesFromAlignedFasta_safetosource.R" )

## Similar to runPoissonFitter, except that the distances are computed within each cluster, then put together.
## Input files are taken to be all files matching the given pattern (default: ${fasta.file.prefix}.cluster\d+.fasta) -- you can change the suffix but the prefix must be fasta.file.prefix.
# Compute Hamming distances, prepare inputs to PFitter.R, call PFitter.R.
runMultiFounderPoissonFitter <- function ( fasta.file.prefix, output.dir = NULL, include.gaps.in.Hamming = FALSE, fasta.file.suffix.pattern = "\\_cluster\\d+\\.fasta", output.dir.suffix = "_MultiFounderPoissonFitterDir", pairwise.hamming.distances.file.suffix = "_multiFounderPairwiseHammingDistances.txt", run.DSPFitter = FALSE, maskOutNonsynonymousCodons = FALSE ) {

    if( length( grep( "^(.*?)\\/[^\\/]+$", fasta.file.prefix ) ) == 0 ) {
        fasta.file.prefix.path <- ".";
    } else {
        fasta.file.prefix.path <-
            gsub( "^(.*?)\\/[^\\/]+$", "\\1", fasta.file.prefix );
    }
    fasta.file.prefix.short <-
        gsub( "^.*?\\/?([^\\/]+?)$", "\\1", fasta.file.prefix, perl = TRUE );
    fasta.file.prefix.short.nosuffix <-
        gsub( "^(.*?)\\.[^\\.]+$", "\\1", fasta.file.prefix.short, perl = TRUE );

  if( is.null( output.dir ) ) {
      output.dir <- fasta.file.prefix.path;
  } else {
      output.dir <- output.dir;
  }
  ## Remove "/" from end of output.dir
  output.dir <-
      gsub( "^(.*?)\\/+$", "\\1", output.dir );

  if( maskOutNonsynonymousCodons ) {
    fasta.file.prefix.short.nosuffix <-
        paste( fasta.file.prefix.short.nosuffix, "_maskNonsynonymousCodons", sep = "" );
    ## TODO: REMOVE
    #cat( "maskOutNonsynonymousCodons. fasta.file.prefix.short.nosuffix is ", fasta.file.prefix.short.nosuffix, fill = T );
  }
    
    output.dir <- paste( output.dir, "/", fasta.file.prefix.short.nosuffix, output.dir.suffix, sep = "" );
    
    unlink( output.dir, recursive = TRUE );
    dir.create( output.dir, showWarnings = TRUE );

    # cat( "fasta.file.prefix.path: ", fasta.file.prefix.path, sep = "", fill = TRUE );
    .dir <- dir( fasta.file.prefix.path, full.names = FALSE );
    #cat( "dir(fasta.file.prefix.path): ", paste( .dir, collapse = ", " ), sep = "", fill = TRUE );
    fasta.files <- grep( paste( "^", fasta.file.prefix.short, "(", fasta.file.suffix.pattern, ")$", sep = "" ), .dir, value = TRUE, perl = TRUE );
    #cat( "grep( \"", paste( "(", fasta.file.suffix.pattern, ")$", sep = "" ), "\", dir(fasta.file.prefix.path) ): ", paste( .dir.withsuffix, collapse = ", " ), sep = "", fill = TRUE );
    #cat( "grep( \"", paste( "^", fasta.file.prefix.short, sep = "" ), "\", grep( \"", paste( "(", fasta.file.suffix.pattern, ")$", sep = "" ), "\", dir(fasta.file.prefix.path) ) ): ", paste( fasta.files, collapse = ", " ), sep = "", fill = TRUE );
    ## TODO: REMOVE
    # warning( paste( "MATCHING suffix:", fasta.file.suffix.pattern, sep = "" ) );
    # warning( paste( "and MATCHING prefix:", fasta.file.prefix.short, sep = "" ) );
    #warning( fasta.files );
    stopifnot( length( fasta.files ) > 1 );

    seq.lengths <- c();
    seq.counts <- c();
    pairwise.distances.as.matrix.flat.by.fasta.file <- lapply( fasta.files, function ( fasta.file ) {
      fasta.file <- paste( fasta.file.prefix.path, fasta.file, sep = "/" );
      if( maskOutNonsynonymousCodons ) {
          fasta.file <-
              maskSynonymousCodonsInAlignedFasta( fasta.file, mask.nonsynonymous = TRUE );
      }
      ## If there are any duplicate sequences, remove them and
      ## incorporate the number of identical sequences into their names (_nnn suffix).
      ## NOTE that the value of include.gaps.in.Hamming is passed to removeDuplicateSequencesFromAlignedFasta, which means that by when this is FALSE (as it is by default), sequences will be considered duplicates even if there are gaps between them.
      # The output has the file name of the consensus file.
      fasta.file.no.duplicates <-
          removeDuplicateSequencesFromAlignedFasta( fasta.file, output.dir, add.copy.number.to.sequence.names = TRUE, include.gaps.in.Hamming = include.gaps.in.Hamming );
      fasta.file.no.duplicates.short <-
          gsub( "^.*?\\/?([^\\/]+?)$", "\\1", fasta.file.no.duplicates, perl = TRUE );
      fasta.file.no.duplicates.short.nosuffix <-
          gsub( "^(.*?)\\.[^\\.]+$", "\\1", fasta.file.no.duplicates.short, perl = TRUE );
      in.fasta.no.duplicates <- read.dna( fasta.file.no.duplicates, format = "fasta" );
  
      ## This is the one with duplicates intact; lazy way to get consensus.
      in.fasta <- read.dna( fasta.file, format = "fasta" );
  
      # Add the consensus to the one with no duplicates.
      .consensus.mat <- matrix( seqinr::consensus( as.character( in.fasta ) ), nrow = 1 );
      consensus <- as.DNAbin( .consensus.mat );
      rownames( consensus ) <- paste( fasta.file.prefix.short.nosuffix, "Consensus", sep = "." );
      #seq.length <- ncol( consensus );
      # Counting everything for which the consensus is not a gap:
      seq.length <- sum( .consensus.mat[ 1, ] != '-' );
      # Counting everything that is not _entirely_ a gap:
      #.consensus.mat.profile <- seqinr::consensus( as.character( in.fasta ), method = "profile" );
      #seq.length <- sum( apply( .consensus.mat.profile, 2, function( .counts ) { !all( .counts[ -1 ] == 0 ) } ) );
      seq.lengths <<- c( seq.lengths, seq.length );
      seq.counts <<- c( seq.counts, nrow( in.fasta ) );
  
      ### NEW, MUCH FASTER:
      fasta.no.duplicates.with.consensus <-
          rbind( consensus, in.fasta.no.duplicates );
      ### OLD, SLOW (not excluding duplicates ):
#       fasta.no.duplicates.with.consensus <-
#           rbind( consensus, in.fasta );
  
      # Remove any columns with a consensus that is a gap, which means
      # that over half of seqs have gaps.  This needs to be removed
      # because it becomes not sensible to consider poisson rates of
      # insertions.  We do however consider rates of deletions, treating
      # them as point mutations (by including the indel counts in the
      # Hamming distance calculation).
      fasta.no.duplicates.with.consensus <- fasta.no.duplicates.with.consensus[ , .consensus.mat[ 1, ] != "-" ];
      
      # The pairwise.deletion = TRUE argument is necessary so that columns with any gaps are not removed.
      # The optional second call adds in the count of sites with gap differences
      # (gap in one sequence but not the other), which are not included
      # in the first call's "raw" count. Adding them back in, we are
      # effectively treating those as differences.
      fasta.no.duplicates.with.consensus.dist <- dist.dna( fasta.no.duplicates.with.consensus, model = "N", pairwise.deletion = TRUE );
      fasta.no.duplicates.with.consensus.dist[ is.nan( fasta.no.duplicates.with.consensus.dist ) ] <- 0;
      if( include.gaps.in.Hamming ) {
          fasta.no.duplicates.with.consensus.dist <- fasta.no.duplicates.with.consensus.dist + dist.dna( fasta.no.duplicates.with.consensus, model = "indel", pairwise.deletion = TRUE );
      }
      
      if( any( is.null( fasta.no.duplicates.with.consensus.dist ) ) || any( is.na( fasta.no.duplicates.with.consensus.dist ) ) || any( !is.finite( fasta.no.duplicates.with.consensus.dist ) ) ) {
        ## TODO: REMOVE
        warning( "UH OH got illegal distance value" );
        print( "UH OH got illegal distance value" );
        print( fasta.no.duplicates.with.consensus.dist );
      }
  
      pairwise.distances.as.matrix <-
          as.matrix( fasta.no.duplicates.with.consensus.dist );
  
      pairwise.distances.as.matrix.flat <-
          matrix( "", nrow = ( ( nrow( pairwise.distances.as.matrix ) * ( ncol( pairwise.distances.as.matrix ) - 1 ) ) / 2 ), ncol = 3 );
      line.i <- 1;
      for( row.i in 1:( nrow( pairwise.distances.as.matrix ) - 1 ) ) {
          for( col.i in ( row.i + 1 ):ncol( pairwise.distances.as.matrix ) ) {
              if( row.i == 1 ) { # consensus, no _1 (multiplicity of the observed sequence)
                  pairwise.distances.as.matrix.flat[ line.i, ] <-
                      c( rownames( pairwise.distances.as.matrix )[ row.i ], colnames( pairwise.distances.as.matrix )[ col.i ], pairwise.distances.as.matrix[ row.i, col.i ] );
              } else {
                  pairwise.distances.as.matrix.flat[ line.i, ] <-
                      c( rownames( pairwise.distances.as.matrix )[ row.i ], colnames( pairwise.distances.as.matrix )[ col.i ], pairwise.distances.as.matrix[ row.i, col.i ] );
              }
              line.i <- line.i + 1;
          }
      }
      return( pairwise.distances.as.matrix.flat );
    } );

    # Ok now merge them.  It suffices to just rbind them.
    pairwise.distances.as.matrix.flat <- do.call( rbind, pairwise.distances.as.matrix.flat.by.fasta.file );
    
    pairwise.distances.as.matrix.file <- paste( output.dir, "/", fasta.file.prefix.short.nosuffix, pairwise.hamming.distances.file.suffix, sep = "" );
    write.table( pairwise.distances.as.matrix.flat, file = pairwise.distances.as.matrix.file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE );

    # Seq lengths may vary.  To make it work with PFitter, which needs just one vlaue here, use the weighted average, where the weight is in terms of evidence for days (which comes only from the intersequence distances).
    .weight.num <- choose( seq.counts, 2 );
    .weight.denom <- sum( .weight.num );
    seq.length <- sum( seq.lengths * ( .weight.num / .weight.denom ) );
    ## TODO: REMOVE
    print( paste( "seq lengths:", paste( seq.lengths, collapse = ", " ), "\\n" ) );
    print( paste( "seq weights:", paste( ( .weight.num / .weight.denom ), collapse = ", " ), "\\n" ) );
    print( paste( "seq.length:", seq.length, "\\n" ) );

    if( all( pairwise.distances.as.matrix.flat[ , 3 ] == 0, ra.rm = TRUE ) ) {
        # Uh oh.  PFitter doesn't handle this case very well.  Instead of running PFitter, write out files communicating that this is a degenerate situation.
        outfile <- paste( output.dir, "/", "LOG_LIKELIHOOD.results.txt", sep="" );
        write( paste("Sample", "Lambda", "St.Dev", "NSeq", "NBases", "MeanHD", "MaxHD","Days(CI)", "Chi2","DF","Goodness_of_pval", sep="\t"), file=outfile, append=FALSE );
        write( paste( fasta.file.prefix.short.nosuffix, format( 0, digits=4 ), format( 0, digits=4 ), sum( seq.counts[ .weight.num > 0 ] ), seq.length, format(0, digits=2), 0, "0 (0, 1)", format(0, digits=4), 0, format(1, digits=4), sep="\t"), file=outfile, append=TRUE );
        outfile2 <- paste( output.dir, "/", "CONVOLUTION.results.txt", sep="" );
        write( paste( fasta.file.prefix.short.nosuffix, "FOLLOWS A STAR-PHYLOGENY", sep = " " ), file=outfile2, append=FALSE );
        .rv <- "0";
    } else {
      R.cmd <- paste( "R CMD BATCH '--vanilla --args", pairwise.distances.as.matrix.file, "2.16e-05", seq.length, "' PFitter.R" );
      .rv <- system( R.cmd );
    }

    if( run.DSPFitter ) {
        DSPFitter.outfile <- paste( output.dir, "/", fasta.file.prefix.short.nosuffix, "_DSPFitter.out", sep = "" );
        #warning( paste( "DSPFitter.outfile:", DSPFitter.outfile ) );
        R.cmd <- paste( "R CMD BATCH '--vanilla --args", pairwise.distances.as.matrix.file, "2.16e-05", seq.length, "' DSPFitter.R", DSPFitter.outfile );
        .rv <- system( R.cmd );
    }
    return( .rv );
} # runMultiFounderPoissonFitter ( fasta.file.prefix, output.dir )

## Here is where the action is.
fasta.file.prefix <- Sys.getenv( "runMultiFounderPoissonFitter_inputFilenamePrefix" ); # alignment
output.dir <- Sys.getenv( "runMultiFounderPoissonFitter_outputDir" ); # NOTE: will create a subdir for the output, named after the fasta file, with "MultiFounderPoissonFitterDir".
suffix.pattern <- Sys.getenv( "runMultiFounderPoissonFitter_suffixPattern" );
if( suffix.pattern == "" ) {
    suffix.pattern <- NULL;
}
run.DSPFitter <- Sys.getenv( "runMultiFounderPoissonFitter_runDSPFitter" ); # NOTE: will create a subdir for the output, named after the fasta file, with "PoissonFitterDir".
if( ( run.DSPFitter == "" ) || ( toupper( run.DSPFitter ) == "F" ) || ( toupper( run.DSPFitter ) == "FALSE" ) || ( run.DSPFitter == "0" ) ) {
    run.DSPFitter = FALSE;
} else {
    run.DSPFitter = TRUE;
}
maskOutNonsynonymousCodons <- Sys.getenv( "runMultiFounderPoissonFitter_maskOutNonsynonymousCodons" );
if( ( maskOutNonsynonymousCodons == "" ) || ( toupper( maskOutNonsynonymousCodons ) == "F" ) || ( toupper( maskOutNonsynonymousCodons ) == "FALSE" ) || ( maskOutNonsynonymousCodons == "0" ) ) {
    maskOutNonsynonymousCodons = FALSE;
} else {
    maskOutNonsynonymousCodons = TRUE;
}

## TODO: REMOVE
#  warning( paste( "alignment input file prefix:", fasta.file.prefix ) );
#  warning( paste( "output dir:", output.dir ) );
#  warning( paste( "suffix.pattern", suffix.pattern ) );
# warning( paste( "run DSPFitter:", run.DSPFitter ) );
# warning( paste( "maskOutNonsynonymousCodons:", maskOutNonsynonymousCodons ) );

if( !is.null( suffix.pattern ) ) {
    print( runMultiFounderPoissonFitter( fasta.file.prefix, output.dir, fasta.file.suffix.pattern = suffix.pattern, output.dir.suffix = "_MultiRegionPoissonFitterDir", pairwise.hamming.distances.file.suffix = "_multiRegionPairwiseHammingDistances.txt", run.DSPFitter = run.DSPFitter, maskOutNonsynonymousCodons = maskOutNonsynonymousCodons ) );
} else {
    print( runMultiFounderPoissonFitter( fasta.file.prefix, output.dir, run.DSPFitter = run.DSPFitter, maskOutNonsynonymousCodons = maskOutNonsynonymousCodons ) );
}

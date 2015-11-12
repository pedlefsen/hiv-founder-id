library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment", "consensus"

## Similar to runPoissonFitter, except that the distances are computed within each cluster, then put together.
## Input files are taken to be all files matching the given pattern (default: ${fasta.file.prefix}.cluster\d+.fasta) -- you can change the suffix but the prefix must be fasta.file.prefix.
# Compute Hamming distances, prepare inputs to PFitter.R, call PFitter.R.
runMultiFounderPoissonFitter <- function ( fasta.file.prefix, output.dir = NULL, include.gaps.in.Hamming = FALSE, fasta.file.suffix.pattern = "\\.cluster\\d+\\.fasta", output.dir.suffix = "_MultiFounderPoissonFitterDir", pairwise.hamming.distances.file.suffix = "_multiFounderPairwiseHammingDistances.txt" ) {

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

    output.dir <- paste( output.dir, "/", fasta.file.prefix.short.nosuffix, output.dir.suffix, sep = "" );
    
    unlink( output.dir, recursive = TRUE );
    dir.create( output.dir, showWarnings = TRUE );

    fasta.files <- grep( paste( "^", fasta.file.prefix.short, sep = "" ), grep( paste( "(", fasta.file.suffix.pattern, ")$", sep = "" ), dir( fasta.file.prefix.path, full.names = FALSE ), value = T, perl = TRUE ), value = T );
    ## TODO: REMOVE
    #print( paste( "MATCHING: ", fasta.file.prefix.short, fasta.file.suffix.pattern, "\n", sep = "" ) );
    print( fasta.files );
    stopifnot( length( fasta.files ) > 1 );

    seq.lengths <- c();
    seq.counts <- c();
    pairwise.distances.as.matrix.flat.by.fasta.file <- lapply( fasta.files, function ( fasta.file ) {
      in.fasta <- read.dna( paste( fasta.file.prefix.path, fasta.file, sep = "/" ), format = "fasta" );
  
      # Add the consensus.
      .consensus.mat <- matrix( seqinr::consensus( as.character( in.fasta ) ), nrow = 1 );
      consensus <- as.DNAbin( .consensus.mat );
      rownames( consensus ) <- paste( "Consensus" );
      #seq.length <- ncol( consensus );
      # Counting everything for which the consensus is not a gap:
      seq.length <- sum( .consensus.mat[ 1, ] != '-' );
      # Counting everything that is not _entirely_ a gap:
      #.consensus.mat.profile <- seqinr::consensus( as.character( in.fasta ), method = "profile" );
      #seq.length <- sum( apply( .consensus.mat.profile, 2, function( .counts ) { !all( .counts[ -1 ] == 0 ) } ) );
      seq.lengths <<- c( seq.lengths, seq.length );
      seq.counts <<- c( seq.counts, nrow( in.fasta ) );
  
      fasta.with.consensus <- rbind( consensus, in.fasta );
  
      # Remove any columns with a consensus that is a gap, which means
      # that over half of seqs have gaps.  This needs to be removed
      # because it becomes not sensible to consider poisson rates of
      # insertions.  We do however consider rates of deletions, treating
      # them as point mutations (by including the indel counts in the
      # Hamming distance calculation).
      fasta.with.consensus <- fasta.with.consensus[ , .consensus.mat[ 1, ] != "-" ];
      
      # The pairwise.deletion = TRUE argument is necessary so that columns with any gaps are not removed.
      # The optional second call adds in the count of sites with gap differences
      # (gap in one sequence but not the other), which are not included
      # in the first call's "raw" count. Adding them back in, we are
      # effectively treating those as differences.
      fasta.with.consensus.dist <- dist.dna( fasta.with.consensus, model = "N", pairwise.deletion = TRUE );
      fasta.with.consensus.dist[ is.nan( fasta.with.consensus.dist ) ] <- 0;
      if( include.gaps.in.Hamming ) {
          fasta.with.consensus.dist <- fasta.with.consensus.dist + dist.dna( fasta.with.consensus, model = "indel", pairwise.deletion = TRUE );
      }
      
      if( any( is.null( fasta.with.consensus.dist ) ) || any( is.na( fasta.with.consensus.dist ) ) || any( !is.finite( fasta.with.consensus.dist ) ) ) {
        ## TODO: REMOVE
        warning( "UH OH got illegal distance value" );
        print( "UH OH got illegal distance value" );
        print( fasta.with.consensus.dist );
      }
  
      pairwise.distances.as.matrix <- as.matrix( fasta.with.consensus.dist );
  
      pairwise.distances.as.matrix.flat <- matrix( "", nrow = ( ( nrow( pairwise.distances.as.matrix ) * ( ncol( pairwise.distances.as.matrix ) - 1 ) ) / 2 ), ncol = 3 );
      line.i <- 1;
      for( row.i in 1:( nrow( pairwise.distances.as.matrix ) - 1 ) ) {
          for( col.i in ( row.i + 1 ):ncol( pairwise.distances.as.matrix ) ) {
              if( row.i == 1 ) { # consensus, no _1 (multiplicity of the observed sequence)
                  pairwise.distances.as.matrix.flat[ line.i, ] <-
                      c( rownames( pairwise.distances.as.matrix )[ row.i ], paste( colnames( pairwise.distances.as.matrix )[ col.i ], "_1", sep = "" ), pairwise.distances.as.matrix[ row.i, col.i ] );
              } else {
                  pairwise.distances.as.matrix.flat[ line.i, ] <-
                      c( paste( rownames( pairwise.distances.as.matrix )[ row.i ], "_1", sep = "" ), paste( colnames( pairwise.distances.as.matrix )[ col.i ], "_1", sep = "" ), pairwise.distances.as.matrix[ row.i, col.i ] );
              }
              line.i <- line.i + 1;
          }
      }
      return( pairwise.distances.as.matrix.flat );
    } );

    # Ok now merge them.  It suffices to just rbind them.
    pairwise.distances.as.matrix.flat <- do.call( rbind, pairwise.distances.as.matrix.flat.by.fasta.file );
    
    pairwise.distances.as.matrix.file <- paste( output.dir, "/", fasta.file.prefix.short, pairwise.hamming.distances.file.suffix, sep = "" );
    write.table( pairwise.distances.as.matrix.flat, file = pairwise.distances.as.matrix.file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE );

    # Seq lengths may vary.  To make it work with PFitter, which needs just one vlaue here, use the weighted average, where the weight is in terms of evidence for days (which comes only from the intersequence distances).
    .weight.num <- choose( seq.counts, 2 );
    .weight.denom <- sum( .weight.num );
    seq.length <- sum( seq.lengths * ( .weight.num / .weight.denom ) );
    ## TODO: REMOVE
    print( paste( "seq lengths:", paste( seq.lengths, collapse = ", " ), "\n" ) );
    print( paste( "seq weights:", paste( ( .weight.num / .weight.denom ), collapse = ", " ), "\n" ) );
    print( paste( "seq.length:", seq.length, "\n" ) );
    
    R.cmd <- paste( "R CMD BATCH '--vanilla --args", pairwise.distances.as.matrix.file, "2.16e-05", seq.length, "' PFitter.R" );
    return( system( R.cmd ) );
} # runMultiFounderPoissonFitter ( fasta.file.prefix, output.dir )

## Here is where the action is.
fasta.file.prefix <- Sys.getenv( "runMultiFounderPoissonFitter_inputFilenamePrefix" ); # alignment
output.dir <- Sys.getenv( "runMultiFounderPoissonFitter_outputDir" ); # NOTE: will create a subdir for the output, named after the fasta file, with "MultiFounderPoissonFitterDir".
suffix.pattern <- Sys.getenv( "runMultiFounderPoissonFitter_suffixPattern" );
if( suffix.pattern == "" ) {
    suffix.pattern <- NULL;
}
## TODO: REMOVE
# warning( paste( "alignment input file:", fasta.file ) );
# warning( paste( "output dir:", output.dir ) );
# warning( paste( "suffix.pattern", suffix.pattern ) );

if( !is.null( suffix.pattern ) ) {
    print( runMultiFounderPoissonFitter( fasta.file.prefix, output.dir, fasta.file.suffix.pattern = suffix.pattern, output.dir.suffix = "_MultiRegionPoissonFitterDir", pairwise.hamming.distances.file.suffix = "_multiRegionPairwiseHammingDistances.txt" ) );
} else {
    if( file.exists( fasta.file.prefix ) ) {
        print( runMultiFounderPoissonFitter( fasta.file.prefix, output.dir ) );
    } else {
        stop( paste( "File does not exist:", fasta.file ) );
    }
}

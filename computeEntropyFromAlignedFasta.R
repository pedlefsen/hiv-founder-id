library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment", "consensus"
library( "entropy" );

# This computes the consensus of the given alignment, writes it to a fasta file, returns the filename.
# consensus.sequence.name = NA means use the name of the input fasta file (not the full path, just the filename, excluding suffix eg ".fasta"), postpended with "_Consensus".
# include.full.alignment == TRUE means that the output file will contain both the consensus and the input alignment (consensus first).
# if use.sequence.numbers.as.names == TRUE, rename non-consensus output sequences using just their order of appearance (eg 1, 2, 3, etc).  This only applies if include.full.alignment is also TRUE.
computeEntropyFromAlignedFasta <- function ( input.fasta.file, output.dir = NULL, output.file = NULL ) {

    if( length( grep( "^(.*?)\\/[^\\/]+$", input.fasta.file ) ) == 0 ) {
        input.fasta.file.path <- ".";
    } else {
        input.fasta.file.path <-
            gsub( "^(.*?)\\/[^\\/]+$", "\\1", input.fasta.file );
    }
    input.fasta.file.short <-
        gsub( "^.*?\\/?([^\\/]+?)$", "\\1", input.fasta.file, perl = TRUE );
    input.fasta.file.short.nosuffix <-
        gsub( "^([^\\.]+)(\\..+)?$", "\\1", input.fasta.file.short, perl = TRUE );
    input.fasta.file.suffix <-
        gsub( "^([^\\.]+)(\\..+)?$", "\\2", input.fasta.file.short, perl = TRUE );

    if( !is.null( output.file ) ) {
      if( length( grep( "^(.*?)\\/[^\\/]+$", output.file ) ) == 0 ) {
          output.file.path <- NULL;
          output.file.path.is.absolute <- NA;
      } else {
          output.file.path <-
              gsub( "^(.*?)\\/[^\\/]+$", "\\1", output.file );
          output.file.path.is.absolute <- ( substring( output.file.path, 1, 1 ) == "/" );
      }
      output.file.short <-
          gsub( "^.*?\\/?([^\\/]+?)$", "\\1", output.file, perl = TRUE );

      if( !is.null( output.file.path ) && output.file.path.is.absolute ) {
          output.dir <- output.file.path;
      } else if( is.null( output.dir ) ) {
        if( is.null( output.file.path ) ) {
            output.dir <- input.fasta.file.path;
        } else {
            output.dir <- output.file.path;
        }
      } else {
          output.dir <- paste( output.dir, output.file.path, sep = "/" );
      }
      output.file <- output.file.short;
    } else { # is.null( output.file )
        output.file <- paste( input.fasta.file.short.nosuffix, ".entropy.txt", sep = "" );
    }
    if( is.null( output.dir ) || ( nchar( output.dir ) == 0 ) ) {
      output.dir <- ".";
    }
    ## Remove "/" from end of output.dir
    output.dir <-
      gsub( "^(.*?)\\/+$", "\\1", output.dir );

    input.fasta <- read.dna( input.fasta.file, format = "fasta" );
    # IUPAC profile (may contain ambiguity chars)
    input.fasta.profile.iupac <- seqinr::consensus( as.character( input.fasta ), method = "profile" );
    input.fasta.profile <- input.fasta.profile.iupac[ intersect( c( "-", "a", "c", "g", "t" ), rownames( input.fasta.profile.iupac ) ), ];
    .iupacs <- setdiff( rownames( input.fasta.profile.iupac ), c( "-", "a", "c", "g", "t" ) );
    # DIVVY up ambiguous counts evenly among the component bases.
    for( .iupac in .iupacs ) {
        .component.bases <- amb( .iupac );
        input.fasta.profile[ .component.bases, ] <-
            input.fasta.profile[ .component.bases, ] + ( input.fasta.profile.iupac[ .iupac, ] / length( .component.bases) );
    }
    input.fasta.profile.nogap <- input.fasta.profile[ setdiff( rownames( input.fasta.profile ), "-" ), ];

    entropies <- apply( input.fasta.profile.nogap, 2, entropy.empirical, unit = "log2" );
    entropies.sd <- sd( entropies, na.rm = T );
    num.entropies <- length( !is.na( entropies ) );
    .results <- summary( entropies[ !is.na( entropies ) ] );
    .results.string <- sprintf( "%0.4f", .results );
    names( .results.string ) <- names( .results );
    .results <- c( N = as.character( nrow( input.fasta ) ), K = as.character( num.entropies ), .results.string, SD = round( entropies.sd, digits = 4 ) );
    .results[ .results == "-0.0000" ] <- "0.0000";
    
    output.file.path <-
        paste( output.dir, "/", output.file, sep = "" );
    write.table( t( as.matrix( .results ) ), file =output.file.path, row.names = F, sep = "\t" );

    # Return the file name.
    return( output.file.path );
} # computeEntropyFromAlignedFasta ( input.fasta.file, ... )

## Here is where the action is.
input.fasta.file <- Sys.getenv( "computeEntropyFromAlignedFasta_inputFilename" );
output.fasta.file <- Sys.getenv( "computeEntropyFromAlignedFasta_outputFilename" );
if( nchar( output.fasta.file ) == 0 ) {
    output.fasta.file <- NULL;
}
output.dir <- Sys.getenv( "computeEntropyFromAlignedFasta_outputDir" );
if( nchar( output.dir ) == 0 ) {
    output.dir <- NULL;
}

## TODO: REMOVE
# warning( paste( "aligned fasta input file:", input.fasta.file ) );
# if( !is.null( output.dir ) ) {
#     warning( paste( "consensus fasta output dir:", output.dir ) );
# }
# if( !is.null( output.fasta.file ) ) {
#     warning( paste( "consensus fasta output file:", output.fasta.file ) );
# }
if( file.exists( input.fasta.file ) ) {
    print( computeEntropyFromAlignedFasta( input.fasta.file, output.dir = output.dir, output.file = output.fasta.file ) );
} else {
    stop( paste( "File does not exist:", input.fasta.file ) );
}

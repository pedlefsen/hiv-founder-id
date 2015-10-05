library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment", "consensus"

                                        # This computes the consensus of the given alignment, writes it to a fasta file, returns the filename.
# consensus.sequence.name = NA means use the name of the input fasta file (not the full path, just the filename).
computeConsensusSequenceFromAlignedFasta <- function ( input.fasta.file, output.dir = NULL, output.file = NULL, consensus.sequence.name = NA, output.fasta.width = 72 ) {

    if( length( grep( "^(.*?)\\/[^\\/]+$", input.fasta.file ) ) == 0 ) {
        input.fasta.file.path <- ".";
    } else {
        input.fasta.file.path <-
            gsub( "^(.*?)\\/[^\\/]+$", "\\1", input.fasta.file );
    }
    input.fasta.file.short <-
        gsub( "^.*?\\/?([^\\/]+?)$", "\\1", input.fasta.file, perl = TRUE );

    if( is.na( consensus.sequence.name ) ) {
        consensus.sequence.name <- input.fasta.file.short;
    }
    
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
      output.file <- paste( input.fasta.file.short, ".cons.fasta", sep = "" );
    }
    if( is.null( output.dir ) || ( nchar( output.dir ) == 0 ) ) {
      output.dir <- ".";
    }
    ## Remove "/" from end of output.dir
    output.dir <-
      gsub( "^(.*?)\\/+$", "\\1", output.dir );

    input.fasta <- read.dna( input.fasta.file, format = "fasta" );

    consensus <- as.DNAbin( matrix( seqinr::consensus( as.character( input.fasta ) ), nrow = 1 ) );
    rownames( consensus ) <- consensus.sequence.name;
    
    output.fasta.path <-
        paste( output.dir, "/", output.file, sep = "" );
    write.dna( consensus, output.fasta.path, format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = output.fasta.width );

    # Return the file name.
    return( output.fasta.path );
} # computeConsensusSequenceFromAlignedFasta ( input.fasta.file, ... )

## Here is where the action is.
input.fasta.file <- Sys.getenv( "computeConsensusSequenceFromAlignedFasta_inputFilename" );
output.fasta.file <- Sys.getenv( "computeConsensusSequenceFromAlignedFasta_outputFilename" );
if( nchar( output.fasta.file ) == 0 ) {
    output.fasta.file <- NULL;
}
output.dir <- Sys.getenv( "computeConsensusSequenceFromAlignedFasta_outputDir" );
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
    print( computeConsensusSequenceFromAlignedFasta( input.fasta.file, output.dir = output.dir, output.file = output.fasta.file ) );
} else {
    stop( paste( "File does not exist:", input.fasta.file ) );
}

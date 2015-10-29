library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment", "consensus"

# This computes the consensus of the given alignment, writes it to a fasta file, returns the filename.
# consensus.sequence.name = NA means use the name of the input fasta file (not the full path, just the filename, excluding suffix eg ".fasta"), postpended with "_Consensus".
# include.full.alignment == TRUE means that the output file will contain both the consensus and the input alignment (consensus first).
# if use.sequence.numbers.as.names == TRUE, rename non-consensus output sequences using just their order of appearance (eg 1, 2, 3, etc).  This only applies if include.full.alignment is also TRUE.
computeConsensusSequenceFromAlignedFasta <- function ( input.fasta.file, output.dir = NULL, output.file = NULL, consensus.sequence.name = NA, output.fasta.width = 72, include.full.alignment = FALSE, include.consensus = TRUE, use.sequence.numbers.as.names = FALSE ) {

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

    if( is.na( consensus.sequence.name ) ) {
        consensus.sequence.name <- paste( input.fasta.file.short.nosuffix, "Consensus", sep = "_" );
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
        if( include.full.alignment ) {
            output.file <- paste( input.fasta.file.short.nosuffix, ".withcons", input.fasta.file.suffix, sep = "" );
        } else {
            output.file <- paste( input.fasta.file.short.nosuffix, ".cons", input.fasta.file.suffix, sep = "" );
        }
    }
    if( is.null( output.dir ) || ( nchar( output.dir ) == 0 ) ) {
      output.dir <- ".";
    }
    ## Remove "/" from end of output.dir
    output.dir <-
      gsub( "^(.*?)\\/+$", "\\1", output.dir );

    input.fasta <- read.dna( input.fasta.file, format = "fasta" );

          if( use.sequence.numbers.as.names ) {
              rownames( input.fasta ) <- 1:nrow( input.fasta );
          }

    if( include.consensus ) {
      consensus <- as.DNAbin( matrix( seqinr::consensus( as.character( input.fasta ) ), nrow = 1 ) );
      rownames( consensus ) <- consensus.sequence.name;
      if( include.full.alignment ) {
        consensus <- rbind( consensus, input.fasta );
      }
    } else {
      consensus <- input.fasta;
    }

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
include.full.alignment <- Sys.getenv( "computeConsensusSequenceFromAlignedFasta_includeFullAlignment" );
include.consensus <- Sys.getenv( "computeConsensusSequenceFromAlignedFasta_includeConsensus" );
if( ( nchar( include.full.alignment ) == 0 ) || ( include.full.alignment == "0" ) || ( toupper( include.full.alignment ) == "F" ) || ( toupper( include.full.alignment ) == "FALSE" ) ) {
    include.full.alignment <- FALSE;
} else {
    include.full.alignment <- TRUE;
}
use.sequence.numbers.as.names <- Sys.getenv( "computeConsensusSequenceFromAlignedFasta_useSeqeunceNumbersAsNames" );
if( ( nchar( use.sequence.numbers.as.names ) == 0 ) || ( use.sequence.numbers.as.names == "0" ) || ( toupper( use.sequence.numbers.as.names ) == "F" ) || ( toupper( use.sequence.numbers.as.names ) == "FALSE" ) ) {
    use.sequence.numbers.as.names <- FALSE;
} else {
    use.sequence.numbers.as.names <- TRUE;
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
    print( computeConsensusSequenceFromAlignedFasta( input.fasta.file, output.dir = output.dir, output.file = output.fasta.file, include.full.alignment = include.full.alignment, include.consensus = include.consensus,  use.sequence.numbers.as.names = use.sequence.numbers.as.names ) );
} else {
    stop( paste( "File does not exist:", input.fasta.file ) );
}

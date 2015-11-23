library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment", "translate"

codonIndexToAAIndex <- function ( codon.index ) {
    floor( ( codon.index - 1 ) / 3 ) + 1
}
AAIndexToCodonIndices <- Vectorize( function ( aa.index ) {
    ( ( aa.index - 1 ) * 3 ) + 1:3
} );

# if the input.aa.fasta.file is not given, it will be computed by translating the input.codon.fasta.file.
# if the input.reference.aa.fasta.file is not given, it will be computed from the input.reference.codon.fasta.file, or, if that is not given, from the dna consensus of the input.codon.fasta.file.
maskSynonymousCodonsInAlignedFasta <- function ( input.codon.fasta.file, input.aa.fasta.file = NULL, input.reference.codon.fasta.file = NULL, input.reference.aa.fasta.file = NULL, output.dir = NULL, output.fasta.file = NULL, output.fasta.width = 72, mask.char = "-", mask.nonsynonymous = FALSE ) {
    if( length( grep( "^(.*?)\\/[^\\/]+$", input.codon.fasta.file ) ) == 0 ) {
        input.codon.fasta.file.path <- ".";
    } else {
        input.codon.fasta.file.path <-
            gsub( "^(.*?)\\/[^\\/]+$", "\\1", input.codon.fasta.file );
    }
    input.codon.fasta.file.short <-
        gsub( "^.*?\\/?([^\\/]+?)$", "\\1", input.codon.fasta.file, perl = TRUE );
    input.codon.fasta.file.short.nosuffix <-
        gsub( "^([^\\.]+)(\\..+)?$", "\\1", input.codon.fasta.file.short, perl = TRUE );
    input.codon.fasta.file.suffix <-
        gsub( "^([^\\.]+)(\\..+)?$", "\\2", input.codon.fasta.file.short, perl = TRUE );

    if( !is.null( output.fasta.file ) ) {
      if( length( grep( "^(.*?)\\/[^\\/]+$", output.fasta.file ) ) == 0 ) {
          output.fasta.file.path <- NULL;
          output.fasta.file.path.is.absolute <- NA;
      } else {
          output.fasta.file.path <-
              gsub( "^(.*?)\\/[^\\/]+$", "\\1", output.fasta.file );
          output.fasta.file.path.is.absolute <- ( substring( output.fasta.file.path, 1, 1 ) == "/" );
      }
      output.fasta.file.short <-
          gsub( "^.*?\\/?([^\\/]+?)$", "\\1", output.fasta.file, perl = TRUE );

      if( !is.null( output.fasta.file.path ) && output.fasta.file.path.is.absolute ) {
          output.dir <- output.fasta.file.path;
      } else if( is.null( output.dir ) ) {
        if( is.null( output.fasta.file.path ) ) {
            output.dir <- input.codon.fasta.file.path;
        } else {
            output.dir <- output.fasta.file.path;
        }
      } else {
          output.dir <-
              paste( output.dir, output.fasta.file.path, sep = "/" );
      }
      output.fasta.file <- output.fasta.file.short;
    } else { # is.null( output.fasta.file )
        if( mask.nonsynonymous ) {
            output.fasta.file <-
                paste( input.codon.fasta.file.short.nosuffix, "_maskNonsynonymousCodons", input.codon.fasta.file.suffix, sep = "" );
        } else {
            output.fasta.file <-
                paste( input.codon.fasta.file.short.nosuffix, "_maskSynonymousCodons", input.codon.fasta.file.suffix, sep = "" );
        }
    }
    if( is.null( output.dir ) || ( nchar( output.dir ) == 0 ) ) {
      output.dir <- ".";
    }
    ## Remove "/" from end of output.dir
    output.dir <-
      gsub( "^(.*?)\\/+$", "\\1", output.dir );

    input.fasta <- read.dna( input.codon.fasta.file, format = "fasta" );
    
    if( is.null( input.aa.fasta.file ) ) {
        input.aa.fasta.mat <- t( apply( as.matrix( 1:nrow( input.fasta ) ), 1, function ( i ) {
            seqinr::translate( as.character( input.fasta[ i, ] ), NAstring = "-" )
        } ) );
    } else {
        input.aa.fasta.in <-
            readAAMultipleAlignment( input.aa.fasta.file, "fasta" );
        input.aa.fasta.mat <- t( apply( as.matrix( 1:nrow( input.aa.fasta.in ) ), 1, function ( i ) {
            s2c( as.character( unmasked( input.aa.fasta.in )[ i, ] ) )
        } ) );
    }

    if( is.null( input.reference.aa.fasta.file ) ) {
        if( is.null( input.reference.codon.fasta.file ) ) {
            input.reference.fasta <-
                as.DNAbin( matrix( seqinr::consensus( as.character( input.fasta ) ), nrow = 1 ) );
            rownames( input.reference.fasta ) <- "Consensus";
        } else {
            input.reference.fasta <-
                read.dna( input.reference.codon.fasta.file, format = "fasta" );
        }
        input.reference.aa.fasta.mat <- t( apply( as.matrix( 1:nrow( input.reference.fasta ) ), 1, function ( i ) {
            seqinr::translate( as.character( input.reference.fasta[ i, ] ), NAstring = "-" )
        } ) );
    } else {
        input.reference.aa.fasta.in <-
            readAAMultipleAlignment( input.reference.aa.fasta.file, "fasta" );
        input.reference.aa.fasta.mat <- t( apply( as.matrix( 1:nrow( input.reference.aa.fasta.in ) ), 1, function ( i ) {
            s2c( as.character( unmasked( input.reference.aa.fasta.in )[ i, ] ) )
        } ) );
    }
    stopifnot( length( input.reference.aa.fasta.mat ) == ncol( input.aa.fasta.mat ) );

    # Get the indices of locations that differ between AAs and ref AA.
    output.fasta.mat <- as.character( input.fasta );
    for( seq.i in 1:nrow( input.fasta ) ) {
        if( mask.nonsynonymous ) {
          differing.aa.indices <-
              which( input.aa.fasta.mat[ seq.i, ] != input.reference.aa.fasta.mat[ 1, ] );
          if( length( differing.aa.indices ) == 0 ) {
              next;
          }
          differing.dna.indices <-
              c( AAIndexToCodonIndices( differing.aa.indices ) );
          ## TODO: REMOVE
          #cat( "MASKING:", differing.dna.indices, fill = T );
          output.fasta.mat[ seq.i, differing.dna.indices ] <- mask.char;
        } else {
          matching.aa.indices <-
              which( input.aa.fasta.mat[ seq.i, ] == input.reference.aa.fasta.mat[ 1, ] );
          if( length( matching.aa.indices ) == 0 ) {
              next;
          }
          matching.dna.indices <-
              c( AAIndexToCodonIndices( matching.aa.indices ) );
          ## TODO: REMOVE
          #cat( "MASKING:", matching.dna.indices, fill = T );
          output.fasta.mat[ seq.i, matching.dna.indices ] <- mask.char;
        }
    }
    output.fasta <- as.DNAbin( output.fasta.mat );
    
    output.fasta.file.path <-
        paste( output.dir, "/", output.fasta.file, sep = "" );
    write.dna( output.fasta, output.fasta.file.path, format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = output.fasta.width );

    # Return the file name.
    return( output.fasta.file.path );
} # maskSynonymousCodonsInAlignedFasta ( input.codon.fasta.file, ... )

## Here is where the action is.
input.codon.fasta.file <- Sys.getenv( "maskSynonymousCodonsInAlignedFasta_inputFilename" );
input.aa.fasta.file <- Sys.getenv( "maskSynonymousCodonsInAlignedFasta_translatedInputFilename" );
output.fasta.file <- Sys.getenv( "maskSynonymousCodonsInAlignedFasta_outputFilename" );
if( nchar( output.fasta.file ) == 0 ) {
    output.fasta.file <- NULL;
}
output.dir <- Sys.getenv( "maskSynonymousCodonsInAlignedFasta_outputDir" );
if( nchar( output.dir ) == 0 ) {
    output.dir <- NULL;
}

## TODO: REMOVE
# warning( paste( "aligned fasta input file:", input.codon.fasta.file ) );
# if( !is.null( output.dir ) ) {
#     warning( paste( "consensus fasta output dir:", output.dir ) );
# }
# if( !is.null( output.fasta.file ) ) {
#     warning( paste( "consensus fasta output file:", output.fasta.file ) );
# }
if( file.exists( input.codon.fasta.file ) ) {
    print( maskSynonymousCodonsInAlignedFasta( input.codon.fasta.file, input.aa.fasta.file, output.dir = output.dir, output.fasta.file = output.fasta.file ) );
} else {
    stop( paste( "File does not exist:", input.codon.fasta.file ) );
}

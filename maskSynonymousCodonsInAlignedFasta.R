library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment", "translate"

maskSynonymousCodonsInAlignedFasta <- function ( input.codon.fasta.file, input.aa.fasta.file = NULL, output.dir = NULL, output.fasta.file = NULL ) {
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
          output.dir <- paste( output.dir, output.fasta.file.path, sep = "/" );
      }
      output.fasta.file <- output.fasta.file.short;
    } else { # is.null( output.fasta.file )
        output.fasta.file <- paste( input.codon.fasta.file.short.nosuffix, ".entropy.txt", sep = "" );
    }
    if( is.null( output.dir ) || ( nchar( output.dir ) == 0 ) ) {
      output.dir <- ".";
    }
    ## Remove "/" from end of output.dir
    output.dir <-
      gsub( "^(.*?)\\/+$", "\\1", output.dir );

    input.fasta <- read.dna( input.codon.fasta.file, format = "fasta" );
    #input.fasta.tr1 <- apply( input.fasta, 1,  );
    #input.aa.fasta <- readAAMultipleAlignment( input.aa.fasta.file, "fasta" );
    .whichs <- which(
        s2c( as.character( unmasked( input.aa.fasta )[ i, ] ) ) !=
        seqinr::translate( as.character( input.fasta[1,]), NAstring = "-" )
    );
    .whichs.codons <-
        sapply( .whichs, function( .which ) { ( ( .which - 1 ) * 3 ) + 1:3 } );
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
    
    output.fasta.file.path <-
        paste( output.dir, "/", output.fasta.file, sep = "" );
    write.table( t( as.matrix( .results ) ), file =output.fasta.file.path, row.names = F, sep = "\t" );

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

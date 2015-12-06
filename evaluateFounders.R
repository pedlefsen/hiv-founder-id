library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment"
library( "Biostrings" ) # for "pairwiseAlignment"

# This compares two nucleotide fasta files, each containing one or more sequences (aligned or not, it doesn't matter; gaps will be stripped internally).  The comparison is conducted in both nucleotide and amino acid space after gene-cutting, codon-aligning, and translating the sequences using GeneCutter at LANL (see runGeneCutterOnline.pl).
evaluateFounders <- function ( estimates.fasta.file, standards.fasta.file, output.dir = NULL, output.file = NULL, output.fasta.width = 72 ) {

    if( length( grep( "^(.*?)\\/[^\\/]+$", estimates.fasta.file ) ) == 0 ) {
        estimates.fasta.file.path <- ".";
    } else {
        estimates.fasta.file.path <-
            gsub( "^(.*?)\\/[^\\/]+$", "\\1", estimates.fasta.file );
    }
    estimates.fasta.file.short <-
        gsub( "^.*?\\/?([^\\/]+?)$", "\\1", estimates.fasta.file, perl = TRUE );
    estimates.fasta.file.short.nosuffix <-
        gsub( "^([^\\.]+)(\\..+)?$", "\\1", estimates.fasta.file.short, perl = TRUE );
    estimates.fasta.file.suffix <-
        gsub( "^([^\\.]+)(\\..+)?$", "\\2", estimates.fasta.file.short, perl = TRUE );

    if( length( grep( "^(.*?)\\/[^\\/]+$", truths.fasta.file ) ) == 0 ) {
        truths.fasta.file.path <- ".";
    } else {
        truths.fasta.file.path <-
            gsub( "^(.*?)\\/[^\\/]+$", "\\1", truths.fasta.file );
    }
    truths.fasta.file.short <-
        gsub( "^.*?\\/?([^\\/]+?)$", "\\1", truths.fasta.file, perl = TRUE );
    truths.fasta.file.short.nosuffix <-
        gsub( "^([^\\.]+)(\\..+)?$", "\\1", truths.fasta.file.short, perl = TRUE );
    truths.fasta.file.suffix <-
        gsub( "^([^\\.]+)(\\..+)?$", "\\2", truths.fasta.file.short, perl = TRUE );

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
            output.dir <- estimates.fasta.file.path;
        } else {
            output.dir <- output.file.path;
        }
      } else {
          output.dir <- paste( output.dir, output.file.path, sep = "/" );
      }
      output.file <- output.file.short;
    } else { # is.null( output.file )
      output.file <- paste( estimates.fasta.file.short.nosuffix, "_evaluateFounders.tab", sep = "" );
    }
    if( is.null( output.dir ) || ( nchar( output.dir ) == 0 ) ) {
      output.dir <- ".";
    }
    ## Remove "/" from end of output.dir
    output.dir <-
      gsub( "^(.*?)\\/+$", "\\1", output.dir );

    estimates.fasta <- read.dna( estimates.fasta.file, format = "fasta" );

    ## Make an ungapped version, if it doesn't already exist.
    estimates.fasta.file.ungapped <- paste( estimates.fasta.file.path, "/", estimates.fasta.file.short.nosuffix, "_ungapped", estimates.fasta.file.suffix, sep = "" );
    if( !file.exists( estimates.fasta.file.ungapped ) ) {
        # Create an ungapped version, and save it.
        .estimates.fasta.as.character <- as.character( estimates.fasta );
        .result.ignored <- lapply( 1:nrow( .estimates.fasta.as.character ), function( .row.i ) {
            .row.seq.chars <- .estimates.fasta.as.character[ .row.i, ];
            .row.seq.ungapped <- as.DNAbin( matrix( .row[ .row != "-" ], nrow = 1 ) );
            rownames( .row.seq.ungapped ) <- rownames( estimates.fasta )[ .row.i ];
            write.dna( .row.seq.ungapped, file = estimates.fasta.file.ungapped, format = "fasta", append = TRUE, colsep = "", indent = "", blocksep = 0, nbcol = 1, colw = output.fasta.width );
            return( NULL );
        } );
    }

    truths.fasta <- read.dna( truths.fasta.file, format = "fasta" );

    ## Make an ungapped version, if it doesn't already exist.
    truths.fasta.file.ungapped <- paste( truths.fasta.file.path, "/", truths.fasta.file.short.nosuffix, "_ungapped", truths.fasta.file.suffix, sep = "" );
    if( !file.exists( truths.fasta.file.ungapped ) ) {
        # Create an ungapped version, and save it.
        .truths.fasta.as.character <- as.character( truths.fasta );
        .result.ignored <- lapply( 1:nrow( .truths.fasta.as.character ), function( .row.i ) {
            .row.seq.chars <- .truths.fasta.as.character[ .row.i, ];
            .row.seq.ungapped <- as.DNAbin( matrix( .row[ .row != "-" ], nrow = 1 ) );
            rownames( .row.seq.ungapped ) <- rownames( truths.fasta )[ .row.i ];
            write.dna( .row.seq.ungapped, file = truths.fasta.file.ungapped, format = "fasta", append = TRUE, colsep = "", indent = "", blocksep = 0, nbcol = 1, colw = output.fasta.width );
            return( NULL );
        } );
    }
    
    ## Put together the two fasta files, the lazy way.
    combined.ungapped.fasta.file.nosuffix <- paste( output.dir, "/", estimates.fasta.file.short.nosuffix, "_ungapped_with_", truths.fasta.file.short.nosuffix, "_ungapped_combined", sep = "" );
    combined.ungapped.fasta.file <- paste( combined.ungapped.fasta.file.nosuffix, truths.fasta.file.suffix, sep = "" );
    system( paste( "cp", truths.fasta.file.ungapped, combined.ungapped.fasta.file ) );
    system( paste( "cat", estimates.fasta.file.ungapped, ">>", combined.ungapped.fasta.file ) );
    ## Run it through GeneCutter.
    ## TODO: Add other comparison seqs first, eg CON_M
    system( paste( "perl runGeneCutterOnline.pl -V ", combined.ungapped.fasta.file, output.dir ) );
    
    ## ERE I AM.
    nucleotides.zipfile <- paste( combined.ungapped.fasta.file.nosuffix, "_allnucs.zip", sep = "" );
    proteins.zipfile <- paste( combined.ungapped.fasta.file.nosuffix, "_allproteins.zip", sep = "" );
    stopifnot( file.exists( nucleotides.zipfile ) );
    stopifnot( file.exists( proteins.zipfile ) );

    ## Now unzip them.
    nucleotides.dir <- paste( combined.ungapped.fasta.file.nosuffix, "_allnucs", sep = "" );
    proteins.dir <- paste( combined.ungapped.fasta.file.nosuffix, "_allproteins", sep = "" );
    system( paste( "rm -rf", nucleotides.dir ) );
    system( paste( "rm -rf", proteins.dir ) );
    system( paste( "unzip", nucleotides.zipfile, "-d", nucleotides.dir, sep = " " ) )
    system( paste( "unzip", proteins.zipfile, "-d", proteins.dir, sep = " " ) )

    # Now process each of the files
    ## ERE I AM
    
    output.table.path <-
        paste( output.dir, "/", output.file, sep = "" );

    # Return the file name.
    return( output.table.path );
} # computeConsensusSequenceFromAlignedFasta ( estimates.fasta.file, ... )

## Here is where the action is.
estimates.fasta.file <- Sys.getenv( "computeConsensusSequenceFromAlignedFasta_estimatesFilename" );
truths.fasta.file <- Sys.getenv( "computeConsensusSequenceFromAlignedFasta_truthsFilename" );
output.table.file <- Sys.getenv( "computeConsensusSequenceFromAlignedFasta_outputFilename" );
if( nchar( output.table.file ) == 0 ) {
    output.table.file <- NULL;
}
output.dir <- Sys.getenv( "computeConsensusSequenceFromAlignedFasta_outputDir" );
if( nchar( output.dir ) == 0 ) {
    output.dir <- NULL;
}

## TODO: REMOVE
# warning( paste( "estimates fasta file:", estimates.fasta.file ) );
# warning( paste( "truths fasta file:", truths.fasta.file ) );
# if( !is.null( output.dir ) ) {
#     warning( paste( "consensus fasta output dir:", output.dir ) );
# }
# if( !is.null( output.table.file ) ) {
#     warning( paste( "consensus fasta output file:", output.table.file ) );
# }
if( file.exists( estimates.fasta.file ) ) {
    if( file.exists( truths.fasta.file ) ) {
        print( computeConsensusSequenceFromAlignedFasta( estimates.fasta.file, truths.fasta.file, output.dir = output.dir, output.file = output.table.file ) );
    } else {
        stop( paste( "'truths' fasta file does not exist:", truths.fasta.file ) );
    }
} else {
    stop( paste( "'estimates' fasta file does not exist:", estimates.fasta.file ) );
}

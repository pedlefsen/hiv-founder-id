library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment", "consensus"

# This removes duplicate sequences from an input multiple alignment in fasta format, writes the output to another fasta file, writes a table of duplicate sequences, and returns the output filename.
# The output table is tab-delimited with two columns: retained and removed, where retained is a single sequence name (the one not removed) and removed is a comma-separated list of removed sequence names.
removeDuplicateSequencesFromAlignedFasta <- function ( input.fasta.file, output.dir = NULL, output.fasta.file = NULL, output.table.file = NULL, output.fasta.width = 72 ) {

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
            output.dir <- input.fasta.file.path;
        } else {
            output.dir <- output.fasta.file.path;
        }
      } else {
          output.dir <- paste( output.dir, output.fasta.file.path, sep = "/" );
      }
      output.fasta.file <- output.fasta.file.short;
    } else { # is.null( output.fasta.file )
      output.fasta.file <- paste( input.fasta.file.short.nosuffix, "_removeDuplicateSequences", input.fasta.file.suffix, sep = "" );
    }
    if( is.null( output.dir ) || ( nchar( output.dir ) == 0 ) ) {
      output.dir <- ".";
    }
    ## Remove "/" from end of output.dir
    output.dir <-
      gsub( "^(.*?)\\/+$", "\\1", output.dir );

    if( is.null( output.table.file ) ) {
        output.table.file <- paste( input.fasta.file.short.nosuffix, "_removeDuplicateSequences.tbl", sep = "" );
    }
      if( length( grep( "^(.*?)\\/[^\\/]+$", output.table.file ) ) == 0 ) {
          output.table.file.path <- output.dir;
          output.table.file.path.is.absolute <- NA;
      } else {
          output.table.file.path <-
              gsub( "^(.*?)\\/[^\\/]+$", "\\1", output.table.file );
          output.table.file.path.is.absolute <- ( substring( output.table.file.path, 1, 1 ) == "/" );
      }
      output.table.file.short <-
          gsub( "^.*?\\/?([^\\/]+?)$", "\\1", output.table.file, perl = TRUE );
      output.table.file <- paste( output.table.file.path, output.table.file.short, sep = "/" );
    
    input.fasta <- read.dna( input.fasta.file, format = "fasta" );

    # identical sequences have no distance between them.
    # The pairwise.deletion = TRUE argument is necessary so that columns with any gaps are not removed.
    # The second call adds in the count of sites with gap differences
    # (gap in one sequence but not the other), which are not included
    # in the first call's "raw" count. Adding them back in, we are
    # effectively treating those as differences.

    in.dist <- dist.dna( input.fasta, model = "N", pairwise.deletion = TRUE );
    in.dist[ is.nan( in.dist ) ] <- 0;
    in.dist <- in.dist + dist.dna( input.fasta, model = "indel", pairwise.deletion = TRUE ); 
    if( any( is.null( in.dist ) ) || any( is.na( in.dist ) ) || any( !is.finite( in.dist ) ) ) {
      ## TODO: REMOVE
      warning( "UH OH got illegal distance value" );
      print( "UH OH got illegal distance value" );
      print( in.dist );
    }
    dist.mat <- as.matrix( in.dist );
    
    # Ok so duplicates are those pairs with dist.mat == 0.
    # Keep the first.
    seq.was.already.removed <- rep( FALSE, nrow( dist.mat ) );
    output.table.as.list <- list();
    for( seq.i in ( 1:nrow( dist.mat ) ) ) {
        if( seq.was.already.removed[ seq.i ] ) {
            next;
        }
        # We should not find this one to be identical to one already
        # removed, since then this one would also have already been
        # removed.
        stopifnot( !any( dist.mat[ seq.i, seq.was.already.removed ] == 0 ) );
        # Are there any identical to this one?
        if( !any( dist.mat[ seq.i, -seq.i ] == 0 ) ) {
            next;
        }
        output.table.as.list[[ rownames( dist.mat )[ seq.i ] ]] <-
            paste( colnames( dist.mat )[ -seq.i ][ which( dist.mat[ seq.i, -seq.i ] == 0 ) ], collapse = "," );
        seq.was.already.removed[ -seq.i ][ which( dist.mat[ seq.i, -seq.i ] == 0 ) ] <- TRUE;
    } # End foreach seq.i

    if( !is.null( output.table.file ) && ( length( output.table.as.list ) > 0 ) ) {
        output.table <- do.call( rbind, output.table.as.list );
        output.table <- cbind( rownames( output.table ), output.table );
        colnames( output.table ) <- c( "retained", "removed" );
        rownames( output.table ) <- NULL;
        write.table( output.table, file = output.table.file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE );
    }
    
    output.fasta <- input.fasta[ !seq.was.already.removed, , drop = FALSE ];
               
    output.fasta.path <-
        paste( output.dir, "/", output.fasta.file, sep = "" );
    write.dna( output.fasta, output.fasta.path, format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = output.fasta.width );

    # Return the file name.
    return( output.fasta.path );
} # removeDuplicateSequencesFromAlignedFasta ( input.fasta.file, ... )

## Here is where the action is.
input.fasta.file <- Sys.getenv( "removeDuplicateSequencesFromAlignedFasta_inputFilename" );
output.fasta.file <- Sys.getenv( "removeDuplicateSequencesFromAlignedFasta_outputFastaFilename" );
if( nchar( output.fasta.file ) == 0 ) {
    output.fasta.file <- NULL;
}
output.dir <- Sys.getenv( "removeDuplicateSequencesFromAlignedFasta_outputDir" );
if( nchar( output.dir ) == 0 ) {
    output.dir <- NULL;
}
output.table.file <- Sys.getenv( "removeDuplicateSequencesFromAlignedFasta_outputTableFilename" );
if( nchar( output.table.file ) == 0 ) {
    output.table.file <- NULL;
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
    print( removeDuplicateSequencesFromAlignedFasta( input.fasta.file, output.dir = output.dir, output.fasta.file = output.fasta.file, output.table.file = output.table.file ) );
} else {
    stop( paste( "File does not exist:", input.fasta.file ) );
}

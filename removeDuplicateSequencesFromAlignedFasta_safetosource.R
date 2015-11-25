library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment", "consensus"

# This removes duplicate sequences from an input multiple alignment in fasta format, writes the output to another fasta file, writes a table of duplicate sequences, and returns the output filename.
# The output table is tab-delimited with two columns: retained and removed, where retained is a single sequence name (the one not removed) and removed is a comma-separated list of removed sequence names.
removeDuplicateSequencesFromAlignedFasta <- function ( input.fasta.file, output.dir = NULL, output.fasta.file = NULL, output.table.file = NULL, output.fasta.width = 72, add.copy.number.to.sequence.names = FALSE ) {

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
    ## NOTE: Because of ambiguity codes, the transitivity property
    ## might _not_ hold, eg you could have dist(A,B) == 0, dist(B,C) == 0, dist(A, C) == 1.
    
    # Ok so duplicates are those pairs with dist.mat == 0.
    # Keep the first.
    seq.removed.since.represented.by.seq <- rep( 0, nrow( dist.mat ) );
    names( seq.removed.since.represented.by.seq ) <- colnames( dist.mat );
    output.table.as.list <- list();
    for( seq.i in ( 1:nrow( dist.mat ) ) ) {
        if( seq.removed.since.represented.by.seq[ seq.i ] != 0 ) {
            next;
        }
        # If the distances were true metric distances, ie if the transitivity property did hold, then it would be the case that we should not find this one to be identical to one already
        # removed, since then this one would also have already been
        # removed.  However, it does not hold (see above).  Thus we might need to add this to a collection of already removed sequences if it is zero-distance from something already removed.
        if( any( seq.removed.since.represented.by.seq != 0 ) && any( dist.mat[ seq.i, ( seq.removed.since.represented.by.seq > 0 ) ] == 0 ) ) {
            # This is part of a group that was already removed. It got
            # missed by being not perfectly identical to every other
            # member due to ambiguities.  We go ahead and include it if it is distance zero from any member.
            .representative.seq <- seq.removed.since.represented.by.seq[ ( seq.removed.since.represented.by.seq > 0 ) & ( dist.mat[ seq.i, ] == 0 ) ];
            stopifnot( length( .representative.seq ) > 0 );
            if( length( .representative.seq ) > 1 ) {
                stopifnot( all( .representative.seq == .representative.seq[ 1 ] ) );
                .representative.seq <- .representative.seq[ 1 ];
            }
            ## TODO: REMOVE
            # print( paste( "adding", colnames( dist.mat )[ seq.i ], "to", rownames( dist.mat )[ .representative.seq ] ) );
            output.table.as.list[[ rownames( dist.mat )[ .representative.seq ] ]][ 2 ] <-
                paste( output.table.as.list[[ rownames( dist.mat )[ .representative.seq ] ]][ 2 ], colnames( dist.mat )[ seq.i ], sep = "," );
            seq.removed.since.represented.by.seq[ seq.i ] <- .representative.seq;
            next;
        }
        # Are there any identical to this one?
        if( !any( dist.mat[ seq.i, -seq.i ] == 0 ) ) {
            next;
        }
        all.the.zeros <-
            colnames( dist.mat )[ -seq.i ][ dist.mat[ seq.i, -seq.i ] == 0 ];
        output.table.as.list[[ rownames( dist.mat )[ seq.i ] ]] <-
            c( colnames( dist.mat )[ seq.i ], paste( all.the.zeros, collapse = "," ) );
        #print( paste( all.the.zeros, collapse = ", " ) );
        seq.removed.since.represented.by.seq[ all.the.zeros ] <- seq.i;
    } # End foreach seq.i

    if( !is.null( output.table.file ) && ( length( output.table.as.list ) > 0 ) ) {
        output.table <- do.call( rbind, output.table.as.list );
        colnames( output.table ) <- c( "retained", "removed" );
        rownames( output.table ) <- NULL;
        write.table( output.table, file = output.table.file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE );
    }
    
    output.fasta <-
        input.fasta[ ( seq.removed.since.represented.by.seq == 0 ), , drop = FALSE ];
    if( add.copy.number.to.sequence.names ) {
        rownames( output.fasta ) <- sapply( rownames( output.fasta ), function ( .retained.seq.name ) {
            if( .retained.seq.name %in% names( output.table.as.list ) ) {
                #print( c( output.table.as.list[[ .retained.seq.name ]], gsub(  "[^,]", "", ( output.table.as.list[[ .retained.seq.name ]] )[ 2 ] ), nchar( gsub(  "[^,]", "", ( output.table.as.list[[ .retained.seq.name ]] )[ 2 ] ) ) + 2 ) );
                paste( .retained.seq.name, nchar( gsub(  "[^,]", "", ( output.table.as.list[[ .retained.seq.name ]] )[ 2 ] ) ) + 2, sep = "_" );
            } else {
                paste( .retained.seq.name, 1, sep = "_" );
            }
        } );
    }
               
    output.fasta.path <-
        paste( output.dir, "/", output.fasta.file, sep = "" );
    write.dna( output.fasta, output.fasta.path, format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = output.fasta.width );

    # Return the file name.
    return( output.fasta.path );
} # removeDuplicateSequencesFromAlignedFasta ( input.fasta.file, ... )

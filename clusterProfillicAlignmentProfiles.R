library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment", "consensus"
library( "dynamicTreeCut" ) # for "cutreeDynamic"

# This extracts the matrix of data from the alignment profile file; note that whenn we convert the strings to R floats we lose the infinite precision floats that profillic has (so we lose precision near 0)
getProfillicAlignmentProfile <- function ( alignment.profile.file ) {

  alignment.profile.file <- gsub( "^\\s*?(.+)\\s*?$", "\\1", alignment.profile.file, perl = TRUE ); # Remove extraneous whitespace
  alignment.profile.file.path <-
      gsub( "^(.*?)\\/[^\\/]+$", "\\1", alignment.profile.file );
  if( alignment.profile.file.path == alignment.profile.file ) {
      alignment.profile.file.path <- ".";
  }
  alignment.profile.file.short <-
      gsub( "^.*?\\/([^\\/]+)$", "\\1", alignment.profile.file );

  all.lines <- readLines( alignment.profile.file );

    parse.one.line <- function ( one.line ) {
        parts <- strsplit( substring( one.line, 3, nchar( one.line ) - 2 ), "\\),\\s*" )[[ 1 ]];        
        state.name <- gsub( "^(.).+$", "\\1", parts );
        parameter.type <- gsub( "^.(:|->).+$", "\\1", parts );

        parameter.values.parts <- strsplit( gsub( "^.(:|->)\\(([^\\)]+)\\)?$", "\\2", parts ), "," );
        names( parameter.values.parts ) <- sapply( 1:length( state.name ), function ( .i ) { paste( state.name[ .i ], parameter.type[ .i ] ); } );
        .rv <- 
        unlist( lapply( parameter.values.parts, function ( key.value.vector ) {
            .l <- strsplit( key.value.vector, "=" );
            .rv <- lapply( .l, function ( .pair ) { .pair[ 2 ] } );
            names( .rv ) <- lapply( .l, function ( .pair ) { .pair[ 1 ] } );
            return( .rv );
        } ) );
        mode( .rv ) <- "numeric";
        return( .rv );
    } # parse.one.line ( one.line )
    parse.all.lines <- Vectorize( parse.one.line, "one.line" );
    
    alignment.profile <- parse.all.lines( all.lines );
    alignment.profile <- t( alignment.profile );

    rownames( alignment.profile ) <- NULL;
    ## Replace those dots with spaces
    colnames( alignment.profile ) <- gsub( "^([^\\.]+)\\.(.)$", "\\1 \\2", colnames( alignment.profile ) );
    
    return( as.data.frame( alignment.profile ) );
} # getProfillicAlignmentProfile ( alignment.profile.file )

# OK; Now do a slew.  All calcaulated using the same profile, so all having the same dimensions.
# span.range.threshold controls the fraction of sites that will be used, based that fraction of the the range of the spans.
# default hclust.method is "average" for UPGMA. see help( "hclust" )
clusterProfillicAlignmentProfiles <- function ( alignment.profiles.filenames, input.fasta.file = NULL, output.dir = NULL, force.one.cluster = FALSE, span.range.threshold = 0.01, dist.method = "euclidean", hclust.method = "average" ) {
    
  if( !is.null( input.fasta.file ) ) {
    if( length( grep( "^(.*?)\\/[^\\/]+$", input.fasta.file ) ) == 0 ) {
        input.fasta.file.path <- ".";
    } else {
        input.fasta.file.path <-
            gsub( "^(.*?)\\/[^\\/]+$", "\\1", input.fasta.file );
    }
    input.fasta.file.short <-
        gsub( "^.*?\\/?([^\\/]+?)$", "\\1", input.fasta.file, perl = TRUE );
    
      if( is.null( output.dir ) ) {
          input.output.dir = input.fasta.file.path;
      } else {
          input.output.dir = output.dir;
      }
    input.fasta <- read.dna( input.fasta.file, format = "fasta" );
  } else {
      input.fasta <- NULL;
      input.output.dir <- output.dir;
  }
  ## Remove "/" from end of output.dir
  input.output.dir <-
    gsub( "^(.*?)\\/+$", "\\1", input.output.dir );
      
  if( force.one.cluster ) {
    # All in "cluster 0"
    clusters <- rep( 0, length( alignment.profiles.filenames ) );
    names( clusters ) <- names( alignment.profiles.filenames );
  } else {
    # Cluster, using the given distance and clustering methods, and with dynamic tree cutting.
    all.the.alignment.profiles <- lapply( alignment.profiles.filenames, getProfillicAlignmentProfile );
    names( all.the.alignment.profiles ) <- names( alignment.profiles.filenames );
    if( is.null( names( all.the.alignment.profiles ) ) ) {
        names( all.the.alignment.profiles ) <- alignment.profiles.filenames;
    }

    # Dumb, for now, try just unlisting them and treating all params the same, for eg euclidean distance
    all.the.vectors <- sapply( all.the.alignment.profiles, unlist );
    # Do not consider rows with no information, where "no" information means having a span less than .01% of the total range of non-zero spans.
    .span <- apply( all.the.vectors, 1, function( .row ) { diff( range( .row ) ) } )
    span.range.threshold.absolute <- ( diff( range( .span ) ) * span.range.threshold ) + base::min( .span );
    .excluded.count <- sum( .span <= span.range.threshold.absolute );
    .excluded.fraction <- mean( .span <= span.range.threshold.absolute );
    print( paste( "Excluding ", .excluded.count, " (", round( ( .excluded.fraction * 100 ), digits = 2 ), "%) of the parameters because their total spans are below the threshold of ", round( ( span.range.threshold * 100 ), digits = 2 ), "% of the range of spans -- they just don't vary enough to bother with." ) );

    subset.of.the.vectors <- all.the.vectors[ .span > span.range.threshold.absolute, , drop = FALSE ];
    
    subset.dist <- dist( t( subset.of.the.vectors ), method = dist.method );
    if( any( is.null( subset.dist ) ) || any( is.na( subset.dist ) ) || any( !is.finite( subset.dist ) ) ) {
      ## TODO: REMOVE
      warning( "UH OH got illegal distance value" );
      print( "UH OH got illegal distance value" );
      print( subset.dist );
    }
      
    ## TODO: REMOVE
    #print( dim( subset.dist ) );
    dendro <- hclust( subset.dist, method = hclust.method ); # "average" for UPGMA
    clusters <- suppressWarnings(
        cutreeDynamic(
            dendro, cutHeight = NULL, minClusterSize = 2,
            method = "hybrid", distM = as.matrix( subset.dist )
        ) );
    names( clusters ) <- dendro$labels;
  }
    
  # Write out one fasta file for each cluster, and a separate one containing just its consensus sequence.
  if( max( clusters ) - min( clusters ) != 0 ) {
      # Ensure that if there is more than one cluster, they are numbered from 1 (for some reason sometimes they are from 0)
      print( "Multiple clusters." );
      stopifnot( min( clusters ) < 2 ); # I'm assuming it's either 0 or 1:
      clusters <- clusters + min( clusters );
  } else {
      stopifnot( all( clusters == 0 ) );
  }
  for( cluster.i in min( clusters ):max( clusters ) ) {
      cluster.input.fasta <- input.fasta[ names( clusters[ clusters == cluster.i ] ), ];
  
      # Write the cluster alignment as a fasta file
      ## TODO: REMOVE
      #print( paste( input.output.dir, "/", input.fasta.file.short, ".cluster", cluster.i, ".fasta", sep = "" ) );
       write.dna( cluster.input.fasta, paste( input.output.dir, "/", input.fasta.file.short, ".cluster", cluster.i, ".fasta", sep = "" ), format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)
        
        # Write the cluster consensus as its own fasta file (note function default is to use majority consensus).
        .input.consensus <- as.DNAbin( matrix( seqinr::consensus( as.character( cluster.input.fasta ) ), nrow = 1 ) );
        rownames( .input.consensus ) <- paste( "Cluster", cluster.i, "consensus sequence" );
        write.dna( .input.consensus, paste( input.output.dir, "/", input.fasta.file.short, ".cluster", cluster.i, ".cons.fasta", sep = "" ), format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)
  } # End foreach cluster.i

  # Return the number of clusters (equivalently the index of the largest cluster)
  return( length( unique( clusters ) ) );
} # clusterProfillicAlignmentProfiles (..)

readFilenamesFromFile <- function ( input.file ) {
    return( readLines( input.file, warn = FALSE ) );
}

## Here is where the action is.
alignment.profile.files.list.file <- Sys.getenv( "clusterProfillicAlignmentProfiles_alignmentProfileFilesListFilename" ); # list of files containing outputs of `profileToAlignmentProfile --individual`
fasta.file <- Sys.getenv( "clusterProfillicAlignmentProfiles_fastaFilename" ); # output of `profileToAlignmentProfile --individual`
output.dir <- Sys.getenv( "clusterProfillicAlignmentProfiles_outputDir" ); # if null, gleaned from input filenames (and may differ for the insites cluster outputs and the original fasta file cluster outputs).
if( output.dir == "" ) {
    output.dir <- NULL;
}
force.one.cluster <- Sys.getenv( "clusterProfillicAlignmentProfiles_forceOneCluster" ); # don't actually cluster?
if( force.one.cluster == "" ) {
    force.one.cluster <- FALSE;
} else {
    force.one.cluster <- TRUE;
}

## TODO: REMOVE
warning( paste( "alignment profile files list input file:", alignment.profile.files.list.file ) );
if( file.exists( alignment.profile.files.list.file ) ) {
    alignment.profiles.files.list <- readFilenamesFromFile( alignment.profile.files.list.file );
    print( clusterProfillicAlignmentProfiles( alignment.profiles.files.list, fasta.file, output.dir, force.one.cluster ) );
} else {
    stop( paste( "File does not exist:", alignment.profiles.files.list.file ) );
}

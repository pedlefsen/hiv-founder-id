library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment", "consensus"
library( "dynamicTreeCut" ) # for "cutreeDynamic"

# Cluster the input sequences using Hamming distance and UPGMA, cutting the tree using "dynamic tree cutting" and use this to create cluster-specific subalignment fasta files and cluster subalignment consensus fasta files.
# Optionally provide an additional alignment with the same sequence names, for wich we will create analogous cluster-specific subalignment fasta files and cluster subalignment consensus fasta files.  Note that the names might have changed slightly - we will assume the sequences are in the same order and use this second file for the names to output if the file is present.
# if output.dir is null, output dirs will gleaned from input filenames (and may differ for the insites cluster outputs and the original fasta file cluster outputs).
# If force.one.cluster is TRUE, won't actually cluster - will just put all seqs in "cluster 0".
clusterSequences <- function ( insites.fasta.file, full.fasta.file = NULL, output.dir = NULL, force.one.cluster = FALSE ) {

    if( length( grep( "^(.*?)\\/[^\\/]+$", insites.fasta.file ) ) == 0 ) {
        insites.fasta.file.path <- ".";
    } else {
        insites.fasta.file.path <-
            gsub( "^(.*?)\\/[^\\/]+$", "\\1", insites.fasta.file );
    }
    insites.fasta.file.short <-
        gsub( "^.*?\\/?([^\\/]+?)$", "\\1", insites.fasta.file, perl = TRUE );

  if( is.null( output.dir ) ) {
      insites.output.dir <- insites.fasta.file.path;
  } else {
    ## Remove "/" from end of output.dir
      output.dir <-
          gsub( "^(.*?)\\/+$", "\\1", output.dir );
      insites.output.dir <- output.dir;
  }
  
  in.fasta <- read.dna( insites.fasta.file, format = "fasta" );
  if( !is.null( full.fasta.file ) ) {
    if( length( grep( "^(.*?)\\/[^\\/]+$", full.fasta.file ) ) == 0 ) {
        full.fasta.file.path <- ".";
    } else {
        full.fasta.file.path <-
            gsub( "^(.*?)\\/[^\\/]+$", "\\1", full.fasta.file );
    }
    full.fasta.file.short <-
        gsub( "^.*?\\/?([^\\/]+?)$", "\\1", full.fasta.file, perl = TRUE );
    
      if( is.null( output.dir ) ) {
          full.output.dir <- full.fasta.file.path;
      } else {
          full.output.dir <- output.dir;
      }
      full.fasta <- read.dna( full.fasta.file, format = "fasta" );
  } else {
      full.fasta <- NULL;
  }

  if( length( in.fasta ) == 0 ) {
    # No informative sites.  That's ok. But it means effectively we force one cluster.
    force.one.cluster <- 1;
    if( !is.null( full.fasta ) ) {
      in.fasta <- full.fasta;
    } else {
      in.fasta <- NULL;
    }
  } else if( !is.null( full.fasta ) ) {
    ## Fix labels.  Assuming input orders are the same.
    if( is.null( dim( full.fasta ) ) ) {
      # Note that apparently InSites only uses the first 1000, so the insites alignment can be shorter than the full one.
      rownames( in.fasta ) <- names( full.fasta )[ 1:nrow( in.fasta ) ];
    } else {
      rownames( in.fasta ) <- rownames( full.fasta )[ 1:nrow( in.fasta ) ];
    }
  }

  if( force.one.cluster ) {
    # All in "cluster 0"
    clusters <- rep( 0, nrow( in.fasta ) );
    names( clusters ) <- rownames( in.fasta );
  } else if( !is.null( in.fasta ) ) {
    # Cluster, using Hamming distance as the distance metric, UPGMA, and dynamic tree cutting.
    # The pairwise.deletion = TRUE argument is necessary so that columns with any gaps are not removed.
    # The second call adds in the count of sites with gap differences
    # (gap in one sequence but not the other), which are not included
    # in the first call's "raw" count. Adding them back in, we are
    # effectively treating those as differences.

    in.dist <- dist.dna( in.fasta, model = "raw", pairwise.deletion = TRUE );
    in.dist[ is.nan( in.dist ) ] <- 0;
    in.dist <- in.dist + dist.dna( in.fasta, model = "indel", pairwise.deletion = TRUE ); 
    if( any( is.null( in.dist ) ) || any( is.na( in.dist ) ) || any( !is.finite( in.dist ) ) ) {
      ## TODO: REMOVE
      warning( "UH OH got illegal distance value" );
      print( "UH OH got illegal distance value" );
      print( in.dist );
    }
      ## TODO: REMOVE
      #print( dim( in.dist ) );
    dendro <- hclust( in.dist, method = "average" ); # UPGMA
    clusters <- suppressWarnings(
        cutreeDynamic(
            dendro, cutHeight = NULL, minClusterSize = 2,
            method = "hybrid", distM = as.matrix( in.dist )
        ) );
    names( clusters ) <- dendro$labels;
  } else {
    die( "There are no sequences to work with; if there are no informative sites then you must supply the full alignment fasta file, and set force.one.cluster to TRUE." );
  }
  
  #clusters <- sort( clusters );
  
  # Write out one fasta file for each cluster, and a separate one containing just its consensus sequence.
  for( cluster.i in min( clusters ):max( clusters ) ) {
      cluster.fasta <- in.fasta[ names( clusters[ clusters == cluster.i ] ), ];

      # Write the cluster alignment as a fasta file
      cluster.fasta.file = paste( insites.output.dir, "/", insites.fasta.file.short, ".cluster", cluster.i, ".fasta", sep = "" );
      # warning( insites.output.dir );
      # warning( insites.fasta.file.short );
      # warning( cluster.fasta.file );
      write.dna( cluster.fasta, cluster.fasta.file, format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)
      
      # Write the cluster consensus as its own fasta file (note function default is to use majority consensus).
      .consensus <- as.DNAbin( matrix( seqinr::consensus( as.character( cluster.fasta ) ), nrow = 1 ) );
      rownames( .consensus ) <- paste( "Cluster", cluster.i, "consensus sequence" );
      write.dna( .consensus, paste( insites.output.dir, "/", insites.fasta.file.short, ".cluster", cluster.i, ".cons.fasta", sep = "" ), format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)

      if( !is.null( full.fasta ) ) {
        cluster.full.fasta <- full.fasta[ names( clusters[ clusters == cluster.i ] ), ];
  
        # Write the cluster alignment as a fasta file
        ## TODO: REMOVE
        #print( paste( full.output.dir, "/", full.fasta.file.short, ".cluster", cluster.i, ".fasta", sep = "" ) );
        write.dna( cluster.full.fasta, paste( full.output.dir, "/", full.fasta.file.short, ".cluster", cluster.i, ".fasta", sep = "" ), format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)
        
        # Write the cluster consensus as its own fasta file (note function default is to use majority consensus).
        .full.consensus <- as.DNAbin( matrix( seqinr::consensus( as.character( cluster.full.fasta ) ), nrow = 1 ) );
        rownames( .full.consensus ) <- paste( "Cluster", cluster.i, "consensus sequence" );
        write.dna( .full.consensus, paste( full.output.dir, "/", full.fasta.file.short, ".cluster", cluster.i, ".cons.fasta", sep = "" ), format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)
      }
  } # End foreach cluster.i

  # Return the number of clusters (equivalently the index of the largest cluster)
  return( length( unique( clusters ) ) );
} # clusterSequences ( insites.fasta.file, full.fasta.file, output.dir )

## Here is where the action is.
insites.fasta.file <- Sys.getenv( "clusterInformativeSites_inputFilename" ); # alignment of just informative sites
original.fasta.file <- Sys.getenv( "clusterInformativeSites_originalFastaFilename" ); # full alignment
if( original.fasta.file == "" ) {
    original.fasta.file <- NULL;
}
output.dir <- Sys.getenv( "clusterInformativeSites_outputDir" ); # if null, gleaned from input filenames (and may differ for the insites cluster outputs and the original fasta file cluster outputs).
if( output.dir == "" ) {
    output.dir <- NULL;
}
force.one.cluster <- Sys.getenv( "clusterInformativeSites_forceOneCluster" ); # don't actually cluster?
if( force.one.cluster == "" ) {
    force.one.cluster <- FALSE;
} else {
    force.one.cluster <- TRUE;
}
## TODO: REMOVE
# warning( paste( "insites alignment input file:", insites.fasta.file ) );
# warning( paste( "original alignment input file:", original.fasta.file ) );
# warning( paste( "output dir:", output.dir ) );
# warning( paste( "force one cluster?", force.one.cluster ) );

if( file.exists( insites.fasta.file ) ) {
    if( is.null( original.fasta.file ) || file.exists( original.fasta.file ) ) {
        print( clusterSequences( insites.fasta.file, original.fasta.file, output.dir, force.one.cluster ) );
    } else if( !is.null( original.fasta.file ) ) {
        stop( paste( "File does not exist:", original.fasta.file ) );
    }
} else {
    stop( paste( "File does not exist:", insites.fasta.file ) );
}

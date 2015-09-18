library("ape") # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr" ) # for "as.alignment", "consensus"
library( "dynamicTreeCut" ) # for "cutreeDynamic"

# Cluster the input sequences using Hamming distance and UPGMA, cutting the tree using "dynamic tree cutting" and use this to create cluster-specific subalignment fasta files and cluster subalignment consensus fasta files.
# Optionally provide an additional alignment with the same sequence names, for wich we will create analogous cluster-specific subalignment fasta files and cluster subalignment consensus fasta files.
# if output.dir is null, output dirs will gleaned from input filenames (and may differ for the insites cluster outputs and the original fasta file cluster outputs).
clusterSequences <- function ( insites.fasta.file, full.fasta.file = NULL, output.dir = NULL ) {

  insites.fasta.file.path <-
      gsub( "^(.*?)\\/[^\\/]+$", "\\1", insites.fasta.file );
  insites.fasta.file.short <-
      gsub( "^.*?\\/([^\\/]+)$", "\\1", insites.fasta.file );

  if( is.null( output.dir ) ) {
      insites.output.dir = insites.fasta.file.path;
  } else {
      insites.output.dir = output.dir;
  }
  
  in.fasta <- read.dna( insites.fasta.file, format = "fasta" );

  if( !is.null( full.fasta.file ) ) {
      full.fasta.file.path <-
          gsub( "/^(.*?)\\/[^\\/]+$/", "\\1", full.fasta.file );
      full.fasta.file.short <-
          gsub( "/^.*?\\/([^\\/]+)$/", "\\1", full.fasta.file );
    
      if( is.null( output.dir ) ) {
          full.output.dir = full.fasta.file.path;
      } else {
          full.output.dir = output.dir;
      }
      full.fasta <- read.dna( full.fasta.file, format = "fasta" );
  } else {
      full.fasta <- NULL;
  }
  
  # Using Hamming distance.
  in.dist <- dist.dna( in.fasta, model = "raw" );
  
  dendro <- hclust( in.dist, method = "average" ); # UPGMA
  
  clusters <-
      cutreeDynamic(
          dendro, cutHeight = NULL, minClusterSize = 2,
          method = "hybrid", distM = as.matrix( in.dist )
      );
  names( clusters ) <- dendro$labels;
  
  #clusters <- sort( clusters );
  
  # Write out one fasta file for each cluster, and a separate one containing just its consensus sequence.
  for( cluster.i in 1:max( clusters ) ) {
      cluster.fasta <- in.fasta[ names( clusters[ clusters == cluster.i ] ), ];

      # Write the cluster alignment as a fasta file
      cluster.fasta.file = paste( insites.output.dir, "/", insites.fasta.file.short, ".cluster", cluster.i, ".fasta", sep = "" );
      warning( insites.output.dir );
      warning( insites.fasta.file.short );
      warning( cluster.fasta.file );
      write.dna( cluster.fasta, cluster.fasta.file, format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)
      
      # Write the cluster consensus as its own fasta file (note function default is to use majority consensus).
      .consensus <- as.DNAbin( matrix( seqinr::consensus( as.character( cluster.fasta ) ), nrow = 1 ) );
      rownames( .consensus ) <- paste( "Cluster", cluster.i, "consensus sequence" );
      write.dna( .consensus, paste( insites.output.dir, "/", insites.fasta.file.short, ".cluster", cluster.i, ".cons.fasta", sep = "" ), format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)

      if( !is.null( full.fasta ) ) {
        cluster.full.fasta <- full.fasta[ names( clusters[ clusters == cluster.i ] ), ];
  
        # Write the cluster alignment as a fasta file
        write.dna( cluster.full.fasta, paste( full.output.dir, "/", full.fasta.file.short, ".cluster", cluster.i, ".fasta", sep = "" ), format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)
        
        # Write the cluster consensus as its own fasta file (note function default is to use majority consensus).
        .full.consensus <- as.DNAbin( matrix( seqinr::consensus( as.character( cluster.full.fasta ) ), nrow = 1 ) );
        rownames( .full.consensus ) <- paste( "Cluster", cluster.i, "consensus sequence" );
        write.dna( .full.consensus, paste( full.output.dir, "/", full.fasta.file.short, ".cluster", cluster.i, ".cons.fasta", sep = "" ), format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)
      }
  } # End foreach cluster.i

  # Return the number of clusters (equivalently the index of the largest cluster)
  return( max( clusters ) );
} # clusterSequences ( insites.fasta.file, full.fasta.file )

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
## TODO: REMOVE
warning( paste( "insites alignment input file:", insites.fasta.file ) );
warning( paste( "original alignment input file:", original.fasta.file ) );
warning( paste( "output dir:", output.dir ) );

if( file.exists( insites.fasta.file ) ) {
    if( is.null( original.fasta.file ) || file.exists( original.fasta.file ) ) {
        print( clusterSequences( insites.fasta.file, original.fasta.file, output.dir ) );
    } else if( !is.null( original.fasta.file ) ) {
        stop( paste( "File does not exist:", original.fasta.file ) );
    }
} else {
    stop( paste( "File does not exist:", insites.fasta.file ) );
}

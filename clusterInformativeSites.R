library("ape") # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr" ) # for "as.alignment", "consensus"
library( "dynamicTreeCut" ) # for "cutreeDynamic"

clusterSequences <- function ( input.fasta.file ) {

  in.fasta <- read.dna( input.fasta.file, format = "fasta" );
  
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
      write.dna( cluster.fasta, paste( input.fasta.file, ".cluster", cluster.i, ".fasta", sep = "" ), format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)
      
      # Write the cluster consensus as its own fasta file (note function default is to use majority consensus).
      .consensus <- as.DNAbin( matrix( seqinr::consensus( as.character( cluster.fasta ) ), nrow = 1 ) );
      rownames( .consensus ) <- paste( "Cluster", cluster.i, "consensus sequence" );
      write.dna( .consensus, paste( input.fasta.file, ".cluster", cluster.i, ".cons.fasta", sep = "" ), format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)
  } # End foreach cluster.i

  return( max( clusters ) );
} # clusterSequences ( input.fasta.file )

## Here is where the action is.
input.fasta.file <- Sys.getenv( "clusterInformativeSites_inputFilename" );
if( file.exists( input.fasta.file ) ) {
  print( clusterSequences( input.fasta.file ) );
} else {
  stop( paste( "File does not exist:", input.fasta.file ) );
}

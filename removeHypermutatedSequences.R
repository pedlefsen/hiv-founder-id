library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment", "consensus"

## This implements the default options of HYPERMUT 2.0 (http://www.hiv.lanl.gov/content/sequence/HYPERMUT/hypermut.html) and its use should cite: Rose, PP and Korber, BT. 2000. Detecting hypermutations in viral sequences with an emphasis on G -> A hypermutation. Bioinformatics 16(4): 400-401.
                                        # Note that for this (unlike the online version) we do not assume that the reference is in the file; instead we compute the consensus and use that.
                                        ## NOTE From Abrahams, 2009: "Sequences were analyzed for evidence of APOBEC3G-induced hypermutation by using the Hypermut 2.0 tool (www.hiv.lanl.gov). Sequences with a P value of 􏰐0.1 were considered enriched for mutations consistent with APOBEC3G sig- natures. In sequence sets showing evidence of enrichment for APOBEC3G- driven G-to-A transitions but with no single significantly hypermutated sequence, hypermutation was tested for after superimposition of all mutations within that sequence set onto a single representative sequence."
## Returns the number of hypermutated (and therefore removed) sequences.
removeHypermutatedSequences <- function ( fasta.file, output.dir = NULL, p.value.threshold = 0.1 ) {

    if( length( grep( "^(.*?)\\/[^\\/]+$", fasta.file ) ) == 0 ) {
        fasta.file.path <- ".";
    } else {
        fasta.file.path <-
            gsub( "^(.*?)\\/[^\\/]+$", "\\1", fasta.file );
    }
    fasta.file.short <-
        gsub( "^.*?\\/?([^\\/]+?)$", "\\1", fasta.file, perl = TRUE );
    fasta.file.short.nosuffix <-
        gsub( "^([^\\.]+)(\\..+)?$", "\\1", fasta.file.short, perl = TRUE );
    fasta.file.short.suffix <-
        gsub( "^([^\\.]+)(\\..+)?$", "\\2", fasta.file.short, perl = TRUE );

  if( is.null( output.dir ) ) {
      output.dir = fasta.file.path;
  }
  if( is.null( output.dir ) ) {
      output.dir = ".";
  }
  ## Remove "/" from end of output.dir
  output.dir <-
      gsub( "^(.*?)\\/+$", "\\1", output.dir );
  
  in.fasta <- read.dna( fasta.file, format = "fasta" );

  .consensus <- as.DNAbin( matrix( seqinr::consensus( as.character( in.fasta ) ), nrow = 1 ) );
  rownames( .consensus ) <- "Consensus sequence";

  compute.hypermut2.p.value <-
    function( seq.i ) {
      ## ERE I AM.  Bizarrely I can't quite reproduce what's on the web site.  I have more "potential" sites than are computed there, and I do not know why (todo: ask LANL)
      num.mut <- 0;
      num.potential.mut <- 0;
      num.control <- 0;
      num.potential.control <- 0;
      for( window.start.i in 1:( ncol( in.fasta ) - 2 ) ) {
          # if the window has any gaps in either sequence, skip it.
          if( any( as.character( in.fasta[ seq.i, window.start.i + 0:2 ] ) == "-" ) || any( as.character( .consensus[ 1, window.start.i + 0:2 ] ) == "-" ) ) {
              next;
          }
          #if( ( as.character( in.fasta[ seq.i, window.start.i + 0 ] ) == "a" ) && ( as.character( in.fasta[ seq.i, window.start.i + 1 ] ) %in% c( "a", "g" ) ) && ( as.character( in.fasta[ seq.i, window.start.i + 2 ] ) %in% c( "a", "g", "t" ) ) ) {
          ## Added that the ref has to be an a or a g.
          if( ( as.character( .consensus[ 1, window.start.i + 0 ] ) %in% c( "a", "g" ) ) && ( as.character( in.fasta[ seq.i, window.start.i + 0 ] ) == "a" ) && ( as.character( in.fasta[ seq.i, window.start.i + 1 ] ) %in% c( "a", "g" ) ) && ( as.character( in.fasta[ seq.i, window.start.i + 2 ] ) %in% c( "a", "g", "t" ) ) && ( in.fasta[ seq.i, window.start.i + 1 ] == .consensus[ 1, window.start.i + 1 ] ) && ( in.fasta[ seq.i, window.start.i + 2 ] == .consensus[ 1, window.start.i + 2 ] ) ) { # unsure whether we need to enforce no change in the "context" sites.
              num.potential.mut <- num.potential.mut + 1;
              #if( ( as.character( .consensus[ window.start.i + 0 ] ) == "g" ) && ( as.character( .consensus[ window.start.i + 1 ] ) %in% c( "a", "g" ) ) && ( as.character( .consensus[ window.start.i + 2 ] ) %in% c( "a", "g", "t" ) ) ) {
              if( ( as.character( .consensus[ 1, window.start.i + 0 ] ) == "g" ) ) { # don't enforce context in reference sequence
                  print( window.start.i );
                  print( as.character( in.fasta[ seq.i, window.start.i + 0:2 ] ) );
                  print( as.character( .consensus[ 1, window.start.i + 0:2 ] ) );
                  num.mut <- num.mut + 1;
              }
          }
                                        #if( ( as.character( in.fasta[ seq.i, window.start.i + 0 ] ) == "a" ) && ( ( as.character( in.fasta[ seq.i, window.start.i + 1 ] ) %in% c( "c", "t" ) ) || ( ( as.character( in.fasta[ seq.i, window.start.i + 1 ] ) %in% c( "a", "g" ) ) && ( as.character( in.fasta[ seq.i, window.start.i + 2 ] ) == "c" ) ) ) ) {
          ## Added that the ref has to be an a or a g.
          if( ( as.character( .consensus[ 1, window.start.i + 0 ] ) %in% c( "a", "g" ) ) && ( as.character( in.fasta[ seq.i, window.start.i + 0 ] ) == "a" ) && ( ( as.character( in.fasta[ seq.i, window.start.i + 1 ] ) %in% c( "c", "t" ) ) || ( ( as.character( in.fasta[ seq.i, window.start.i + 1 ] ) %in% c( "a", "g" ) ) && ( as.character( in.fasta[ seq.i, window.start.i + 2 ] ) == "c" ) ) ) && ( in.fasta[ seq.i, window.start.i + 1 ] == .consensus[ 1, window.start.i + 1 ] ) && ( in.fasta[ seq.i, window.start.i + 2 ] == .consensus[ 1, window.start.i + 2 ] ) ) { # unsure whether we need to enforce no change in the "context" sites.
              num.potential.control <- num.potential.control + 1;
              #if( ( as.character( .consensus[ window.start.i + 0 ] ) == "g" ) && ( as.character( .consensus[ window.start.i + 1 ] ) %in% c( "a", "g" ) ) && ( as.character( .consensus[ window.start.i + 2 ] ) %in% c( "a", "g", "t" ) ) ) {
              if( ( as.character( .consensus[ 1, window.start.i + 0 ] ) == "g" ) ) { # don't enforce context in reference sequence
                  print( window.start.i );
                  print( as.character( in.fasta[ seq.i, window.start.i + 0:2 ] ) );
                  print( as.character( .consensus[ 1, window.start.i + 0:2 ] ) );
                  num.control <- num.control + 1;
              }
          }
      }
      p.value <- fisher.test( ( matrix( c( num.control, ( num.potential.control - num.control ), num.mut, ( num.potential.mut - num.mut ) ), nrow = 2, byrow = T ) ) )$p.value;
      ## TODO: REMOVE
      # print( c( num.mut = num.mut, num.potential.mut = num.potential.mut ) );
      return( p.value );
    } # compute.hypermut2.p.value (..)

  exclude.sequence <- rep( FALSE, nrow( in.fasta ) );
  for( seq.i in 1:nrow( in.fasta ) ) {
    p.value <- compute.hypermut2.p.value( seq.i );
    if( p.value < p.value.threshold ) {
        print( paste( "Excluding", rownames( in.fasta )[ seq.i ], "because the pseudo-HYPERMUT2.0 p-value is", p.value, "." ) );
        exclude.sequence[ seq.i ] <- TRUE;
    }
  } # End foreach seq.i

  out.fasta <- in.fasta[ !exclude.sequence, ];

  # Write the subalignment as a fasta file
  out.fasta.file = paste( output.dir, "/", fasta.file.short.nosuffix, "_removeHypermutatedSequences", fasta.file.short.suffix, sep = "" );

  write.dna( out.fasta, out.fasta.file, format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)
    
  return( sum( exclude.sequence ) );
} # removeHypermutatedSequences (..)

## Here is where the action is.
fasta.file <- Sys.getenv( "removeHypermutatedSequences_inputFilename" ); # alignment of just informative sites
output.dir <- Sys.getenv( "removeHypermutatedSequences_outputDir" ); # if null, gleaned from input filenames (and may differ for the insites cluster outputs and the original fasta file cluster outputs).
if( output.dir == "" ) {
    output.dir <- NULL;
}
p.value.threshold <- Sys.getenv( "removeHypermutatedSequences_pValueThreshold" ); # how sensitive to be?
if( p.value.threshold == "" ) {
    p.value.threshold <- "0.1"; # Abrahams et al used liberal threshold of 0.1.
}

## TODO: REMOVE
# warning( paste( "alignment input file:", fasta.file ) );
# warning( paste( "output dir:", output.dir ) );
if( file.exists( fasta.file ) ) {
    print( removeHypermutatedSequences( fasta.file, output.dir, p.value.threshold = p.value.threshold ) );
} else {
    stop( paste( "File does not exist:", fasta.file ) );
}

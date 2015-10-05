library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment", "consensus"

# This parses the summaryTable output of the RAP tool; see runRAPOnline.pl.
removeRecombinedSequences <- function ( fasta.file, RAP.summaryTable.file, output.dir = NULL, p.value.threshold = 0.0007 ) {

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

    RAP.summaryTable <- readLines( RAP.summaryTable.file );
    if( length( RAP.summaryTable ) < 3 ) {
        stop( paste( "Unexpected: Got RAP.summaryTable:", RAP.summaryTable, sep = "\n" ) );
    }
    
  # The first two are header lines.
  exclude.sequence <- rep( FALSE, nrow( in.fasta ) );
  names( exclude.sequence ) <- rownames( in.fasta );
    for( line.i in 3:length( RAP.summaryTable ) ) {
        #warning( line.i );
        #warning( RAP.summaryTable[ line.i ] );
    p.value <- gsub( "^Set \\d+ \\S+ \\S+ ([0-9.]+) .+$", "\\1", RAP.summaryTable[ line.i ] );
    if( !is.na( p.value ) && ( p.value < p.value.threshold ) ) {
        seq.name <- gsub( "^Set \\d+ (\\S+) \\S+ [0-9.]+ .+$", "\\1", RAP.summaryTable[ line.i ] );
        seq.parents <- gsub( "^Set \\d+ \\S+ (\\S+) [0-9.]+ .+$", "\\1", RAP.summaryTable[ line.i ] );
        ## TODO: REMOVE
        warning( paste( "Excluding ", seq.name, " because the RAP p-value is ", p.value, ". It is a combination of ", seq.parents, ".", sep = "" ) );
        exclude.sequence[ seq.name ] <- TRUE;
    }
  } # End foreach line.i

  out.fasta <- in.fasta[ !exclude.sequence, ];

  # Write the subalignment as a fasta file
  out.fasta.file = paste( output.dir, "/", fasta.file.short.nosuffix, "_removeRecombinedSequences", fasta.file.short.suffix, sep = "" );

  write.dna( out.fasta, out.fasta.file, format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)

  return( sum( exclude.sequence ) );
} # removeRecombinedSequences (..)

## Here is where the action is.
fasta.file <- Sys.getenv( "removeRecombinedSequences_inputFilename" );
RAP.summaryTable.file <- Sys.getenv( "removeRecombinedSequences_RAPOutputFile" );
output.dir <- Sys.getenv( "removeRecombinedSequences_outputDir" );
if( output.dir == "" ) {
    output.dir <- NULL;
}
p.value.threshold <- Sys.getenv( "removeRecombinedSequences_pValueThreshold" ); # how sensitive to be?
if( p.value.threshold == "" ) {
    p.value.threshold <- "0.0007"; # Appears to be the suggestion from the output file "(summaryTable)"'s column header, which reads "Pvalues<0.0007".
}

## TODO: REMOVE
#warning( paste( "alignment input file:", fasta.file ) );
#warning( paste( "output dir:", output.dir ) );
if( file.exists( fasta.file ) ) {
    print( removeRecombinedSequences( fasta.file, RAP.summaryTable.file, output.dir, p.value.threshold = p.value.threshold ) );
} else {
    stop( paste( "File does not exist:", fasta.file ) );
}

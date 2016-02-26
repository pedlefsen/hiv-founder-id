library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment", "consensus"

# for removeDuplicateSequencesFromAlignedFasta(..)
source( "removeDuplicateSequencesFromAlignedFasta_safetosource.R" )

## 

#' Removes hyper mutated sequences
#'
#' This function implements the default options of HYPERMUT 2.0
#' (http://www.hiv.lanl.gov/content/sequence/HYPERMUT/hypermut.html) and its
#' use should cite: Rose, PP and Korber, BT. 2000. Detecting hypermutations in
#' viral sequences with an emphasis on G -> A hypermutation. Bioinformatics
#' 16(4): 400-401.
#'
#' Instead of removing hypermutated sequences, this function can also either
#' 'correct' them by replaceing the hypermutated As with Rs. (R IUPAC
#' ambiguiety character for A or G).
#'
#' Unlike the LANL implementation, the first sequence is not assumed to be the
#' reference. The consensus is computed and that is taken as the reference.
#'
#' Sequences are deduplicated before trying to detect hypermutation. This has
#' the side effect of creating some output files in the output directory
#' associated with the deduplication process.
#'
#' An alternative way to detect hypermutation is described in Abrahams, 2009:
#' "Sequences were analyzed for evidence of APOBEC3G-induced hypermutation by
#' using the Hypermut 2.0 tool (www.hiv.lanl.gov). Sequences with a P value of
#' 0.1 were considered enriched for mutations consistent with APOBEC3G sig-
#' natures. In sequence sets showing evidence of enrichment for APOBEC3G-
#' driven G-to-A transitions but with no single significantly hypermutated
#' sequence, hypermutation was tested for after superimposition of all
#' mutations within that sequence set onto a single representative sequence."
#'
#' @param fasta.file The input sequnce file name.
#' @param output.dir Target output directory.
#' @param p.value.threshold The cutoff for the Fisher Exact test to decide
#' whether or not hypermutation is present.
#' @param fix.instead.of.remove If TRUE, then hypermutated sequences will not
#' be removed. The hyper mutated bases will be replaced with the value in
#' 'fix.with'.
#' @param fix.with The letter to replace hypermutated bases with. Default r.
#' @param .consensus If left NULL, the consensus will be computed from
#' fasta.file, else the first sequence is assumed to the be reference sequence
#'
#' @return Returns the number of hypermutated (and therefore removed or fixed)
#' sequences. It also produces output files with the hypermutated sequences
#' removed.
#' @export

removeHypermutatedSequences <- function ( fasta.file, output.dir = NULL, p.value.threshold = 0.1, 
                                          fix.instead.of.remove = FALSE, fix.with = "r",
                                          .consensus = NULL ) {

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

  if( fix.instead.of.remove ){
      fix.with <- tolower(fix.with)
  }
  
  in.fasta <- read.dna( fasta.file, format = "fasta" );

  if (is.null(.consensus)){
    .consensus <- as.DNAbin( matrix( seqinr::consensus( as.character( in.fasta ) ), nrow = 1 ) )
    rownames( .consensus ) <- "Consensus sequence";
    in.fasta.tmp <- in.fasta
  } else {
    .consensus <- in.fasta[1,]
    in.fasta.tmp <- in.fasta[-1,]
  }

  # File copy hack to make removeDuplicateSequencesFromAlignedFasta use the
  # modified fasta file
  fasta.file.tmp <- paste(fasta.file.path, '/', fasta.file.short.nosuffix, '_tmp', fasta.file.short.suffix, sep = '')
  write.dna(in.fasta, fasta.file.tmp, format='fasta')
  write.dna(in.fasta.tmp, fasta.file, format='fasta')

    ## REMOVE DUPLICATES FIRST.
    ## If there are any duplicate sequences, remove them and
    ## incorporate the number of identical sequences into their names (_nnn suffix).
    # The output has the file name of the consensus file.
    fasta.file.no.duplicates <-
        removeDuplicateSequencesFromAlignedFasta( fasta.file, output.dir );
    fasta.file.no.duplicates.short <-
        gsub( "^.*?\\/?([^\\/]+?)$", "\\1", fasta.file.no.duplicates, perl = TRUE );
    fasta.file.no.duplicates.short.nosuffix <-
        gsub( "^(.*?)\\.[^\\.]+$", "\\1", fasta.file.no.duplicates.short, perl = TRUE );
    in.fasta.no.duplicates <-
        read.dna( fasta.file.no.duplicates, format = "fasta" );
    duplicate.sequences.tbl.file <-
        paste( output.dir, "/", fasta.file.no.duplicates.short.nosuffix, ".tbl", sep = "" );

  write.dna(in.fasta, fasta.file, format='fasta')
  file.remove(fasta.file.tmp)

  if( file.exists( duplicate.sequences.tbl.file ) ) { # won't exist if nothing was removed.
    duplicate.sequences.tbl.in <- read.table( file = duplicate.sequences.tbl.file, sep = "\t", header = TRUE );
    duplicate.sequences.tbl <- apply( duplicate.sequences.tbl.in, 1:2, function( .seq.name ) {
          # Fix our "escape" of the BAR and SLASH.
          .seq.name <- gsub( "-x-BAR-x-", "|", .seq.name );
          .seq.name <- gsub( "-x-SLASH-x-", "/", .seq.name );
          .seq.name <- gsub( "-x-BACKSLASH-x-", "\\", .seq.name );
          .seq.name <- gsub( "-x-DOT-x-", ".", .seq.name );
          .seq.name <- gsub( "-x-DOTDOT-x-", "..", .seq.name );
          return( .seq.name );
      } );
   } else {
       duplicate.sequences.tbl <- NULL;
   }
    
  compute.hypermut2.p.value <-
    function( seq.i, fix.sequence = FALSE ) {
      num.mut <- 0;
      num.potential.mut <- 0;
      num.control <- 0;
      num.potential.control <- 0;
      potential.pos <- data.frame(seq.name = character(0),
                                  pos = numeric(0),
                                  base.in.query = character(0),
                                  stringsAsFactors = F)
      for( window.start.i in 1:( ncol( in.fasta.no.duplicates ) - 2 ) ) {
          # if the window has any gaps in either sequence, skip it.
          # print(window.start.i)
          context.indx1 <- 1
          context.indx2 <- 2
          # Do not have to check position 0 for gaps since they will fail the == 'g' or == 'a' tests

          # Move the context forward if gaps are encountered to ensure that the
          # context pattern is matched to the sequence that the APOBEC enzyme
          # would have encountered.
          # Lots of IFs to prevent attempting to access values outside of the
          # valid range
          # I am sure that this can be done more elegantly...
          while( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + context.indx1 ] ) == "-" ){
              context.indx1 <- context.indx1 + 1
              context.indx2 <- context.indx2 + 1

              if (window.start.i + context.indx2 > ncol(in.fasta.no.duplicates)){
                break
              }
          }
          if (window.start.i + context.indx2 > ncol(in.fasta.no.duplicates)){
            next
          }
          while( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + context.indx2 ] ) == "-" ){
              context.indx2 <- context.indx2 + 1
              if (window.start.i + context.indx2 > ncol(in.fasta.no.duplicates)){
                break
              }
          }
          if (window.start.i + context.indx2 > ncol(in.fasta.no.duplicates)){
            next
          }

          if( ( as.character( .consensus[ 1, window.start.i + 0 ] ) == "g" ) && # Reference must mutate from G
              ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + context.indx1 ] ) %in% c( "a", "g" ) ) && # Context position 1 must match R = [AG] in query
              ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + context.indx2 ] ) %in% c( "a", "g", "t" ) ) ){ # Context position 2 must match D = [AGT] in query
              num.potential.mut <- num.potential.mut + 1;
              potential.pos <- rbind(potential.pos,
                data.frame(seq.name = row.names(in.fasta.no.duplicates)[seq.i],
                           pos = window.start.i,
                           base.in.query = as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 0 ] ),
                           stringsAsFactors = F))
              if( ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 0 ] ) == "a" ) ) { # If G -> A mutation occurred
                  #print( window.start.i );
                  #print( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 0:2 ] ) );
                  #print( as.character( .consensus[ 1, window.start.i + 0:2 ] ) );
                  num.mut <- num.mut + 1;
                  # NOTE: rbind is VERY slow - check here for performance issues
                  if( fix.sequence ) {
                      in.fasta.no.duplicates[ seq.i, window.start.i ] <<- as.DNAbin( fix.with );
                  }
              }
          }
          if( ( as.character( .consensus[ 1, window.start.i + 0 ] ) == "g" ) && # Reference must mutate from G
              ( ( ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + context.indx1 ] ) %in% c( "c", "t" ) ) && # Option 1 Context position 1 must match Y = [CT] in query
                  ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + context.indx2 ] ) %in% c( "a", "c", "g", "t" ) )) || # Option 1 Context position 2 must match N = [ACGT] in query
                ( ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + context.indx1 ] ) %in% c( "a", "g" ) ) && # Option 2 Context position 1 must match R = [AG] in query
                  ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + context.indx2 ] ) == "c" ) ) ) ){ # Option 2 Context position 2 must match C in query
              num.potential.control <- num.potential.control + 1;
              if( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 0 ] ) == "a" ) { # If G -> A mutation occureed
                  #print( window.start.i );
                  #print( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 0:2 ] ) );
                  #print( as.character( .consensus[ 1, window.start.i + 0:2 ] ) );
                  num.control <- num.control + 1;
              }
          }
      }
      p.value <- fisher.test( ( matrix( c( num.control, ( num.potential.control - num.control ), num.mut, ( num.potential.mut - num.mut ) ), nrow = 2, byrow = T ) ) )$p.value;
      ## TODO: REMOVE
      # print( row.names(in.fasta.no.duplicates)[seq.i])
      # print( c( num.mut = num.mut, num.potential.mut = num.potential.mut, num.control = num.control, num.potential.control = num.potential.control ) );
      row.names(potential.pos) <- NULL
      return( list(p.value = p.value,
                   potential.pos = potential.pos) );
    }; # compute.hypermut2.p.value (..)
    
  exclude.sequence <- rep( FALSE, nrow( in.fasta ) );
  names( exclude.sequence ) <- rownames( in.fasta );
  fixed.sequence <- rep( FALSE, nrow( in.fasta ) );
  names( fixed.sequence ) <- rownames( in.fasta );
  all.potential.pos.no.duplicates <- data.frame(seq.name = character(0),
                                                pos = numeric(0),
                                                base.in.query = character(0),
                                                stringsAsFactors = F)
  for( seq.i in 1:nrow( in.fasta.no.duplicates ) ) {
    .result.to.parse <- compute.hypermut2.p.value( seq.i )
    p.value <- .result.to.parse$p.value
    all.potential.pos.no.duplicates <- rbind(all.potential.pos.no.duplicates,
                                             .result.to.parse$potential.pos)
    rm(.result.to.parse)
    if( p.value < p.value.threshold ) {
        if( fix.instead.of.remove ) {
            .message <- paste( "Fixing", rownames( in.fasta.no.duplicates )[ seq.i ], "because the pseudo-HYPERMUT2.0 p-value is", p.value, "." );
            # Run it again but this time fix it.
            .result.ignored <-
                compute.hypermut2.p.value( seq.i, fix.sequence = TRUE );
            fixed.sequence[ rownames( in.fasta.no.duplicates )[ seq.i ] ] <- TRUE;
        } else {
            .message <- paste( "Excluding", rownames( in.fasta.no.duplicates )[ seq.i ], "because the pseudo-HYPERMUT2.0 p-value is", p.value, "." );
            exclude.sequence[ rownames( in.fasta.no.duplicates )[ seq.i ] ] <- TRUE;
        }
        ## TODO: REMOVE
        cat( .message, fill = TRUE );
    }
  } # End foreach seq.i
  all.potential.pos <- all.potential.pos.no.duplicates

  if( is.null( duplicate.sequences.tbl ) ) {
    out.fasta <- in.fasta;
  } else {
    if( fix.instead.of.remove ) {
      out.fasta <- in.fasta;
    }
    #print( names( exclude.sequence )[ exclude.sequence ] );
    #print( duplicate.sequences.tbl[ , "retained" ] );
    excluded.sequences.with.duplicates <- intersect( names( exclude.sequence )[ exclude.sequence ], duplicate.sequences.tbl[ , "retained" ] );
    newly.excluded.sequence <- rep( FALSE, nrow( in.fasta ) );
    names( newly.excluded.sequence ) <- rownames( in.fasta );

    # detroying Paul's beautiful apply-based design with some nested for loops
    for (.retained.sequence in duplicate.sequences.tbl[,1]){
      for (.duplicate.sequence in strsplit(duplicate.sequences.tbl[,2][duplicate.sequences.tbl[,1] == .retained.sequence], ',')){
        all.potential.pos <- rbind(all.potential.pos,
                                   data.frame(seq.name = .duplicate.sequence,
                                              pos = all.potential.pos.no.duplicates$pos[all.potential.pos.no.duplicates$seq.name == .retained.sequence],
                                              base.in.query = all.potential.pos.no.duplicates$base.in.query[all.potential.pos.no.duplicates$seq.name == .retained.sequence],
                                              stringsAsFactors = F))
      }
    }

    while( length( excluded.sequences.with.duplicates ) > 0 ) {
        newly.excluded.sequence[ newly.excluded.sequence ] <- FALSE;
        rownames( duplicate.sequences.tbl ) <- duplicate.sequences.tbl[ , "retained" ];
        .duplicates.strings <- duplicate.sequences.tbl[ excluded.sequences.with.duplicates, "removed" ];
        .duplicates.list <- strsplit( .duplicates.strings, "," );
        names( .duplicates.list ) <- excluded.sequences.with.duplicates;
        .result.ignored <- 
        lapply( excluded.sequences.with.duplicates, function( .retained.sequence ) {
            if( length( .duplicates.list[[ .retained.sequence ]] ) > 0 ) {
                if( fix.instead.of.remove ) {
                  ## TODO: REMOVE
                  cat( paste( "Fixing '", .duplicates.list[[ .retained.sequence ]], "' because it is a duplicate of ", .retained.sequence, ".", sep = "", collapse = "\n" ), fill = TRUE );
                  # Copy the fixed sequence.
                  .result.ignored <- 
                  sapply( .duplicates.list[[ .retained.sequence ]], function ( .retained ) {
                      out.fasta[ .retained, ] <<- out.fasta[ .retained.sequence, ];
                      fixed.sequence[ .duplicates.list[[ .retained.sequence ]] ] <<- TRUE;
                      return( NULL );
                  } );
              } else {
                  ## TODO: REMOVE
                  cat( paste( "Excluding '", .duplicates.list[[ .retained.sequence ]], "' because it is a duplicate of ", .retained.sequence, ".", sep = "", collapse = "\n" ), fill = TRUE );
              }
            }
            exclude.sequence[ .duplicates.list[[ .retained.sequence ]] ] <<- TRUE;
            newly.excluded.sequence[ .duplicates.list[[ .retained.sequence ]] ] <<- TRUE;
            return( NULL );
        } );
        # Recurse in case some of the sequences just excluded were recursively removed.
        excluded.sequences.with.duplicates <- intersect( names( newly.excluded.sequence )[ newly.excluded.sequence ], duplicate.sequences.tbl[ , "retained" ] );
    } # End while( length( excluded.sequences.with.duplicates ) > 0 )
    if( !fix.instead.of.remove ) {
        out.fasta <- in.fasta[ !exclude.sequence, , drop = FALSE ];
    }
  } # End if is.null( duplicate.sequences.tbl ) .. else ..

  # Write the subalignment as a fasta file
  if( fix.instead.of.remove ) {
      out.fasta.file = paste( output.dir, "/", fasta.file.short.nosuffix, "_fixHypermutatedSequencesWith", fix.with, fasta.file.short.suffix, sep = "" );
  } else {
      out.fasta.file = paste( output.dir, "/", fasta.file.short.nosuffix, "_removeHypermutatedSequences", fasta.file.short.suffix, sep = "" );
  }

  write.dna( out.fasta, out.fasta.file, format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)
  write.csv(all.potential.pos, paste( output.dir, "/", fasta.file.short.nosuffix, "_potentialHypermutatedPositions.csv", sep = ''), row.names=F)
    
  if( fix.instead.of.remove ) {
      .rv <- sum( fixed.sequence );
  } else {
      .rv <- sum( exclude.sequence );
  }
    
  return( .rv );
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
fix.instead.of.remove <- Sys.getenv( "removeHypermutatedSequences_fixInsteadOfRemove" );
if( ( nchar( fix.instead.of.remove ) == 0 ) || ( fix.instead.of.remove == "0" ) || ( toupper( fix.instead.of.remove ) == "F" ) || ( toupper( fix.instead.of.remove ) == "FALSE" ) ) {
    fix.instead.of.remove <- FALSE;
} else {
    fix.instead.of.remove <- TRUE;
}
fix.with <- Sys.getenv( "removeHypermutatedSequences_fixWith" ); # with what should we fix? "R" seems right but "-" might be best for algorithms...
if( fix.with == "" ) {
    fix.with <- "R";
}

## TODO: REMOVE
# warning( paste( "alignment input file:", fasta.file ) );
# warning( paste( "output dir:", output.dir ) );
if( file.exists( fasta.file ) ) {
    print( removeHypermutatedSequences( fasta.file, output.dir, p.value.threshold = p.value.threshold, fix.instead.of.remove = fix.instead.of.remove, fix.with = fix.with ) );
} else {
    stop( paste( "File does not exist:", fasta.file ) );
}

## TODO: remove
## Testing scribbles
#fasta.file = '/tmp/hypermut_lanl_produces_expected_results_dups.fasta'
#removeHypermutatedSequences('/tmp/hypermut_lanl_produces_expected_results.fasta', .consensus = 'first')
#removeHypermutatedSequences('/tmp/hypermut_lanl_produces_expected_results_dups.fasta', .consensus = 'first')

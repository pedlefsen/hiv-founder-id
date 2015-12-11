library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment"
library( "Biostrings" ) # for "pairwiseAlignment"

# This compares two nucleotide fasta files, each containing one or more sequences (aligned or not, it doesn't matter; gaps will be stripped internally).  The comparison is conducted in both nucleotide and amino acid space after gene-cutting, codon-aligning, and translating the sequences using GeneCutter at LANL (see runGeneCutterOnline.pl).
# If the output file exists and append is true, the output will be appended without printing out a new header row.
# set recreate.ungapped.fastas = TRUE to force recreation of ungapped versions of the input file (in the same dir as the input file) in the event that it's already there (by default it'll use the existing file).
evaluateFounders <- function ( estimates.fasta.file, truths.fasta.file, output.dir = NULL, output.file = NULL, output.file.append = FALSE, output.fasta.width = 72, recreate.ungapped.fastas = FALSE, genecutter.proteins.list = "-GAG-POL-VIF-VPR-TAT-REV-VPU-ENV-NEF", genecutter.genome.region = "ALL" ) {

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
            output.dir <- ".";
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

    combined.ungapped.fasta.file.nosuffix <- paste( output.dir, "/", estimates.fasta.file.short.nosuffix, "_ungapped_with_", truths.fasta.file.short.nosuffix, "_ungapped_combined", sep = "" );
    nucleotides.dir <- paste( combined.ungapped.fasta.file.nosuffix, "_allnucs", sep = "" );
    proteins.dir <- paste( combined.ungapped.fasta.file.nosuffix, "_allproteins", sep = "" );
    
    estimates.fasta <- NULL;
    try.count <- 0;
    TOO.MANY.TRIES <- 100;
    while( ( try.count < TOO.MANY.TRIES ) && ( is.null( estimates.fasta ) || is.null( nrow( estimates.fasta ) ) || ( nrow( estimates.fasta ) == 0 ) ) ) {
       if( try.count > 0 ) {
           warning( paste( "Trouble reading results in file", estimates.fasta.file, ".. trying again in a second." ) );
           Sys.sleep( 1 );
       }
       estimates.fasta <- read.dna( estimates.fasta.file, format = "fasta" );
    }
    if( try.count >= TOO.MANY.TRIES ) {
        stop( paste( "COULD NOT READ results in file", estimates.fasta.file, ".. GIVING UP." ) );
    }
    ## Make an ungapped version, if it doesn't already exist. (in output.dir)
    estimates.fasta.file.ungapped <-
      paste( output.dir, "/", estimates.fasta.file.short.nosuffix, "_ungapped", estimates.fasta.file.suffix, sep = "" );
    if( recreate.ungapped.fastas || !file.exists( estimates.fasta.file.ungapped ) ) {
          # Create an ungapped version, and save it.
          .estimates.fasta.as.character <- as.character( estimates.fasta );
          .result.ignored <- lapply( 1:nrow( .estimates.fasta.as.character ), function( .row.i ) {
              .row.seq.chars <- .estimates.fasta.as.character[ .row.i, ];
              .row.seq.ungapped <- as.DNAbin( matrix( .row.seq.chars[ .row.seq.chars != "-" ], nrow = 1 ) );
              # While we're at it, remove spaces from the header because GeneCutter might get confused if the seqnames (before the first space) are non-unique.
              rownames( .row.seq.ungapped ) <- gsub( "\\s", "_", rownames( estimates.fasta )[ .row.i ] );
              write.dna( .row.seq.ungapped, file = estimates.fasta.file.ungapped, format = "fasta", append = ( .row.i > 1 ), colsep = "", indent = "", blocksep = 0, nbcol = 1, colw = output.fasta.width );
              return( NULL );
          } );
    }

    ## This is a hex to ward off the bad juju in the filesystem (#voodoo)
    truths.fasta <- NULL;
    try.count <- 0;
    TOO.MANY.TRIES <- 100;
    while( ( try.count < TOO.MANY.TRIES ) && ( is.null( truths.fasta ) || is.null( nrow( truths.fasta ) ) || ( nrow( truths.fasta ) == 0 ) ) ) {
       if( try.count > 0 ) {
           warning( paste( "Trouble reading results in file", truths.fasta.file, ".. trying again in a second." ) );
           Sys.sleep( 1 );
       }
       truths.fasta <- read.dna( truths.fasta.file, format = "fasta" );
    }
    if( nrow( truths.fasta ) == 0 ) {
        stop( paste( "There are NO TRUE FOUNDERS in file", truths.fasta.file ) );
    }
  
    ## Make an ungapped version, if it doesn't already exist.
    truths.fasta.file.ungapped <- paste( output.dir, "/", truths.fasta.file.short.nosuffix, "_ungapped", truths.fasta.file.suffix, sep = "" );
    if( recreate.ungapped.fastas || !file.exists( truths.fasta.file.ungapped ) ) {
        # Create an ungapped version, and save it.
        .truths.fasta.as.character <- as.character( truths.fasta );
        .result.ignored <- lapply( 1:nrow( .truths.fasta.as.character ), function( .row.i ) {
            .row.seq.chars <- .truths.fasta.as.character[ .row.i, ];
            .row.seq.ungapped <- as.DNAbin( matrix( .row.seq.chars[ .row.seq.chars != "-" ], nrow = 1 ) );
            # While we're at it, remove spaces from the header because GeneCutter might get confused if the seqnames (before the first space) are non-unique.
            rownames( .row.seq.ungapped ) <- gsub( "\\s", "_", rownames( truths.fasta )[ .row.i ] );
            write.dna( .row.seq.ungapped, file = truths.fasta.file.ungapped, format = "fasta", append = ( .row.i > 1 ), colsep = "", indent = "", blocksep = 0, nbcol = 1, colw = output.fasta.width );
            return( NULL );
        } );
    }
      
    ## Put together the two fasta files, the lazy way.
    combined.ungapped.fasta.file <-
        paste( combined.ungapped.fasta.file.nosuffix, truths.fasta.file.suffix, sep = "" );
    system( paste( "cp", truths.fasta.file.ungapped, combined.ungapped.fasta.file ) );
    system( paste( "cat", estimates.fasta.file.ungapped, ">>", combined.ungapped.fasta.file ) );
    if( !file.exists( nucleotides.dir ) || !file.exists( proteins.dir ) ) {
      nucleotides.zipfile <- paste( combined.ungapped.fasta.file.nosuffix, "_allnucs.zip", sep = "" );
      proteins.zipfile <- paste( combined.ungapped.fasta.file.nosuffix, "_allproteins.zip", sep = "" );
      if( !file.exists( nucleotides.zipfile ) || !file.exists( proteins.zipfile ) ) {
        ## Run it through GeneCutter.
        ## TODO: Add other comparison seqs first, eg CON_M
        system( paste( "perl runGeneCutterOnline.pl -P", genecutter.proteins.list, "-R", genecutter.genome.region, "-V", combined.ungapped.fasta.file, output.dir ) );
        stopifnot( file.exists( nucleotides.zipfile ) );
        stopifnot( file.exists( proteins.zipfile ) );
      }

      ## Now unzip them.
      system( paste( "rm -rf", nucleotides.dir ) );
      system( paste( "rm -rf", proteins.dir ) );
      ## NOTE: IF the protein list has only one element, this isn't a zipfile, it's just a fasta file.
      if( sum( sapply( strsplit( genecutter.proteins.list, "-" )[[1]], function( .x ) { nchar( .x ) > 0 } ) ) == 1 ) {
          dir.create( nucleotides.dir );
          cat( unlist( system( paste( "cp", nucleotides.zipfile, paste( nucleotides.dir, "/", genecutter.proteins.list, ".NA.FASTA", sep = "" ), sep = " " ) ) ), fill = F );
          dir.create( proteins.dir );
          cat( unlist( system( paste( "cp", proteins.zipfile, paste( proteins.dir, "/", genecutter.proteins.list, ".AA.FASTA", sep = "" ), sep = " " ) ) ), fill = F );
      } else {
          system( paste( "unzip", nucleotides.zipfile, "-d", nucleotides.dir, sep = " " ) )
          system( paste( "unzip", proteins.zipfile, "-d", proteins.dir, sep = " " ) )
      }
    } # End if the dirs don't already exist ...
    # We use the fact that the first num.reference.sequences sequences are the HXB2 and then the "truths" (rather than use sequence names, which can get mangled).
    num.references <- 1 + nrow( truths.fasta );
    num.estimates <- nrow( estimates.fasta );
    evaluateOneAlignedPair <- function ( pair.of.chars, include.gaps = FALSE ) {
        stopifnot( length( pair.of.chars ) == 2 );
        if( all( pair.of.chars == "-" ) ) {
            return( NA );
        }                            
        if( any( pair.of.chars == "-" ) ) {
            if( include.gaps ) {
                return( 1 );
            } else {
                return( NA );
            }
        }
        if( pair.of.chars[ 1 ] == pair.of.chars[ 2 ] ) {
            return( 0 );
        } else {
            return( 1 );
        }
    } # evaluateOneAlignedPair (..)
    
    computeRelevantDistancesForOneFile <- function ( filename, file.is.aa ) {
        if( !file.exists( filename ) ) {
            stop( "ERROR: FILE MISSING" );
            return( list( HD.includingGaps = NA, denominator.includingGaps = NA, HD.ignoringGaps = NA, denominator.ignoringGaps = NA ) );
        }
        
        fasta.in <- NULL;
        tryCatch( {
        if( file.is.aa ) {
            fasta.in <-
                readAAMultipleAlignment( filename, "fasta" );
        } else {
            fasta.in <-
                readDNAMultipleAlignment( filename, "fasta" );
        } }, error = function( e ) {
            #cat( unlist( e ), fill = TRUE );
            return( list( HD.includingGaps = NA, denominator.includingGaps = NA, HD.ignoringGaps = NA, denominator.ignoringGaps = NA ) );
        } );

        if( is.null( fasta.in ) || ( nrow( fasta.in ) == 0 ) ) {
          # no sequences in the file?!
          return( list( HD.includingGaps = NA, denominator.includingGaps = NA, HD.ignoringGaps = NA, denominator.ignoringGaps = NA ) );
        }
        fasta.mat <- t( apply( as.matrix( 1:nrow( fasta.in ) ), 1, function ( i ) {
            s2c( as.character( unmasked( fasta.in )[ i, ] ) )
        } ) );
        ## Sometimes for some reason genecutter puts some extra (translated) seqs in there.
        if( length( unique( rownames( fasta.in ) ) ) != length( ( rownames( fasta.in ) ) ) ) {
            # There are non-unique entries.  Use the first one, for now.
            ## TODO: ?
            warning( paste( "Duplicate sequence names in fasta file", fasta.in ) );
            ## TODO: Check that the first one is consistently the best; seems so.
        }
        #stopifnot( nrow( fasta.mat ) == ( num.references + num.estimates ) );
        # for each combination of reference and target, compute the Hamming distance, with and without gaps.
        hamming.distances.ignoring.gaps <-
            matrix( nrow = num.estimates, ncol = num.references );
        rownames( hamming.distances.ignoring.gaps ) <-
            rownames( estimates.fasta );
        colnames( hamming.distances.ignoring.gaps ) <-
            c( "HXB2", rownames( truths.fasta ) );
        denominator.ignoring.gaps <-
            matrix( nrow = num.estimates, ncol = num.references );
        rownames( denominator.ignoring.gaps ) <-
            rownames( estimates.fasta );
        colnames( denominator.ignoring.gaps ) <-
            c( "HXB2", rownames( truths.fasta ) );
        hamming.distances.including.gaps <-
            matrix( nrow = num.estimates, ncol = num.references );
        rownames( hamming.distances.including.gaps ) <-
            rownames( estimates.fasta );
        colnames( hamming.distances.including.gaps ) <-
            c( "HXB2", rownames( truths.fasta ) );
        denominator.including.gaps <-
            matrix( nrow = num.estimates, ncol = num.references );
        rownames( denominator.including.gaps ) <-
            rownames( estimates.fasta );
        colnames( denominator.including.gaps ) <-
            c( "HXB2", rownames( truths.fasta ) );
        for( est.i in 1:num.estimates ) {
            for( ref.j in 1:num.references ) {
                .ignoring.gaps <- apply( fasta.mat[ c( ref.j, num.references + est.i ), ], 2, evaluateOneAlignedPair, include.gaps = FALSE );
                .including.gaps <- apply( fasta.mat[ c( ref.j, num.references + est.i ), ], 2, evaluateOneAlignedPair, include.gaps = TRUE );
                hamming.distances.ignoring.gaps[ est.i, ref.j ] <-
                    sum( .ignoring.gaps, na.rm = T );
                denominator.ignoring.gaps[ est.i, ref.j ] <-
                    sum( !is.na( .ignoring.gaps ) );
                hamming.distances.including.gaps[ est.i, ref.j ] <-
                    sum( .including.gaps, na.rm = T );
                denominator.including.gaps[ est.i, ref.j ] <-
                    sum( !is.na( .including.gaps ) );
            } # End foreach ref.j
        } # End foreach est.i
        return( list( HD.includingGaps = hamming.distances.including.gaps, denominator.includingGaps = denominator.including.gaps, HD.ignoringGaps = hamming.distances.ignoring.gaps, denominator.ignoringGaps = denominator.ignoring.gaps ) );
    } # computeRelevantDistancesForOneFile (..)

    ## When there is no non-gap sequence data for a gene or another error occurs, the specific fasta file will not contain fasta but instead will contain the error message.
    ## This returns TRUE if the file is a real fasta file (that is, if its first character is ">")
    checkFastaFileIsReal <- function ( filename ) {
        first.line <- readLines( filename, n = 1 );
        if( length( grep( "^>", first.line ) ) > 0 ) {
            return( TRUE );
        } else {
            return( FALSE );
        }
    } # checkFastaFileIsReal (..)
    
    evaluateOneFile <- function ( filename, file.is.aa ) {
        relevant.distances <-
            computeRelevantDistancesForOneFile( filename, file.is.aa );

        ## includingGaps:
         average.HD.includingGaps <-
             relevant.distances$HD.includingGaps[ , -1, drop = FALSE ] /
                 relevant.distances$denominator.includingGaps[ , -1, drop = FALSE ];
        ## ignoringGaps:
         average.HD.ignoringGaps <-
             relevant.distances$HD.ignoringGaps[ , -1, drop = FALSE ] /
                 relevant.distances$denominator.ignoringGaps[ , -1, drop = FALSE ];
  
        ## includingGaps:
        ## These are the two "perspectives": average of nearest estimate over truths, or average of nearest truth over estimates.
        nearest.founder.index.by.reference <- apply( average.HD.includingGaps, 2, base::which.min );
        nearest.founder.index.by.estimate <- apply( average.HD.includingGaps, 1, base::which.min );
        # For this, skip col 1, HXB2.
        truths.nearest.founder.HD.sum.includingGaps <-
            sum( sapply( 2:ncol( relevant.distances$HD.includingGaps ), function( .i ) { relevant.distances$HD.includingGaps[ nearest.founder.index.by.reference[ .i - 1 ], .i ] } ) );
        truths.nearest.founder.denominator.sum.includingGaps <-
            sum( sapply( 2:ncol( relevant.distances$denominator.includingGaps ), function( .i ) { relevant.distances$denominator.includingGaps[ nearest.founder.index.by.reference[ .i - 1 ], .i ] } ) );
        truths.nearest.founder.average.HD.includingGaps <-
            truths.nearest.founder.HD.sum.includingGaps / truths.nearest.founder.denominator.sum.includingGaps;
        estimates.nearest.founder.HD.sum.includingGaps <-
            sum( sapply( 1:nrow( relevant.distances$HD.includingGaps ), function( .i ) { relevant.distances$HD.includingGaps[ .i, nearest.founder.index.by.estimate[ .i ] ] } ) );
        estimates.nearest.founder.denominator.sum.includingGaps <-
            sum( sapply( 1:nrow( relevant.distances$denominator.includingGaps ), function( .i ) { relevant.distances$denominator.includingGaps[ .i, nearest.founder.index.by.estimate[ .i ] ] } ) );
        estimates.nearest.founder.average.HD.includingGaps <-
            estimates.nearest.founder.HD.sum.includingGaps / estimates.nearest.founder.denominator.sum.includingGaps;
        
        ## TODO: REMOVE
        # print( truths.nearest.founder.average.HD.includingGaps );
        # print( estimates.nearest.founder.average.HD.includingGaps );
        # print( nearest.founder.average.HD.includingGaps );
        
        ## ignoringGaps:
        ## These are the two "perspectives": average of nearest estimate over truths, or average of nearest truth over estimates.
        nearest.founder.index.by.reference <- apply( average.HD.ignoringGaps, 2, base::which.min );
        nearest.founder.index.by.estimate <- apply( average.HD.ignoringGaps, 1, base::which.min );
        # For this, skip col 1, HXB2.
        truths.nearest.founder.HD.sum.ignoringGaps <-
            sum( sapply( 2:ncol( relevant.distances$HD.ignoringGaps ), function( .i ) { relevant.distances$HD.ignoringGaps[ nearest.founder.index.by.reference[ .i - 1 ], .i ] } ) );
        truths.nearest.founder.denominator.sum.ignoringGaps <-
            sum( sapply( 2:ncol( relevant.distances$denominator.ignoringGaps ), function( .i ) { relevant.distances$denominator.ignoringGaps[ nearest.founder.index.by.reference[ .i - 1 ], .i ] } ) );
        truths.nearest.founder.average.HD.ignoringGaps <-
            truths.nearest.founder.HD.sum.ignoringGaps / truths.nearest.founder.denominator.sum.ignoringGaps;
        estimates.nearest.founder.HD.sum.ignoringGaps <-
            sum( sapply( 1:nrow( relevant.distances$HD.ignoringGaps ), function( .i ) { relevant.distances$HD.ignoringGaps[ .i, nearest.founder.index.by.estimate[ .i ] ] } ) );
        estimates.nearest.founder.denominator.sum.ignoringGaps <-
            sum( sapply( 1:nrow( relevant.distances$denominator.ignoringGaps ), function( .i ) { relevant.distances$denominator.ignoringGaps[ .i, nearest.founder.index.by.estimate[ .i ] ] } ) );
        estimates.nearest.founder.average.HD.ignoringGaps <-
            estimates.nearest.founder.HD.sum.ignoringGaps / estimates.nearest.founder.denominator.sum.ignoringGaps;

        ## ERE I AM TRYING TO TRACK DOWN APPARENTLY GREATER assymetry than I would have anticipated; shouldn't the perspectives give the same results (as each other) when there is a single-vs-single founder comparison?
        
        ## TODO: REMOVE
        # print( truths.nearest.founder.average.HD.ignoringGaps );
        # print( estimates.nearest.founder.average.HD.ignoringGaps );
        # print( nearest.founder.average.HD.ignoringGaps );

        .list <- c( truths.perspective.HD.average.ignoringGaps = truths.nearest.founder.average.HD.ignoringGaps, truths.perspective.HD.sum.ignoringGaps = truths.nearest.founder.HD.sum.ignoringGaps, truths.perspective.denominator.sum.ignoringGaps = truths.nearest.founder.denominator.sum.ignoringGaps, estimates.perspective.HD.average.ignoringGaps = estimates.nearest.founder.average.HD.ignoringGaps, estimates.perspective.HD.sum.ignoringGaps = estimates.nearest.founder.HD.sum.ignoringGaps, estimates.perspective.denominator.sum.ignoringGaps = estimates.nearest.founder.denominator.sum.ignoringGaps,
                   truths.perspective.HD.average.includingGaps = truths.nearest.founder.average.HD.includingGaps, truths.perspective.HD.sum.includingGaps = truths.nearest.founder.HD.sum.includingGaps, truths.perspective.denominator.sum.includingGaps = truths.nearest.founder.denominator.sum.includingGaps, estimates.perspective.HD.average.includingGaps = estimates.nearest.founder.average.HD.includingGaps, estimates.perspective.HD.sum.includingGaps = estimates.nearest.founder.HD.sum.includingGaps, estimates.perspective.denominator.sum.includingGaps = estimates.nearest.founder.denominator.sum.includingGaps );
        return( sapply( .list, function( .avg ) { sprintf( "%1.8f", .avg ) } ) );
    } # evaluateOneFile (..)

    evaluate.results.by.protein.file <- 
        lapply( dir( proteins.dir, full.names = T ), function ( .file ) {
            print( .file );
            if( checkFastaFileIsReal( .file ) ) {
                return( evaluateOneFile( .file, file.is.aa = TRUE ) );
            } else {
                return( c(
                    truths.perspective.HD.average.ignoringGaps = NA, truths.perspective.HD.sum.ignoringGaps = NA, truths.perspective.denominator.sum.ignoringGaps = NA, estimates.perspective.HD.average.ignoringGaps = NA, estimates.perspective.HD.sum.ignoringGaps = NA, estimates.perspective.denominator.sum.ignoringGaps = NA,
                    truths.perspective.HD.average.includingGaps = NA, truths.perspective.HD.sum.includingGaps = NA, truths.perspective.denominator.sum.includingGaps = NA, estimates.perspective.HD.average.includingGaps = NA, estimates.perspective.HD.sum.includingGaps = NA, estimates.perspective.denominator.sum.includingGaps = NA
                ) );
            }
        } );
    names( evaluate.results.by.protein.file ) <-
        gsub( ".FASTA", "", dir( proteins.dir, full.names = F ) );

    evaluate.results.by.nucleotide.file <- 
        lapply( dir( nucleotides.dir, full.names = T ), function ( .file ) {
            if( checkFastaFileIsReal( .file ) ) {
                return( evaluateOneFile( .file, file.is.aa = FALSE ) );
            } else {
                return( c(
                    truths.perspective.HD.average.ignoringGaps = NA, truths.perspective.HD.sum.ignoringGaps = NA, truths.perspective.denominator.sum.ignoringGaps = NA, estimates.perspective.HD.average.ignoringGaps = NA, estimates.perspective.HD.sum.ignoringGaps = NA, estimates.perspective.denominator.sum.ignoringGaps = NA,
                    truths.perspective.HD.average.includingGaps = NA, truths.perspective.HD.sum.includingGaps = NA, truths.perspective.denominator.sum.includingGaps = NA, estimates.perspective.HD.average.includingGaps = NA, estimates.perspective.HD.sum.includingGaps = NA, estimates.perspective.denominator.sum.includingGaps = NA
                ) );
            }
        } );
    names( evaluate.results.by.nucleotide.file ) <-
        gsub( ".FASTA", "", dir( nucleotides.dir, full.names = F ) );
    
    output.table.row.columns <- c( "estimates.file" = estimates.fasta.file.short, "truths.file" = truths.fasta.file.short, unlist( evaluate.results.by.protein.file ), unlist( evaluate.results.by.nucleotide.file ) );
    
    output.table.path <-
        paste( output.dir, "/", output.file, sep = "" );

    write.table( t( as.matrix( output.table.row.columns ) ), file = output.table.path, append = output.file.append, row.names = FALSE, col.names = ( !output.file.append || !file.exists( output.table.path ) ), sep = "\t", quote = FALSE );

    # Return the file name.
    return( output.table.path );
} # evaluateFounders ( estimates.fasta.file, ... )

## Here is where the action is.
estimates.fasta.file <- Sys.getenv( "evaluateFounders_estimatesFilename" );
truths.fasta.file <- Sys.getenv( "evaluateFounders_truthsFilename" );
output.table.file <- Sys.getenv( "evaluateFounders_outputFilename" );
if( nchar( output.table.file ) == 0 ) {
    output.table.file <- NULL;
}
output.dir <- Sys.getenv( "evaluateFounders_outputDir" );
if( nchar( output.dir ) == 0 ) {
    output.dir <- NULL;
}
append.to.output.file <- Sys.getenv( "evaluateFounders_append" );
if( ( nchar( append.to.output.file ) == 0 ) || ( append.to.output.file == "0" ) || ( toupper( append.to.output.file ) == "F" ) || ( toupper( append.to.output.file ) == "FALSE" ) ) {
    append.to.output.file <- FALSE;
} else {
    append.to.output.file <- TRUE;
}
genecutter.proteins.list <- Sys.getenv( "evaluateFounders_proteinsList" );
if( nchar( genecutter.proteins.list ) == 0 ) {
    genecutter.proteins.list <- "-GAG-POL-VIF-VPR-TAT-REV-VPU-ENV-NEF";
}
genecutter.region <- Sys.getenv( "evaluateFounders_genomeRegion" );
if( nchar( genecutter.region ) == 0 ) {
    genecutter.region <- "ALL";
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
        print( evaluateFounders( estimates.fasta.file, truths.fasta.file, output.dir = output.dir, output.file = output.table.file, output.file.append = append.to.output.file, genecutter.proteins.list = genecutter.proteins.list, genecutter.genome.region = genecutter.region ) );
    } else {
        stop( paste( "'truths' fasta file does not exist:", truths.fasta.file ) );
    }
} else {
    stop( paste( "'estimates' fasta file does not exist:", estimates.fasta.file ) );
}

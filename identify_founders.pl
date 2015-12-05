#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) identify_founders
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##  Description:
##
##      Script for using DiveIn tools (at the Mullins lab web site) to
##      estimate the number of founders and perform sequence
##      reconstruction of the founders.
##
##      Note that this creates output files in subdirectories named
##      after the input fasta file name (unless you specify an output dir).
##
##      R packages you'll need: (see installRPackages.R).
##      ade4
##      ape
##      dynamicTreeCut
##      entropy
##      ggplot2
##      binom
## source("https://bioconductor.org/biocLite.R")
##      biocLite("Biostrings")
##      biocLite("seqinr")
##
##      Try: rm -r rv217_1W_gold_standard-hiv-founder-id_-fr_resultDir/; mkdir rv217_1W_gold_standard-hiv-founder-id_-fr_resultDir/; perl -w ./identify_founders.pl -fr -O rv217_1W_gold_standard-hiv-founder-id_-fr_resultDir/ ~/src/from-git/projects/tholzman/MorgansFounderIDMethod/rv217_1W_gold_standard.list > rv217_1W_gold_standard-hiv-founder-id_-fr_resultDir/identify-founders.out
##      Or: rm -r caprisa002_1W_gold_standard-hiv-founder-id_-rP_resultDir/; mkdir caprisa002_1W_gold_standard-hiv-founder-id_-rP_resultDir/; perl ./identify_founders.pl -rP -O caprisa002_1W_gold_standard-hiv-founder-id_-rP_resultDir/ caprisa002_1W_gold_standard.list  > caprisa002_1W_gold_standard-hiv-founder-id_-rP_resultDir/identify-founders.out
##      This next one skips RAP because it seems to be exceptionally slow with these large numbers of sequences.  Unsure how to proceed - maybe iterately evaluate subsets? Or use a different program.
##      Or: mkdir caprisa002_1W_gold_standard-hiv-founder-id_resultDir/; perl ./identify_founders.pl -V -R -P -O caprisa002_1W_gold_standard-hiv-founder-id_resultDir/ caprisa002_1W_gold_standard.list  > caprisa002_1W_gold_standard-hiv-founder-id_resultDir/identify-founders.out
##      Or: mkdir CAPRISA002_ft_seqs-hiv-founder-id_resultDir/; perl ./identify_founders.pl -O CAPRISA002_ft_seqs-hiv-founder-id_resultDir/ ~/src/from-git/projects/tholzman/MorgansFounderIDMethod/CAPRISA002_ft_seqs.txt  > CAPRISA002_ft_seqs-hiv-founder-id_resultDir/identify-founders.out
##      Or: rm -r Abrahams-2009aa-hiv-founder-id_-fr_resultDir/; mkdir Abrahams-2009aa-hiv-founder-id_-fr_resultDir/; perl ./identify_founders.pl -fr -O Abrahams-2009aa-hiv-founder-id_-fr_resultDir/ Abrahams-2009aa/preparedFor_hiv-identify-founders.list > Abrahams-2009aa-hiv-founder-id_-fr_resultDir/identify-founders.out 
###******************************************************************************

use Getopt::Std; # for getopts
use File::Path qw( make_path );
use Path::Tiny;
require Sort::Fields; # for ??

use strict;
use vars qw( $opt_D $opt_V $opt_o $opt_O $opt_C $opt_P $opt_R $opt_F $opt_E $opt_H $opt_f $opt_w $opt_n $opt_I $opt_i $opt_v $opt_r );
use vars qw( $VERBOSE $DEBUG );

sub splitFastaFileOnHeaderPatterns {
## NOTE: If there is only one sequence in a group, we don't bother writing out the file.
  my $output_path_dir  = shift @_;
  my $input_fasta_file = shift @_;
  my @header_patterns  = @_;

  ## TODO: REMOVE/FIX MAGIC # $newline_column.
  our $newline_column = 40;
  ## TODO: REMOVE/FIX MAGIC # $minimum_sequences_for_output_file.
  our $minimum_sequences_for_output_file = 2;

  my ( $input_fasta_file_path, $input_fasta_file_short ) =
    ( $input_fasta_file =~ /^(.*?)\/([^\/]+)$/ );
  unless( $input_fasta_file_short ) {
    $input_fasta_file_short = $input_fasta_file;
    $input_fasta_file_path = ".";
  }
  my ( $input_fasta_file_short_nosuffix, $input_fasta_file_suffix ) =
    ( $input_fasta_file_short =~ /^([^\.]+)(\..+)?$/ );

  if( $VERBOSE ) { print "Opening input Fasta file \"$input_fasta_file\".."; }
  open( INPUT_FASTA_FH, $input_fasta_file ) or
    die "Unable to open fasta file \"$input_fasta_file\": $!";
  if( $VERBOSE ) { print ".done.\n"; }

  if( $VERBOSE ) { print "Parsing Fasta file.."; }

  my $header_patterns_together = join( "|", @header_patterns );
  my $sequence;
  my $sequence_header;
  my $have_a_match = 0;
  my %sequences_by_pattern;
  while( <INPUT_FASTA_FH> ) {
    if( /^>/ ) {
      if( defined( $sequence ) ) {
        foreach my $header_pattern ( @header_patterns ) {
          my $header_pattern_copy = $header_pattern;
          $header_pattern_copy =~ s/\|/\\\|/; # escape any bars
          if( $sequence_header =~ m/$header_pattern_copy/ ) {
            if( !exists( $sequences_by_pattern{ $header_pattern } ) ) {
              $sequences_by_pattern{ $header_pattern } = {};
            }
            if( exists( $sequences_by_pattern{ $header_pattern }{ $sequence_header } ) ) {
              die( "Found multiple sequences named '$sequence_header'!" );
            }
            $sequences_by_pattern{ $header_pattern }{ $sequence_header } = $sequence;
            last; # So it can't match more than one pattern.
          }
        } # End foreach $header_pattern
        undef $sequence;
      }
      $sequence_header = $_;
      $have_a_match = ( $sequence_header =~ m/$header_patterns_together/ );
      if( $DEBUG ) {
        print "Found sequence header: $sequence_header is a match to a pattern:";
        foreach my $header_pattern ( @header_patterns ) {
          my $header_pattern_copy = $header_pattern;
          $header_pattern_copy =~ s/\|/\\\|/; # escape any bars
          if( $sequence_header =~ m/$header_pattern_copy/ ) {
            print "$header_pattern_copy\n";
            last; # So it can't match more than one pattern.
          }
        }
      }
    } elsif( $have_a_match ) {
      # Get rid of all preexisting whitespace
      $_ =~ s/\s//g;
  
      $sequence .= $_;
    }
  } # End foreach $_ <input_fasta_fileFH>
      if( defined( $sequence ) ) {
        foreach my $header_pattern ( @header_patterns ) {
          my $header_pattern_copy = $header_pattern;
          $header_pattern_copy =~ s/\|/\\\|/; # escape any bars
          if( $sequence_header =~ m/$header_pattern_copy/ ) {
            if( !exists( $sequences_by_pattern{ $header_pattern } ) ) {
              $sequences_by_pattern{ $header_pattern } = {};
            }
            if( exists( $sequences_by_pattern{ $header_pattern }{ $sequence_header } ) ) {
              die( "Found multiple sequences named '$sequence_header'!" );
            }
            $sequences_by_pattern{ $header_pattern }{ $sequence_header } = $sequence;
            last; # So it can't match more than one pattern.
          }
        } # End foreach $header_pattern
        undef $sequence;
      }
  if( $VERBOSE ) { print ".done.\n"; }

  if( $VERBOSE  ) { print "Closing file \"$input_fasta_file\".."; }
  close INPUT_FASTA_FH;
  if( $VERBOSE ) { print ".done.\n"; }

  sub print_sequence {
    my $OUTFH = shift;
    my $sequence = shift;

    my $seq_len = length( $sequence );
    for( my $seq_pos = 0; $seq_pos < $seq_len; $seq_pos += $newline_column ) {
      if( ( $seq_pos + $newline_column ) > $seq_len ) {
        print $OUTFH substr( $sequence, $seq_pos, ( $seq_len - $seq_pos ) ), "\n\n";
      } else {
        print $OUTFH substr( $sequence, $seq_pos, $newline_column ), "\n";
      }
    }
  } # print_sequence (..)

  my @output_files;
  foreach my $header_pattern ( @header_patterns ) {
    if( $DEBUG ) {
      print "--- $header_pattern ---\n";
    }

    next unless ( exists( $sequences_by_pattern{ $header_pattern } ) );

    if( scalar( keys %{ $sequences_by_pattern{ $header_pattern } } ) < $minimum_sequences_for_output_file ) {
      if( $VERBOSE ) {
        print "Skipping $header_pattern because only ", scalar( keys %{ $sequences_by_pattern{ $header_pattern } } ), " sequence matches it, which is not enough (we require a minimum of $minimum_sequences_for_output_file).\n";
      }
      next;
    }

    my $header_pattern_copy = $header_pattern;
    $header_pattern_copy =~ s/\|//; # remove any bars from the pattern

    my $output_file = $output_path_dir . '/' . $input_fasta_file_short_nosuffix . '_' . $header_pattern_copy . $input_fasta_file_suffix;

    if( $VERBOSE ) { print "Opening file \"$output_file\" for writing.."; }
    my $OUTPUT_FH;
    unless( open $OUTPUT_FH, ">$output_file" ) {
      warn "Unable to open output file \"$output_file\": $!";
      return 1;
    }
    if( $VERBOSE ) { print ".done.\n"; }

    if( $VERBOSE ) { print "Printing sequences with headers matching pattern \"$header_pattern\""; }
    foreach my $sequence_header ( keys %{ $sequences_by_pattern{ $header_pattern } } ) {
      if( $VERBOSE ) { print "."; }
      print $OUTPUT_FH $sequence_header; # newline is still there.
      print_sequence( $OUTPUT_FH, $sequences_by_pattern{ $header_pattern }{ $sequence_header } );
    }
    if( $VERBOSE ) { print ".done.\n"; }

    if( $VERBOSE && $output_file ) { print "Closing file \"$output_file\".."; }
    close $OUTPUT_FH;
    if( $VERBOSE && $output_file ) { print ".done.\n"; }

    push @output_files, $output_file;
  } # End foreach $header_pattern
         
  return( @output_files );
} # sub splitFastaFileOnHeaderPatterns (..)

sub identify_founders {
  @ARGV = @_;

  sub identify_founders_usage {
    print "\tidentify_founders [-DV] [-CPRFEHfnIsr] [-i <insites_threshold>] [-v <insites_threshold_v3>] [-(o|O) <output_dir>] <input_fasta_file1 or file_listing_input_files> [<input_fasta_file2> ...] \n";
    exit;
  }

  # This means -D, -o, -O, -V, ... are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_o is an optional directory to put the output; default depends on input filename.
  # opt_O is just like opt_o.  Same thing.
  # opt_V means be verbose.
  # opt_C means skip (do not run) the PhyML diversity or InSites analysis and the primary founder call (which depends on those)
  # opt_P means skip (do not run) profillic
  # opt_R means skip (do not run) RAP (alignment-internal recombination detection program)
  # opt_F means skip (do not run) PFitter
  # opt_E means skip (do not run) entropy statistics calculation
  # opt_H means skip (do not run) HYPERMUT 2.0 hypermutation detection
  # opt_f means fix hypermutated sequences, instead of removing them.
  # opt_w is what we should fix hypermutated sequences with (default: R)
  # opt_n means mask out nonsynonymous codons (compared to input file consensus) when running PFitter.
  # opt_I means run inSites online (instead of offline)
  # opt_i is the insites threshold to use for all regions except v3 (default: 0.85).
  # opt_v is the insites threshold to use for v3 (default: 0.33).
  # opt_r means recursively operate on clusters identified using the Informtive Sites method.
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_o, $opt_O, $opt_C, $opt_P, $opt_R, $opt_F, $opt_E, $opt_H, $opt_f, $opt_w, $opt_n, $opt_I, $opt_i, $opt_v, $opt_s, $opt_r ) = ();
  if( not getopts('DVo:O:CPRFEHfw:nIi:v:sr') ) {
    identify_founders_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;

  my $run_InSites_and_PhyML = !$opt_C;
  my $run_profillic = !$opt_P;
  my $run_RAP = !$opt_R;
  my $run_PFitter = !$opt_F;
  my $run_entropy = !$opt_E;
  my $run_Hypermut = !$opt_H;
  my $fix_hypermutated_sequences = $opt_f || 0;
  my $fix_hypermutated_sequences_with = $opt_w || "R";
  my $mask_out_nonsynonymous_codons_in_PFitter = $opt_n || 0;
  my $in_sites_ratio_threshold = $opt_i || 0.85; # For Abrahams and RV217
  my $in_sites_ratio_threshold_v3 = $opt_v || 0.33; # For caprisa002
  my $runInSites_online = $opt_I || 0;
  my $recurse_on_clusters = $opt_r || 0;

  ## TODO: make into an arg
  my $run_DSPFitter = 1;
  
  ## TODO: DEHACKIFY MAGIC #s
  my $mean_diversity_threshold = 0.001;
  my $RAP_pValueThreshold = 0.0007; # Appears to be the suggestion from the output file "(summaryTable)"'s column header, which reads "Pvalues<0.0007".
  my $hypermut2_pValueThreshold = 0.1; # Matches Abrahams 2009

  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }
  
  my @fasta_files_remaining_to_process = @ARGV;
  my @all_input_fasta_files = @fasta_files_remaining_to_process;
  unless( scalar @fasta_files_remaining_to_process ) { identify_founders_usage(); }

  # Special case: if there is only one file, it might be a textfile
  # containing names of files.
  if( scalar( @fasta_files_remaining_to_process ) == 1 ) {
    #print $fasta_files_remaining_to_process[ 0 ], "\n";
    unless( ( $fasta_files_remaining_to_process[ 0 ] ) =~ /\.fast?a?$/ ) {
      # Try opening it, see if it's a list of files.
      if( $VERBOSE ) {
        print "Reading file names from file \"", $fasta_files_remaining_to_process[ 0 ], "\"..";
      }
      my $input_fasta_file_contents = path( $fasta_files_remaining_to_process[ 0 ] )->slurp();
      if( $VERBOSE ) {
        print ".done\n";
      }
      if( $DEBUG ) {
        print $input_fasta_file_contents;
      }
      @fasta_files_remaining_to_process = split( "\n", $input_fasta_file_contents );
      @all_input_fasta_files = @fasta_files_remaining_to_process;
    }
  }

  my $output_path_dir = $opt_o || $opt_O || undef;
  # Remove the trailing "/" if any
  if( defined( $output_path_dir ) ) {
    ( $output_path_dir ) = ( $output_path_dir =~ /^(.+)\/*$/ );
  }
  if( defined $output_path_dir ) {
    make_path( $output_path_dir );
  } else {
    $output_path_dir = ".";
  }

  my $extra_flags = "";
  if( $DEBUG ) {
    $extra_flags .= "-D ";
  }
  if( $VERBOSE ) {
    $extra_flags .= "-V ";
  }

  my $output_table_file = $output_path_dir . '/' . "identify_founders.tab";

  if( $VERBOSE ) { print "Opening file \"$output_table_file\" for writing.."; }
  unless( open OUTPUT_TABLE_FH, ">$output_table_file" ) {
    warn "Unable to open output file \"$output_table_file\": $!";
    return 1;
  }
  # Turn on autoflush for this too.
  my $old_fh = select( OUTPUT_TABLE_FH );
  $| = 1;
  select( $old_fh );
  if( $VERBOSE ) { print ".done.\n"; }

  if( $VERBOSE ) { print "Writing results table to file \"$output_table_file\".."; }
  my @table_column_headers =
    (
     "infile"
    );
  if( $run_Hypermut ) {
    if( $fix_hypermutated_sequences ) {
      push @table_column_headers, "fixed-hypermut";
    } else {
      push @table_column_headers, "removed-hypermut";
    }
  }
  if( $run_RAP ) {
    push @table_column_headers, "removed-recomb";
  }

  push @table_column_headers,
  (
   "file"
   );
  if( $run_InSites_and_PhyML ) {
    push @table_column_headers,
      (
       "num.seqs", "num.diversity.seqs", "diversity", "inf.to.priv.ratio",
       "exceeds.diversity.threshold", "exceeds.ratio.threshold", "ratio.threshold",
       "is.one.founder"
      );
  }

  if( $run_entropy ) {
    push @table_column_headers,
    (
     "mean.entropy", "sd.entropy"
    );
  }
  
  if( $run_PFitter ) {
    push @table_column_headers,
    (
     "PFitter.lambda",
     "PFitter.se",
     "PFitter.nseq",
     "PFitter.nbases",
     "PFitter.mean.hd",
     "PFitter.max.hd",
     "PFitter.time.est",
     "PFitter.time.ci.low",
     "PFitter.time.ci.high",
     "PFitter.chi.sq.stat",
     "PFitter.chi.sq.df",
     "PFitter.chi.sq.p.value",
     "is.poisson",
     "is.starlike",
     "Bayesian.PFitter.lambda.est",
     "Bayesian.PFitter.lambda.ci.low",
     "Bayesian.PFitter.lambda.ci.high",
     "Bayesian.PFitter.days.est",
     "Bayesian.PFitter.days.ci.low",
     "Bayesian.PFitter.days.ci.high",
     "DS.PFitter.lambda.est",
     "DS.PFitter.lambda.ci.low",
     "DS.PFitter.lambda.ci.high",
     "DS.PFitter.days.est",
     "DS.PFitter.days.ci.low",
     "DS.PFitter.days.ci.high",
     "DS.PFitter.distance.est",
     "DS.PFitter.distance.ci.low",
     "DS.PFitter.distance.ci.high",
     "DS.PFitter.assertion.low",
     "DS.PFitter.assertion.high",
     "DS.PFitter.fits",
     "DS.PFitter.R",
     "DS.PFitter.is.starlike",
     "DS.PFitter.starlike.pvalue",
     "DS.PFitter.lower.is.starlike",
     "DS.PFitter.lower.starlike.pvalue",
     "DS.PFitter.upper.is.starlike",
     "DS.PFitter.upper.starlike.pvalue",
     "is.one.founder.alt"
    );
    
    push @table_column_headers,
    (
     "multifounder.PFitter.lambda",
     "multifounder.PFitter.se",
     "multifounder.PFitter.nseq",
     "multifounder.PFitter.nbases",
     "multifounder.PFitter.mean.hd",
     "multifounder.PFitter.max.hd",
     "multifounder.PFitter.time.est",
     "multifounder.PFitter.time.ci.low",
     "multifounder.PFitter.time.ci.high",
     "multifounder.PFitter.chi.sq.stat",
     "multifounder.PFitter.chi.sq.df",
     "multifounder.PFitter.chi.sq.p.value",
     "multifounder.is.poisson",
#     "multifounder.is.starlike",
     "multifounder.Bayesian.PFitter.lambda.est",
     "multifounder.Bayesian.PFitter.lambda.ci.low",
     "multifounder.Bayesian.PFitter.lambda.ci.high",
     "multifounder.Bayesian.PFitter.days.est",
     "multifounder.Bayesian.PFitter.days.ci.low",
     "multifounder.Bayesian.PFitter.days.ci.high",
     "multifounder.DS.PFitter.lambda.est",
     "multifounder.DS.PFitter.lambda.ci.low",
     "multifounder.DS.PFitter.lambda.ci.high",
     "multifounder.DS.PFitter.days.est",
     "multifounder.DS.PFitter.days.ci.low",
     "multifounder.DS.PFitter.days.ci.high",
     "multifounder.DS.PFitter.distance.est",
     "multifounder.DS.PFitter.distance.ci.low",
     "multifounder.DS.PFitter.distance.ci.high",
     "multifounder.DS.PFitter.assertion.low",
     "multifounder.DS.PFitter.assertion.high",
     "multifounder.DS.PFitter.fits",
     "multifounder.DS.PFitter.R",
     "multifounder.DS.PFitter.is.starlike",
     "multifounder.DS.PFitter.starlike.pvalue",
     "multifounder.DS.PFitter.lower.is.starlike",
     "multifounder.DS.PFitter.lower.starlike.pvalue",
     "multifounder.DS.PFitter.upper.is.starlike",
     "multifounder.DS.PFitter.upper.starlike.pvalue"
    );
  } # End if $run_PFitter

  if( $run_profillic ) {
    push @table_column_headers, "profillic.clusters", "profillic.founders";
  } # End if $run_profillic

  if( $run_InSites_and_PhyML || $run_PFitter ) {
    push @table_column_headers, "insites.clusters";
  }
  if( $run_InSites_and_PhyML ) {
    push @table_column_headers, "insites.founders";
  }
  if( $run_PFitter ) {
    push @table_column_headers, "starphy.founders";
  }
  
  my $table_header = join( "\t", @table_column_headers );
  print OUTPUT_TABLE_FH $table_header, "\n";

  # After removing/fixing hypermutated sequences and removing recombined sequences, run InSites (online) and get the informative:private site ratio stat; also run PhyML (locally) to get the mean pairwise diversity statistic, and also cluster the informative sites subalignments and compute the full consensus sequences for each cluster.
  my $id_string = "";
  my $output_path_dir_for_input_fasta_file;
  my $R_output;
  my $input_fasta_file;
  my $fasta_file;
  my $original_short_name_nosuffix_stripped;
  my %final_input_fasta_file_short_names_by_original_short_name_stripped;
  while ( scalar( @fasta_files_remaining_to_process ) > 0 ) {
    $input_fasta_file = shift @fasta_files_remaining_to_process;
    $fasta_file = $input_fasta_file;
    if( $VERBOSE ) {
      print $fasta_file, "\n";
    }
    my ( $fasta_file_path, $fasta_file_short ) =
      ( $fasta_file =~ /^(.*?)\/([^\/]+)$/ );
    unless( $fasta_file_short ) {
      $fasta_file_short = $fasta_file;
      $fasta_file_path = ".";
    }
    my ( $fasta_file_short_nosuffix, $fasta_file_suffix ) =
      ( $fasta_file_short =~ /^([^\.]+)(\..+)?$/ );

    if( defined $output_path_dir ) {
      $output_path_dir_for_input_fasta_file = $output_path_dir;
    } else {
      $output_path_dir_for_input_fasta_file =
        $fasta_file_path . "/" . $fasta_file_short_nosuffix . "_hiv-founder-id_resultsDir";
    }
    # remove trailing "/"
    ( $output_path_dir_for_input_fasta_file ) = ( $output_path_dir_for_input_fasta_file =~ /^(.*[^\/])\/*$/ );

    if( $fasta_file_short_nosuffix =~ /^(.+)(?:(?:_[LR]H)|(?:_NFLG)|(?:_env))/ ) {
      ( $original_short_name_nosuffix_stripped ) =
        ( $fasta_file_short_nosuffix =~ /^(.+)(?:(?:_[LR]H)|(?:_NFLG)|(?:_env))/ );
    } else {
      $original_short_name_nosuffix_stripped =
        $fasta_file_short_nosuffix;
    }

    # Duplicate the input file, possibly modifying it.
    my $fasta_file_contents = path( $fasta_file )->slurp();
    my ( @seq_headers ) = ( $fasta_file_contents =~ /^>(.*)$/mg );
    if( scalar( @seq_headers ) == 1 ) {
      # THE INPUT FASTA FILE CONTAINS ONLY ONE SEQUENCE.  SKIP IT.
      #if( $VERBOSE ) {
        print( "\nThere's only one sequence in input file \"$fasta_file\".  Skipping it.\n" );
      #}
      next;
    }
    ## HACK: make sure there are no bangs in the input file (since sometimes there are, right now).
    if( $fasta_file_contents =~ /\!/ ) {
      if( $VERBOSE ) {
        print( "Input file \"$fasta_file\" contains illegal characters \"!\"; changing them to gaps, saving to output directory.\n" );
        ## TODOL REMOVE
        #print( $output_path_dir_for_input_fasta_file );
      }
      $fasta_file_contents =~ s/\!/-/g;
    } else {
      if( $VERBOSE ) {
        print( "Copying input file \"$fasta_file\" to output directory.\n" );
      }
    }
    # Now write it out to a temporary location in the output dir.
    $fasta_file_path = $output_path_dir_for_input_fasta_file;
    $fasta_file = "${fasta_file_path}/${fasta_file_short}";
    if( $VERBOSE ) {
      if( $fasta_file_contents =~ /\!/ ) {
        print( "Writing out fixed input file \"$fasta_file\".." );
      } else {
        print( "Writing out copy of input file to \"$fasta_file\".." );
      }
    }
    if( $VERBOSE ) { print "Opening file \"$fasta_file\" for writing..\n"; }
    unless( open input_fasta_fileFH, ">$fasta_file" ) {
      warn "Unable to open output file \"$fasta_file\": $!\n";
      return 1;
    }
    ( @seq_headers ) = map { s/\!/-/g } @seq_headers;
    print input_fasta_fileFH $fasta_file_contents;
    close( input_fasta_fileFH );
    
    if( ( $fasta_file_contents =~ /RH\|/ ) && ( $fasta_file_contents =~ /LH\|/ ) ) {
      ## SPECIAL CASE: If the input file contains left-half and right-half genomes, separate them into separate input files and do both.
      if( $VERBOSE ) {
        print( "Input file contains left- and right- halves of the genome; separating them.\n" );
      }
      my @new_input_fasta_files =
        splitFastaFileOnHeaderPatterns( $output_path_dir_for_input_fasta_file, $fasta_file, , "env\|", "NFLG\|", "LH\|", "RH\|" );
      if( $VERBOSE ) {
        print( "Queueing the new files ( ", join( ", ", @new_input_fasta_files ), " ).\n" );
      }
      push @all_input_fasta_files, @new_input_fasta_files;
      unshift @fasta_files_remaining_to_process, @new_input_fasta_files;
      next; # That's all we do with the original input file.
    }

    # This is the original unaltered fasta file (path removed), printed to the table out:
    print OUTPUT_TABLE_FH $fasta_file_short;
  
    print "\nInput Fasta file: $fasta_file_short\n";

    my $fasta_file_mask_out_nonsynonymous_codons_in_PFitter = $mask_out_nonsynonymous_codons_in_PFitter;
    if( $fasta_file_short_nosuffix =~ /^(.+)_NFLG/ ) {
      # NFLG seqs are not in one reading frame, so codons don't translate, so we don't do it.
      $fasta_file_mask_out_nonsynonymous_codons_in_PFitter = 0;
    }
    
    # First fix/remove hypermutated sequences, using an implementation of the HYPERMUT 2.0 algorithm.
    if( $run_Hypermut ) {
      if( $VERBOSE ) {
          if( $fix_hypermutated_sequences ) {
            print "Calling R to fix hypermutated sequences (with $fix_hypermutated_sequences_with)..";
          } else {
            print "Calling R to remove hypermutated sequences..";
          }
      }
      $R_output = `export removeHypermutatedSequences_fixWith="$fix_hypermutated_sequences_with"; export removeHypermutatedSequences_fixInsteadOfRemove="$fix_hypermutated_sequences"; export removeHypermutatedSequences_pValueThreshold="$hypermut2_pValueThreshold"; export removeHypermutatedSequences_inputFilename="$fasta_file"; export removeHypermutatedSequences_outputDir="$output_path_dir_for_input_fasta_file"; R -f removeHypermutatedSequences.R --vanilla --slave`;
      if( $VERBOSE ) {
        print( $R_output );
      }
      ## extract the number fixed/removed from the output
      my ( $num_hypermut_sequences ) = ( $R_output =~ /^\[1\]\s*(\d+)\s*$/m );
  #    if( $VERBOSE ) {
          if( $fix_hypermutated_sequences ) {
            print( "The number of hypermutated sequences fixed is: $num_hypermut_sequences\n" );
          } else {
            print( "The number of hypermutated sequences removed is: $num_hypermut_sequences\n" );
          }
  #    }

      print OUTPUT_TABLE_FH "\t", $num_hypermut_sequences;

      # Now use the output from that..
      if( $num_hypermut_sequences > 0 ) {
        $fasta_file_path = $output_path_dir_for_input_fasta_file;
        if( $fix_hypermutated_sequences ) {
          $fasta_file_short = "${fasta_file_short_nosuffix}_fixHypermutatedSequencesWith${fix_hypermutated_sequences_with}${fasta_file_suffix}";
        } else {
          $fasta_file_short = "${fasta_file_short_nosuffix}_removeHypermutatedSequences${fasta_file_suffix}";
          # Recompute the seq_headers.
          $fasta_file_contents = path( "${fasta_file_path}/${fasta_file_short}" )->slurp();
          ( @seq_headers ) = ( $fasta_file_contents =~ /^>(.*)$/mg );
        }
        ( $fasta_file_short_nosuffix, $fasta_file_suffix ) =
          ( $fasta_file_short =~ /^([^\.]+)(\..+)?$/ );
        $fasta_file = "${fasta_file_path}/${fasta_file_short}";
      } # End if we need to change files..
      
      if( $VERBOSE ) {
        print ".done.\n";
      }
    } # End if $run_Hypermut

    if( $run_RAP ) {
      ## Run RAP on LANL, which gives individual
      ## sequences that are recombinants of other individual sequences,
      ## allowing those to be flagged for removal just like hypermutated
      ## ones.
      if( $VERBOSE ) {
        print "Running RAP at LANL to compute recombined sequences..";
      }
      my $removed_recombined_sequences = 0;
      my $RAP_result_stdout = `perl runRAPOnline.pl $extra_flags $fasta_file $output_path_dir_for_input_fasta_file`;
      if( $VERBOSE ) {
        print $RAP_result_stdout;
      }
      my ( $removeDuplicateSequencesFromAlignedFasta_output_file ) = ( $RAP_result_stdout =~ /^Table of duplicates removed: (.+)\s*$/m );
      if( !defined( $removeDuplicateSequencesFromAlignedFasta_output_file ) ) { # Do duplicates.
        $removeDuplicateSequencesFromAlignedFasta_output_file = "";
      }
      ## TODO: REMOVE
      #print( "\$removeDuplicateSequencesFromAlignedFasta_output_file is $removeDuplicateSequencesFromAlignedFasta_output_file" );
      if( $RAP_result_stdout =~ /Recombinants identified/ ) {
        my ( $RAP_output_file ) = ( $RAP_result_stdout =~ /Recombinants identified \((.+)\)\s*$/ );
        $R_output = `export removeRecombinedSequences_pValueThreshold="$RAP_pValueThreshold"; export removeRecombinedSequences_RAPOutputFile="$RAP_output_file"; export removeRecombinedSequences_removeDuplicateSequencesFromAlignedFastaOutputFile="$removeDuplicateSequencesFromAlignedFasta_output_file"; export removeRecombinedSequences_inputFilename="$fasta_file"; export removeRecombinedSequences_outputDir="$output_path_dir_for_input_fasta_file"; R -f removeRecombinedSequences.R --vanilla --slave`;
        if( $VERBOSE ) {
          print $R_output;
        }
        ## extract the number fixed/removed from the output
        ( $removed_recombined_sequences ) = ( $R_output =~ /^.*\[1\]\s*(\d+)\s*$/m );
        if( $VERBOSE ) {
          print( ".done." );
        }
        # Now use the output from that..
        $fasta_file_path = $output_path_dir_for_input_fasta_file;
        $fasta_file_short = "${fasta_file_short_nosuffix}_removeRecombinedSequences${fasta_file_suffix}";
        ( $fasta_file_short_nosuffix, $fasta_file_suffix ) =
          ( $fasta_file_short =~ /^([^\.]+)(\..+)?$/ );
        $fasta_file = "${fasta_file_path}/${fasta_file_short}";
        # Recompute the seq_headers.
        $fasta_file_contents = path( $fasta_file )->slurp();
        ( @seq_headers ) = ( $fasta_file_contents =~ /^>(.*)$/mg );
      } else {
        if( $VERBOSE ) {
          print( ".done." );
        }
        $removed_recombined_sequences = 0;
      } # End if any recombinants were identified. .. else ..
      print( "The number of recombined sequences removed is: $removed_recombined_sequences\n" );

      print OUTPUT_TABLE_FH "\t", $removed_recombined_sequences;
    } # End if $run_RAP

    print "Fasta file for analysis: $fasta_file_short\n";
    # This is the file for analysis (again, path stripped); maybe altered from the original input file.
    print OUTPUT_TABLE_FH "\t", $fasta_file_short;

    ## TODO: REMOVE
    #print( "SETTING \$final_input_fasta_file_short_names_by_original_short_name_stripped{ $original_short_name_nosuffix_stripped } = $fasta_file_short_nosuffix\n" );
    if( exists( $final_input_fasta_file_short_names_by_original_short_name_stripped{ $original_short_name_nosuffix_stripped } ) ) {
      push @{ $final_input_fasta_file_short_names_by_original_short_name_stripped{ $original_short_name_nosuffix_stripped } }, $fasta_file_short_nosuffix;
    } else {
      $final_input_fasta_file_short_names_by_original_short_name_stripped{ $original_short_name_nosuffix_stripped } = [ $fasta_file_short_nosuffix ];
    }

    my $morgane_calls_one_cluster = 1;
    my $diversity_threshold_exceeded = 0;
    my $in_sites_ratio_threshold_exceeded = 0;
    my $in_sites_founders_call = 1;
    my $ds_pfitter_founders_call = 1;
    my $num_phyml_seqs = undef;
    my $mean_diversity = undef;
    if( $run_InSites_and_PhyML ) {
      my $runInSites_output;
      if( $runInSites_online ) {
        $runInSites_output = `perl runInSitesOnline.pl $extra_flags $fasta_file $output_path_dir_for_input_fasta_file`;
      } else {
        $runInSites_output = `perl runInSitesOffline.pl $extra_flags $fasta_file $output_path_dir_for_input_fasta_file`;
      }
      if( $VERBOSE ) {
        print $runInSites_output;
      }
      my $getInSitesStat_output = `perl getInSitesStat.pl $extra_flags ${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}_informativeSites.txt ${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}_privateSites.txt ${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}_inSitesRatioStat.txt`;
      if( $VERBOSE ) {
        print $getInSitesStat_output;
      }
      my $in_sites_ratio = `cat ${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}_inSitesRatioStat.txt`;
      ( $in_sites_ratio ) = ( $in_sites_ratio =~ /^\s*(\S+)\s*$/ ); # Strip trailing newline, and any other flanking whitespace.
      
      # Run phyML, get stats
      my $phyml_out = `perl runPhyML.pl $extra_flags ${fasta_file} ${output_path_dir_for_input_fasta_file}`;
      if( $VERBOSE ) {
        print $phyml_out;
      }
      my $pairwise_diversity_stats = `cat ${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}.phylip_phyml_pwdiversity.txt`;
      ( $num_phyml_seqs, $mean_diversity ) =
        ( $pairwise_diversity_stats =~ /Max\s+(\d+)\s+([e\-\.\d]+)\s+/ );
      unless( defined( $num_phyml_seqs ) ) {
        warn( "UH OH: $fasta_file\nGOT:\n$pairwise_diversity_stats\n" );
        $num_phyml_seqs = 0;
        $mean_diversity = 0;
      }
      # unless( $num_phyml_seqs == scalar( @seq_headers ) ) {
      #   ## THIS IS BECAUSE phyml apparently only counts unique sequences.
      #   # print "WHY ARE THERE $num_seqs SEQS, when there are ", scalar( @seq_headers ), " sequences in the input file?\n";
      # }
      my $num_seqs = scalar( @seq_headers );
      print "Number of sequences: $num_seqs\n";
      if( $num_seqs != $num_phyml_seqs ) { print "BUT NOTE: Number of PhyML sequences: $num_phyml_seqs\n" };
      print "Mean pairwise diversity: $mean_diversity\n";
      print "Informative sites to private sites ratio: $in_sites_ratio\n"; # Newline is no longer on there
      
      #print "ABOUT TO WRITE OUT \$num_seqs ($num_seqs) and \$num_phyml_seqs ($num_phyml_seqs)\n";
      printf OUTPUT_TABLE_FH "\t%d\t%d\t%1.4f\t%1.4f", ( $num_seqs, $num_phyml_seqs, $mean_diversity, $in_sites_ratio );
      #print "did it\n";
      
      my $the_relevant_in_sites_ratio_threshold = $in_sites_ratio_threshold;
      if( $fasta_file_short_nosuffix =~ /^(.+)(?:_v3)/ ) {
        if( $VERBOSE ) {
          print "Using v3-specific InSites ratio threshold: $in_sites_ratio_threshold_v3\n";
        }
        $the_relevant_in_sites_ratio_threshold = $in_sites_ratio_threshold_v3;
      }
      
      # Now cluster the informative sites (only relevant if one or both of the above exceeds a threshold.
      if( $mean_diversity > $mean_diversity_threshold ) {
        $diversity_threshold_exceeded = 1;
        if( $in_sites_ratio > $the_relevant_in_sites_ratio_threshold ) {
          $in_sites_ratio_threshold_exceeded = 1;
          print( "Diversity threshold ($mean_diversity_threshold) exceeded\n" );
          print( "Ratio threshold ($the_relevant_in_sites_ratio_threshold) exceeded too\n" );
          $morgane_calls_one_cluster = 0;
        } else {
          print( "Diversity threshold ($mean_diversity_threshold) exceeded\n" );
        }
      } elsif( $in_sites_ratio > $the_relevant_in_sites_ratio_threshold ) {
          $in_sites_ratio_threshold_exceeded = 1;
          print( "Ratio threshold ($the_relevant_in_sites_ratio_threshold) exceeded\n" );
      }
      print OUTPUT_TABLE_FH "\t", $diversity_threshold_exceeded;
      print OUTPUT_TABLE_FH "\t", $in_sites_ratio_threshold_exceeded;
      print OUTPUT_TABLE_FH "\t", $the_relevant_in_sites_ratio_threshold;
      ## "is.one.founder":
      print OUTPUT_TABLE_FH "\t", $morgane_calls_one_cluster;
      
      if( $morgane_calls_one_cluster ) {
        print "Number of founders estimated using informative:private sites ratio and diversity thresholding is: 1\n";
        $in_sites_founders_call = 1; # go with Morgane's call if it is 1.
      } else {
        print "Number of founders estimated using informative:private sites ratio and diversity thresholding is: greater than 1\n";
      }
    } # End if $run_InSites_and_PhyML

    if( $run_entropy ) {
      $R_output = `export computeEntropyFromAlignedFasta_inputFilename="$fasta_file"; export computeEntropyFromAlignedFasta_outputDir="$output_path_dir_for_input_fasta_file"; R -f computeEntropyFromAlignedFasta.R --vanilla --slave`;
      if( $VERBOSE ) {
        print( $R_output );
      }
      my ( $entropy_stats_filename ) = ( $R_output =~ /^\[1\]\s*\"(.+)\"\s*$/m );
      if( !-e $entropy_stats_filename ) {
        warning( "UH OH: MISSING ENTROPY RESULT FILE \"$entropy_stats_filename\"\n" );
        #stopifnot( -e $entropy_stats_filename );
      }
      my $entropy_stats = path( $entropy_stats_filename )->slurp();

      my ( $entropy_seqs, $entropy_sites, $entropy_min, $entropy_q1, $entropy_median, $mean_entropy, $entropy_q3, $entropy_max, $sd_entropy ) =
        ( $entropy_stats =~ /\"SD\"\s+\"(\d*)\"\s+\"(\d*)\"\s+\"([^\"]*)\"\s+\"([^\"]*)\"\s+\"([^\"]*)\"\s+\"([^\"]*)\"\s+\"([^\"]*)\"\s+\"([^\"]*)\"\s+\"([^\"]*)\"\s*$/ );
      my $entropy_IQR = sprintf( "%0.4f", ( $entropy_q3 - $entropy_q1 ) );
      if( $VERBOSE ) {
        print "Entropy: Mean = $mean_entropy (SD = $sd_entropy)\n";
        #print "Entropy: Median = $entropy_median (IQR = $entropy_IQR)\n";
      }
      print OUTPUT_TABLE_FH "\t", $mean_entropy;
      print OUTPUT_TABLE_FH "\t", $sd_entropy;
    } # End if $run_entropy
    
    ## Now run PoissonFitter.
    my $PFitter_lambda = 0;
    my $PFitter_se = 0;
    my $PFitter_nseq = 0;
    my $PFitter_nbases = 0;
    my $PFitter_mean_hd = 0;
    my $PFitter_max_hd = 0;
    my $PFitter_chi_sq_stat = 0;
    my $PFitter_chi_sq_df = 0;
    my $PFitter_chi_sq_p_value = 1;
    my $PFitter_days_est_and_ci = "0 (0,0)";
    my $PFitter_days_est = 0;
    my $PFitter_days_ci_low = 0;
    my $PFitter_days_ci_high = 0;
    my $is_poisson = 1;
    my $is_starlike = 1;
    my $Bayesian_PFitter_lambda_est = 0;
    my $Bayesian_PFitter_lambda_ci_low = 0;
    my $Bayesian_PFitter_lambda_ci_high = 0;
    my $Bayesian_PFitter_days_est = 0;
    my $Bayesian_PFitter_days_ci_low = 0;
    my $Bayesian_PFitter_days_ci_high = 0;
    my $DS_PFitter_lambda_est = 0;
    my $DS_PFitter_lambda_ci_low = 0;
    my $DS_PFitter_lambda_ci_high = 0;
    my $DS_PFitter_days_est = 0;
    my $DS_PFitter_days_ci_low = 0;
    my $DS_PFitter_days_ci_high = 0;
    my $DS_PFitter_distance_mean = 0;
    my $DS_PFitter_distance_ci_low = 0;
    my $DS_PFitter_distance_ci_high = 0;
    my $DS_PFitter_fitstext = "DSPFitter test that intersequence rate = 2 x seq-consensus rate: OK";
    my $DS_PFitter_fits = 1;
    my $DS_PFitter_is_starlike = 1;
    my $DS_PFitter_starlike_pvalue = 1;
    my $DS_PFitter_lower_is_starlike = 1;
    my $DS_PFitter_lower_starlike_pvalue = 1;
    my $DS_PFitter_upper_is_starlike = 1;
    my $DS_PFitter_upper_starlike_pvalue = 1;
    my $DS_PFitter_assertion_low = 1.5;
    my $DS_PFitter_assertion_high = 2.5;
    my $DS_PFitter_R = "1.0";
    my ( $DS_starlike_text, $DS_lower_starlike_text, $DS_upper_starlike_text );
    my $paul_calls_one_cluster = 1;
    my $maybe_masked = "";
    if( $run_PFitter ) {
      #if( $mean_diversity == 0 ) {
      #  if( $VERBOSE ) {
      #    print "\nNOT running PoissonFitter because there is no variation whatsoever, and PFitter doesn't handle that case.\n";
      #  }
      #  print "PoissonFitter Determination: Degenerate Phylogeny (all seqs are the same)\n";
      #  print "PoissonFitter Poisson Fit: NA\n";
      #  print "DS Poisson Fit: $DS_PFitter_fitstext (R=$DS_PFitter_R).\n";
      #  print "Average distance to nearest Poisson CDF (2.5%, 97.5% quantiles): $DS_PFitter_distance_mean ($DS_PFitter_distance_ci_low, $DS_PFitter_distance_ci_high)\n";
      #  print "PoissonFitter Poisson time estimate (95\% CI): $PFitter_days_est ($PFitter_days_ci_low, $PFitter_days_ci_high)\n";
      #  #print "\n$PFitter_fitter_stats_raw\n";
      #  print "DS Poisson time estimate (95\% CI): $DS_PFitter_days_est ($DS_PFitter_days_ci_low, $DS_PFitter_days_ci_high)\n";
      #  print "Bayesian Poisson time estimate (95\% CI): $Bayesian_PFitter_days_est ($Bayesian_PFitter_days_ci_low, $Bayesian_PFitter_days_ci_high)\n";
      #} else {
        if( $VERBOSE ) {
          print "\nCalling R to run PoissonFitter..";
        }
        $R_output = `export runPoissonFitter_inputFilename="$fasta_file"; export runPoissonFitter_outputDir="$output_path_dir_for_input_fasta_file"; export runPoissonFitter_runDSPFitter="$run_DSPFitter"; export runPoissonFitter_maskOutNonsynonymousCodons="$fasta_file_mask_out_nonsynonymous_codons_in_PFitter"; R -f runPoissonFitter.R --vanilla --slave`;
        if( $VERBOSE ) {
          print( $R_output );
          print( "done.\n" );
        }
        if( $fasta_file_mask_out_nonsynonymous_codons_in_PFitter ) {
          $maybe_masked = "maskNonsynonymousCodons_";
        }
        my $PFitter_fitter_stats_raw =
          `cat ${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}_${maybe_masked}PoissonFitterDir/LOG_LIKELIHOOD.results.txt`;
        ( $PFitter_lambda, $PFitter_se, $PFitter_nseq, $PFitter_nbases, $PFitter_mean_hd, $PFitter_max_hd, $PFitter_days_est_and_ci, $PFitter_chi_sq_stat, $PFitter_chi_sq_df, $PFitter_chi_sq_p_value  ) =
          (
           $PFitter_fitter_stats_raw =~ /\n[^\t]+\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(\S+)\s*$/ );
        ( $PFitter_days_est, $PFitter_days_ci_low, $PFitter_days_ci_high ) =
          ( $PFitter_days_est_and_ci =~ /(\S+) \((\S+), (\S+)\)/ );

        $is_poisson = ( defined( $PFitter_chi_sq_p_value ) && ( $PFitter_chi_sq_p_value > 0.05 ) ) || 0;
        my $starlike_raw =
          `cat ${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}_${maybe_masked}PoissonFitterDir/CONVOLUTION.results.txt`;
        if( $DEBUG ) {
          ## TODO: REMOVE
          # print "PoissonFitter RAW: $starlike_raw\n";
        }
        my ( $starlike_text ) = ( $starlike_raw =~ /^.+(FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY\s*$/m );
        $is_starlike = "0";
        if( $starlike_text eq "FOLLOWS" ) {
          $is_starlike = "1";
        }

        # DS results
        my $DSPFitter_outfile = "${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}_${maybe_masked}PoissonFitterDir/${fasta_file_short_nosuffix}_${maybe_masked}DSPFitter.out";
        if( !-e $DSPFitter_outfile ) {
          print( "UH OH: Missing DSPFitter output file \"$DSPFitter_outfile\"\n" );
        } else {
          my $DSPFitter_fitter_stats_raw =
            `cat $DSPFitter_outfile`;
          ( $Bayesian_PFitter_lambda_est, $Bayesian_PFitter_lambda_ci_low, $Bayesian_PFitter_lambda_ci_high ) =
            ( $DSPFitter_fitter_stats_raw =~ /Bayesian PFitter Estimated Lambda is (\S+) \(95% CI (\S+) to (\S+)\)/ );
          ( $Bayesian_PFitter_days_est, $Bayesian_PFitter_days_ci_low, $Bayesian_PFitter_days_ci_high ) =
            ( $DSPFitter_fitter_stats_raw =~ /Bayesian PFitter Estimated Days: (\S+) \((\S+), (\S+)\)/ );
          ( $DS_PFitter_lambda_est, $DS_PFitter_lambda_ci_low, $DS_PFitter_lambda_ci_high ) =
            ( $DSPFitter_fitter_stats_raw =~ /DS PFitter Estimated Lambda is (\S+) \(95% CI (\S+) to (\S+)\)/ );
          ( $DS_PFitter_days_est, $DS_PFitter_days_ci_low, $DS_PFitter_days_ci_high ) =
            ( $DSPFitter_fitter_stats_raw =~ /DS PFitter Estimated Days: (\S+) \((\S+), (\S+)\)/ );
          ( $DS_PFitter_distance_mean, $DS_PFitter_distance_ci_low, $DS_PFitter_distance_ci_high ) =
            ( $DSPFitter_fitter_stats_raw =~ /It seems that the CDF of the closest Poisson distribution is roughly (\S+)% away from the pepr-sampled empirical CDFs \(middle 95% (\S+) to (\S+)\)./ );
          ( $DS_PFitter_fitstext ) =
            ( $DSPFitter_fitter_stats_raw =~ /^DSPFitter test that intersequence rate = 2 x seq-consensus rate: (BAD|OK)$/m );
          $DS_PFitter_fits =
            ( ( $DS_PFitter_fitstext =~ /^OK$/ ) ? "1" : "0" );
          ( $DS_PFitter_assertion_low, $DS_PFitter_assertion_high, $DS_PFitter_R ) =
            ( $DSPFitter_fitter_stats_raw =~ /There is .*evidence against the assertion that the Poisson rate between sequences is between (\S+) and (\S+) times the rate of sequences to the consensus \(R \<?= (\S+)\)/ );
    
          ( $DS_starlike_text, $DS_PFitter_starlike_pvalue ) = ( $DSPFitter_fitter_stats_raw =~ /^DSPFitter convolution test: (FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY\s*\(P (?:= )?([^\(]+)\)\s*$/m );
          $DS_PFitter_is_starlike = "0";
          if( $DS_starlike_text eq "FOLLOWS" ) {
            $DS_PFitter_is_starlike = "1";
          }
          ( $DS_lower_starlike_text, $DS_PFitter_lower_starlike_pvalue ) = ( $DSPFitter_fitter_stats_raw =~ /^DSPFitter lower convolution test: (FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY\s*\(P (?:= )?([^\(]+)\)\s*$/m );
          $DS_PFitter_lower_is_starlike = "0";
          if( $DS_lower_starlike_text eq "FOLLOWS" ) {
            $DS_PFitter_lower_is_starlike = "1";
          }
          ( $DS_upper_starlike_text, $DS_PFitter_upper_starlike_pvalue ) = ( $DSPFitter_fitter_stats_raw =~ /^DSPFitter upper convolution test: (FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY\s*\(P (?:= )?([^\(]+)\)\s*$/m );
          $DS_PFitter_upper_is_starlike = "0";
          if( $DS_upper_starlike_text eq "FOLLOWS" ) {
            $DS_PFitter_upper_is_starlike = "1";
          }
        } # End if the DSPFitter output file exists.

        print "PoissonFitter Determination: ";
        if( $is_starlike ) {
          print "Star-Like Phylogeny";
        } else {
          print "Non-Star-Like Phylogeny";
        }
        print "\nPoissonFitter Poisson Fit: ";
        if( $is_poisson ) {
          print "OK\n";
        } else {
          print "BAD (p = $PFitter_chi_sq_p_value)\n";
        }
        print "Average distance to nearest Poisson CDF (2.5%, 97.5% quantiles): $DS_PFitter_distance_mean ($DS_PFitter_distance_ci_low, $DS_PFitter_distance_ci_high)\n";
        print "DS PoissonFitter Determination: ";
        if( $DS_PFitter_is_starlike ) {
          print "Star-Like Phylogeny (p = $DS_PFitter_starlike_pvalue)\n";
        } else {
          print "Non-Star-Like Phylogeny (p = $DS_PFitter_starlike_pvalue)\n";
        }
        print "DS lower PoissonFitter Determination (using 2.5th percentile of epsilon): ";
        if( $DS_PFitter_lower_is_starlike ) {
          print "Star-Like Phylogeny (p = $DS_PFitter_lower_starlike_pvalue)\n";
        } else {
          print "Non-Star-Like Phylogeny (p = $DS_PFitter_lower_starlike_pvalue)\n";
        }
        print "DS upper PoissonFitter Determination (using 97.5th percentile of epsilon): ";
        if( $DS_PFitter_upper_is_starlike ) {
          print "Star-Like Phylogeny (p = $DS_PFitter_upper_starlike_pvalue)\n";
        } else {
          print "Non-Star-Like Phylogeny (p = $DS_PFitter_upper_starlike_pvalue)\n";
        }
        print "DS Poisson Fit: $DS_PFitter_fitstext (R=$DS_PFitter_R).\n";
        print "PoissonFitter Poisson time estimate (95\% CI): $PFitter_days_est ($PFitter_days_ci_low, $PFitter_days_ci_high)\n";
        #print "\n$PFitter_fitter_stats_raw\n";
        print "DS Poisson time estimate (95\% CI): $DS_PFitter_days_est ($DS_PFitter_days_ci_low, $DS_PFitter_days_ci_high)\n";
        print "Bayesian Poisson time estimate (95\% CI): $Bayesian_PFitter_days_est ($Bayesian_PFitter_days_ci_low, $Bayesian_PFitter_days_ci_high)\n";
      #} # End if( $mean_diversity > 0 );
      print OUTPUT_TABLE_FH "\t", $PFitter_lambda;
      print OUTPUT_TABLE_FH "\t", $PFitter_se;
      print OUTPUT_TABLE_FH "\t", $PFitter_nseq;
      print OUTPUT_TABLE_FH "\t", $PFitter_nbases;
      print OUTPUT_TABLE_FH "\t", $PFitter_mean_hd;
      print OUTPUT_TABLE_FH "\t", $PFitter_max_hd;
      print OUTPUT_TABLE_FH "\t", $PFitter_days_est;
      print OUTPUT_TABLE_FH "\t", $PFitter_days_ci_low;
      print OUTPUT_TABLE_FH "\t", $PFitter_days_ci_high;
      print OUTPUT_TABLE_FH "\t", $PFitter_chi_sq_stat;
      print OUTPUT_TABLE_FH "\t", $PFitter_chi_sq_df;
      print OUTPUT_TABLE_FH "\t", $PFitter_chi_sq_p_value;
      print OUTPUT_TABLE_FH "\t", $is_poisson;
      print OUTPUT_TABLE_FH "\t", $is_starlike;
      print OUTPUT_TABLE_FH "\t", $Bayesian_PFitter_lambda_est;
      print OUTPUT_TABLE_FH "\t", $Bayesian_PFitter_lambda_ci_low;
      print OUTPUT_TABLE_FH "\t", $Bayesian_PFitter_lambda_ci_high;
      print OUTPUT_TABLE_FH "\t", $Bayesian_PFitter_days_est;
      print OUTPUT_TABLE_FH "\t", $Bayesian_PFitter_days_ci_low;
      print OUTPUT_TABLE_FH "\t", $Bayesian_PFitter_days_ci_high;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_lambda_est;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_lambda_ci_low;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_lambda_ci_high;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_days_est;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_days_ci_low;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_days_ci_high;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_distance_mean;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_distance_ci_low;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_distance_ci_high;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_assertion_low;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_assertion_high;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_fits;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_R;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_is_starlike;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_starlike_pvalue;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_lower_is_starlike;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_lower_starlike_pvalue;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_upper_is_starlike;
      print OUTPUT_TABLE_FH "\t", $DS_PFitter_upper_starlike_pvalue;

      ## I'll call it more than one cluster if it doesn't conform to the model by one of the two DS versions of the PFitter tests, and its diversity is sufficiently high.
      if( ( ( defined( $mean_diversity ) ? ( ( $mean_diversity > $mean_diversity_threshold ) && ( $mean_diversity > 0 ) ) : 1 ) ) && ( ( $DS_PFitter_fits eq "0" ) || ( ( $DS_PFitter_is_starlike eq "0" ) ) ) ) {
        $paul_calls_one_cluster = 0;
      }
      if( $paul_calls_one_cluster ) {
        print "Number of founders estimated using p-values from DSPfitter is: 1\n";
      } else {
        print "Number of founders estimated using p-values from DSPfitter is: greater than 1\n";
      }

      ## "is.one.founder.alt":
      print OUTPUT_TABLE_FH "\t", $paul_calls_one_cluster;
    } # End if( $run_PFitter )

    my $num_clusters = undef;
    if( $run_InSites_and_PhyML || $run_PFitter ) {
      foreach my $force_one_cluster ( 0, 1 ) {
        my $tmp_extra_flags = $extra_flags;
        if( $force_one_cluster ) {
          $tmp_extra_flags .= "-f ";
        }
        my $clusterInformativeSites_output = `perl clusterInformativeSites.pl $tmp_extra_flags $fasta_file ${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}_informativeSites.txt $output_path_dir_for_input_fasta_file`;
        if( $VERBOSE ) {
          print $clusterInformativeSites_output;
        }
        if( !$force_one_cluster ) {
          # There might be extra text in there.  Look for the telltale "[1]" output
          ( $num_clusters ) =
            ( $clusterInformativeSites_output =~ /\[1\] (\d+?)\s*/ );
          if( $run_InSites_and_PhyML ) {
            if( $morgane_calls_one_cluster ) {
              $in_sites_founders_call = 1;
            } else {
              $in_sites_founders_call = $num_clusters;
              print "Number of clusters found when clustering informative sites: $num_clusters\n";
            }
          }
          if( $run_PFitter ) {
            if( $paul_calls_one_cluster ) {
              $ds_pfitter_founders_call = 1;
            } else {
              $ds_pfitter_founders_call = $num_clusters;
              print "Number of founders estimated by the Informative Sites method using DS PFitter thresholding: $ds_pfitter_founders_call\n";
            }
          }
        } # End if !$force_one_cluster
      } # End foreach $force_one_cluster in ( 0, 1 )
    } # End if $run_InSites_and_PhyML || $run_PFitter

    if( $run_PFitter ) {
      if( $num_clusters == 1 ) {
        ## Avoid NA in the table output.  Multifounder results default to single-founder results.
        print OUTPUT_TABLE_FH "\t", $PFitter_lambda;
        print OUTPUT_TABLE_FH "\t", $PFitter_se;
        print OUTPUT_TABLE_FH "\t", $PFitter_nseq;
        print OUTPUT_TABLE_FH "\t", $PFitter_nbases;
        print OUTPUT_TABLE_FH "\t", $PFitter_mean_hd;
        print OUTPUT_TABLE_FH "\t", $PFitter_max_hd;
        print OUTPUT_TABLE_FH "\t", $PFitter_days_est;
        print OUTPUT_TABLE_FH "\t", $PFitter_days_ci_low;
        print OUTPUT_TABLE_FH "\t", $PFitter_days_ci_high;
        print OUTPUT_TABLE_FH "\t", $PFitter_chi_sq_stat;
        print OUTPUT_TABLE_FH "\t", $PFitter_chi_sq_df;
        print OUTPUT_TABLE_FH "\t", $PFitter_chi_sq_p_value;
        print OUTPUT_TABLE_FH "\t", $is_poisson;
#        print OUTPUT_TABLE_FH "\t", $is_starlike;
        print OUTPUT_TABLE_FH "\t", $Bayesian_PFitter_lambda_est;
        print OUTPUT_TABLE_FH "\t", $Bayesian_PFitter_lambda_ci_low;
        print OUTPUT_TABLE_FH "\t", $Bayesian_PFitter_lambda_ci_high;
        print OUTPUT_TABLE_FH "\t", $Bayesian_PFitter_days_est;
        print OUTPUT_TABLE_FH "\t", $Bayesian_PFitter_days_ci_low;
        print OUTPUT_TABLE_FH "\t", $Bayesian_PFitter_days_ci_high;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_lambda_est;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_lambda_ci_low;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_lambda_ci_high;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_days_est;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_days_ci_low;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_days_ci_high;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_distance_mean;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_distance_ci_low;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_distance_ci_high;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_assertion_low;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_assertion_high;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_fits;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_R;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_is_starlike;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_starlike_pvalue;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_lower_is_starlike;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_lower_starlike_pvalue;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_upper_is_starlike;
        print OUTPUT_TABLE_FH "\t", $DS_PFitter_upper_starlike_pvalue;
      } else {
        die unless( $num_clusters > 1 );
        ## Now run PoissonFitter on the clusters.
        if( $VERBOSE ) {
          print "Calling R to run MultiFounderPoissonFitter..";
        }
        $R_output = `export runMultiFounderPoissonFitter_inputFilenamePrefix="${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}"; export runMultiFounderPoissonFitter_outputDir="$output_path_dir_for_input_fasta_file"; export runMultiFounderPoissonFitter_suffixPattern=""; export runMultiFounderPoissonFitter_runDSPFitter="TRUE"; export runMultiFounderPoissonFitter_maskOutNonsynonymousCodons="$fasta_file_mask_out_nonsynonymous_codons_in_PFitter"; R -f runMultiFounderPoissonFitter.R --vanilla --slave`;
        if( $VERBOSE ) {
          print( $R_output );
          print( "done.\n" );
        }
        my $multifounder_PFitter_fitter_stats_raw =
          `cat ${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}_${maybe_masked}MultiFounderPoissonFitterDir/LOG_LIKELIHOOD.results.txt`;
        #print "\$multifounder_PFitter_fitter_stats_raw: $multifounder_PFitter_fitter_stats_raw\n";
        my ( $multifounder_PFitter_lambda, $multifounder_PFitter_se, $multifounder_PFitter_nseq, $multifounder_PFitter_nbases, $multifounder_PFitter_mean_hd, $multifounder_PFitter_max_hd, $multifounder_PFitter_days_est_and_ci, $multifounder_PFitter_chi_sq_stat, $multifounder_PFitter_chi_sq_df, $multifounder_PFitter_chi_sq_p_value  ) =
          (
           $multifounder_PFitter_fitter_stats_raw =~ /\n[^\t]+\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(\S+)\s*$/ );
        my ( $multifounder_PFitter_days_est, $multifounder_PFitter_days_ci_low, $multifounder_PFitter_days_ci_high ) =
          ( $multifounder_PFitter_days_est_and_ci =~ /(\S+) \((\S+), (\S+)\)/ );
        my $multifounder_is_poisson = ( defined( $multifounder_PFitter_chi_sq_p_value ) && ( $multifounder_PFitter_chi_sq_p_value > 0.05 ) ) || 0;
         ## NOTE THAT the convolution is not set up to handle multi-founder data because the convolution should be done within each founder; so for now we just exclude these results.  TODO: implement multi-founder version of the convolution.
         # my $multifounder_starlike_raw =
         #   `cat ${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}_${maybe_masked}MultiFounderPoissonFitterDir/CONVOLUTION.results.txt`;
         # if( $DEBUG ) {
         #   ## TODO: REMOVE
         #   #print "MULTI-FOUNDER PoissonFitter RAW: $multifounder_starlike_raw\n";
         # }
         # my ( $multifounder_starlike_text ) = ( $multifounder_starlike_raw =~ m/(FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY/ );
         # my $multifounder_is_starlike = ( $multifounder_starlike_text eq "FOLLOWS" );
 
        # DS results
        my $multifounder_DSPFitter_fitter_stats_raw =
          `cat ${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}_${maybe_masked}MultiFounderPoissonFitterDir/${fasta_file_short_nosuffix}_${maybe_masked}DSPFitter.out`;
        my ( $multifounder_Bayesian_PFitter_lambda_est, $multifounder_Bayesian_PFitter_lambda_ci_low, $multifounder_Bayesian_PFitter_lambda_ci_high ) =
           ( $multifounder_DSPFitter_fitter_stats_raw =~ /Bayesian PFitter Estimated Lambda is (\S+) \(95% CI (\S+) to (\S+)\)/ );
        my ( $multifounder_Bayesian_PFitter_days_est, $multifounder_Bayesian_PFitter_days_ci_low, $multifounder_Bayesian_PFitter_days_ci_high ) =
           ( $multifounder_DSPFitter_fitter_stats_raw =~ /Bayesian PFitter Estimated Days: (\S+) \((\S+), (\S+)\)/ );
        my ( $multifounder_DS_PFitter_lambda_est, $multifounder_DS_PFitter_lambda_ci_low, $multifounder_DS_PFitter_lambda_ci_high ) =
           ( $multifounder_DSPFitter_fitter_stats_raw =~ /DS PFitter Estimated Lambda is (\S+) \(95% CI (\S+) to (\S+)\)/ );
        my ( $multifounder_DS_PFitter_days_est, $multifounder_DS_PFitter_days_ci_low, $multifounder_DS_PFitter_days_ci_high ) =
           ( $multifounder_DSPFitter_fitter_stats_raw =~ /DS PFitter Estimated Days: (\S+) \((\S+), (\S+)\)/ );
        my ( $multifounder_DS_PFitter_distance_mean, $multifounder_DS_PFitter_distance_ci_low, $multifounder_DS_PFitter_distance_ci_high ) =
           ( $multifounder_DSPFitter_fitter_stats_raw =~ /It seems that the CDF of the closest Poisson distribution is roughly (\S+)% away from the pepr-sampled empirical CDFs \(middle 95% (\S+) to (\S+)\)./ );
        my ( $multifounder_DS_PFitter_fitstext ) =
           ( $multifounder_DSPFitter_fitter_stats_raw =~ /^DSPFitter test that intersequence rate = 2 x seq-consensus rate: (BAD|OK)$/m );
        my $multifounder_DS_PFitter_fits = "0";
        if( $multifounder_DS_PFitter_fitstext =~ /^OK$/ ) {
          $multifounder_DS_PFitter_fits = "1";
        }
        my ( $multifounder_DS_PFitter_assertion_low, $multifounder_DS_PFitter_assertion_high, $multifounder_DS_PFitter_R ) =
          ( $multifounder_DSPFitter_fitter_stats_raw =~ /There is .*evidence against the assertion that the Poisson rate between sequences is between (\S+) and (\S+) times the rate of sequences to the consensus \(R \<?= (\S+)\)/ );

          my ( $multifounder_DS_starlike_text, $multifounder_DS_PFitter_starlike_pvalue ) = ( $multifounder_DSPFitter_fitter_stats_raw =~ /^DSPFitter convolution test: (FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY\s*\(P (?:= )?([^\(]+)\)\s*$/m );
        my $multifounder_DS_PFitter_is_starlike = "0";
        if( $multifounder_DS_starlike_text eq "FOLLOWS" ) {
          $multifounder_DS_PFitter_is_starlike = "1";
        }
          my ( $multifounder_DS_lower_starlike_text, $multifounder_DS_PFitter_lower_starlike_pvalue ) = ( $multifounder_DSPFitter_fitter_stats_raw =~ /^DSPFitter lower convolution test: (FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY\s*\(P (?:= )?([^\(]+)\)\s*$/m );
        my $multifounder_DS_PFitter_lower_is_starlike = "0";
        if( $multifounder_DS_lower_starlike_text eq "FOLLOWS" ) {
          $multifounder_DS_PFitter_lower_is_starlike = "1";
        }
        my ( $multifounder_DS_upper_starlike_text, $multifounder_DS_PFitter_upper_starlike_pvalue ) = ( $multifounder_DSPFitter_fitter_stats_raw =~ /^DSPFitter upper convolution test: (FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY\s*\(P (?:= )?([^\(]+)\)\s*$/m );
        my $multifounder_DS_PFitter_upper_is_starlike = "0";
        if( $multifounder_DS_upper_starlike_text eq "FOLLOWS" ) {
          $multifounder_DS_PFitter_upper_is_starlike = "1";
        }

         # print "Multi-Founder PoissonFitter Determination: ";
         # if( $multifounder_is_starlike ) {
         #   print "Star-Like Phylogenies within clusters";
         # } else {
         #   print "Non-Star-Like Phylogenies within clusters";
        # }
        print "Multi-Founder Average distance to nearest Poisson CDF (2.5%, 97.5% quantiles): $multifounder_DS_PFitter_distance_mean ($multifounder_DS_PFitter_distance_ci_low, $multifounder_DS_PFitter_distance_ci_high)\n";
        print "Multi-Founder DS PoissonFitter Determination: ";
        if( $multifounder_DS_PFitter_is_starlike ) {
          print "Star-Like Phylogeny (p = $multifounder_DS_PFitter_starlike_pvalue)\n";
        } else {
          print "Non-Star-Like Phylogeny (p = $multifounder_DS_PFitter_starlike_pvalue)\n";
        }
        print "Multi-Founder DS lower PoissonFitter Determination (using 2.5th percentile of epsilon): ";
        if( $multifounder_DS_PFitter_lower_is_starlike ) {
          print "Star-Like Phylogeny (p = $multifounder_DS_PFitter_lower_starlike_pvalue)\n";
        } else {
          print "Non-Star-Like Phylogeny (p = $multifounder_DS_PFitter_lower_starlike_pvalue)\n";
        }
        print "Multi-Founder DS upper PoissonFitter Determination (using 97.5th percentile of epsilon): ";
        if( $multifounder_DS_PFitter_upper_is_starlike ) {
          print "Star-Like Phylogeny (p = $multifounder_DS_PFitter_upper_starlike_pvalue)\n";
        } else {
          print "Non-Star-Like Phylogeny (p = $multifounder_DS_PFitter_upper_starlike_pvalue)\n";
        }
        print "Multi-Founder DS Poisson Fit: $multifounder_DS_PFitter_fitstext (R=$multifounder_DS_PFitter_R).\n";
        print "Multi-Founder PoissonFitter Poisson time estimate (95\% CI): $multifounder_PFitter_days_est ($multifounder_PFitter_days_ci_low, $multifounder_PFitter_days_ci_high)\n";
        print "Multi-Founder DS Poisson time estimate (95\% CI): $multifounder_DS_PFitter_days_est ($multifounder_DS_PFitter_days_ci_low, $multifounder_DS_PFitter_days_ci_high)\n";
        print "Multi-Founder Bayesian Poisson time estimate (95\% CI): $multifounder_Bayesian_PFitter_days_est ($multifounder_Bayesian_PFitter_days_ci_low, $multifounder_Bayesian_PFitter_days_ci_high)\n";
        #print "\n$multifounder_PFitter_fitter_stats_raw\n";
        
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_lambda;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_se;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_nseq;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_nbases;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_mean_hd;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_max_hd;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_days_est;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_days_ci_low;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_days_ci_high;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_chi_sq_stat;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_chi_sq_df;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_chi_sq_p_value;
        print OUTPUT_TABLE_FH "\t", $multifounder_is_poisson;
 #       print OUTPUT_TABLE_FH "\t", $multifounder_is_starlike;
        print OUTPUT_TABLE_FH "\t", $multifounder_Bayesian_PFitter_lambda_est;
        print OUTPUT_TABLE_FH "\t", $multifounder_Bayesian_PFitter_lambda_ci_low;
        print OUTPUT_TABLE_FH "\t", $multifounder_Bayesian_PFitter_lambda_ci_high;
        print OUTPUT_TABLE_FH "\t", $multifounder_Bayesian_PFitter_days_est;
        print OUTPUT_TABLE_FH "\t", $multifounder_Bayesian_PFitter_days_ci_low;
        print OUTPUT_TABLE_FH "\t", $multifounder_Bayesian_PFitter_days_ci_high;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_lambda_est;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_lambda_ci_low;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_lambda_ci_high;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_days_est;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_days_ci_low;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_days_ci_high;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_distance_mean;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_distance_ci_low;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_distance_ci_high;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_assertion_low;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_assertion_high;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_fits;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_R;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_is_starlike;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_starlike_pvalue;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_lower_is_starlike;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_lower_starlike_pvalue;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_upper_is_starlike;
        print OUTPUT_TABLE_FH "\t", $multifounder_DS_PFitter_upper_starlike_pvalue;
      } # End if $num_clusters > 1
    } # End if $run_PFitter

    ## Now try it the more profillic way.  This is a hybrid approach that makes profiles only of the informative sites.
    if( $run_profillic ) {
        my $alignment_profiles_output_file = "${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}_profileToAlignmentProfile.alignmentprofs";
        if( $VERBOSE ) {
          print "Running Profillic..\n";
        }
        my $alignment_profiles_output_files_list_file = "${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}_profillic_AlignmentProfilesList.txt";
        my $runProfillic_output = `perl runProfillic.pl $extra_flags $fasta_file $alignment_profiles_output_files_list_file $output_path_dir_for_input_fasta_file`;
        if( $VERBOSE ) {
          print( $runProfillic_output );
        }

        if( $VERBOSE ) {
          print "Clustering..\n";
        }
        $R_output = `perl clusterProfillicAlignmentProfiles.pl $extra_flags $fasta_file $alignment_profiles_output_files_list_file $output_path_dir_for_input_fasta_file`;
        if( $VERBOSE ) {
          print( $R_output );
        }
        my ( $num_profillic_clusters ) = ( $R_output =~ /^\[1\]\s*(\d+)\s*$/m );
        # Print out the number of clusters
        print "Number of clusters found using profillic: $num_profillic_clusters\n";
        print OUTPUT_TABLE_FH "\t", $num_profillic_clusters;
        if( $paul_calls_one_cluster ) {
          print "Number of founders estimated by the Informative Sites Profillic method: 1\n";
          print OUTPUT_TABLE_FH "\t", 1;
        } else {
          print "Number of founders estimated by the Informative Sites Profillic method: $num_profillic_clusters\n";
          print OUTPUT_TABLE_FH "\t", $num_profillic_clusters;
        }
    } # End if $run_profillic

    if( $run_InSites_and_PhyML || $run_PFitter ) {
      print OUTPUT_TABLE_FH "\t", $num_clusters;
    }
    if( $run_InSites_and_PhyML ) {
      print OUTPUT_TABLE_FH "\t", $in_sites_founders_call;
    }
    if( $run_PFitter ) {
      print OUTPUT_TABLE_FH "\t", $ds_pfitter_founders_call;
    }
    print OUTPUT_TABLE_FH "\n";
    
    # If this is a half-genome dataset, should we call
    # runMultiFounderPoissonFitter to put together the two datasets
    # for better Poisson estimation?
    if( $run_PFitter ) {
      # NOTE that this strips off any modifiers like remove/fix HypermutatedSequences and removeRecombinedSequences: (this is dealt with in runMultiFounderPoissonFitter.R)
      my ( $fasta_file_very_short, $fasta_file_region ) =
         ( $fasta_file_short =~ /^(.+)(_LH|_NFLG|_RH|_env)/ );
      my @region_fasta_files = ();
      # Is this the last one in the order of the input list?
      my $fasta_file_is_last_alternate_region_in_list = 0;
      if( defined( $fasta_file_very_short ) ) {
        ## TODO: REMOVE
        #print "\$fasta_file_very_short is $fasta_file_very_short\n";
        ## TODO: REMOVE
        #print "\$input_fasta_file is $input_fasta_file\n";
        ## TODO: REMOVE
        #print "\@all_input_fasta_files: ( ", join( ", ", @all_input_fasta_files ), " )\n";
        #print "GREP: /$fasta_file_very_short(_LH|_NFLG|_RH|_env)/\n";
        @region_fasta_files = grep( /$fasta_file_very_short(_LH|_NFLG|_RH|_env)/, @all_input_fasta_files );
        @region_fasta_files = map { if( $_ !~ /\// ) { "./$_" } else { $_ } } @region_fasta_files;
      }
      ## TODO: REMOVE
      #print "\@region_fasta_files: ( ", join( ", ", @region_fasta_files ), " )\n";
      if( scalar( @region_fasta_files ) >= 2 ) {
        $fasta_file_is_last_alternate_region_in_list =
          ( $input_fasta_file eq $region_fasta_files[ $#region_fasta_files ] );
      }
      #print "\$fasta_file: $fasta_file\n";
      #print "\$fasta_file_is_last_alternate_region_in_list: $fasta_file_is_last_alternate_region_in_list\n";
      if( $fasta_file_is_last_alternate_region_in_list ) {
          # OK, do it.
          ## Now run PoissonFitter on the regions.
          if( $VERBOSE ) {
            print "Calling R to run MultiRegionPoissonFitter..";
          }
          #print "\$fasta_file_very_short: $fasta_file_very_short\n"; 
          if( !exists ( $final_input_fasta_file_short_names_by_original_short_name_stripped{ $fasta_file_very_short } ) ) {
            die( "\$final_input_fasta_file_short_names_by_original_short_name_stripped{ $fasta_file_very_short } doesn't exist!" );
          }
          if( $VERBOSE ) {
            print "[Combining these sequences: ", join( ",", @{ $final_input_fasta_file_short_names_by_original_short_name_stripped{ $fasta_file_very_short } } ), "]..";
          }
          
          ## BUT WAIT.  If there are clusters, gather the cluster files instead.
          my @files_for_regions = @{ $final_input_fasta_file_short_names_by_original_short_name_stripped{ $fasta_file_very_short } };
          ## TODO: ONLY DO THIS IF THE (respective) INFORMATIVE_SITES CALL WAS > 1 [for consistent behavior when using -s flag]
          ## For each one, see if there are cluster files.
          opendir OUTPUT_DIR, $output_path_dir_for_input_fasta_file;
          my @all_files = readdir( OUTPUT_DIR );
          closedir OUTPUT_DIR;
          ## TODO: REMOVE
          #print "\nALL: ", join( ", ", @all_files ), "\n";
          ## These are anonymous functions because in Perl the rules for scope differ for anonymous functions.  Access to "my" variables defined (lexically) outside the scope of anonymous functions from within them are as expected (changes to the variable are seen by the inner function); access to "my" variables defined in outer functions and accessed within named functions (which in c parlance are static) access their values at first call to the outer function.  See http://stackoverflow.com/questions/4048248/variable-foo-will-not-stay-shared-warning-error-in-perl-while-calling-subrout
              my $foo = sub {
                my $fileprefix = shift;
                my $file_name = shift;
                my $inwithcluster = "^" . $fileprefix . "_cluster\\d+\\.fasta\$";
                #print "\$inwithcluster is $inwithcluster\n";
                my $isamatch = ( $file_name =~ m/$inwithcluster/ );
                if( $isamatch ) {
                  #print "MATCH: $file_name\n";
                } else {
                  #print "MISMATCH: $file_name\n";
                }
                return( $isamatch );
              }; # End sub $foo
              my $baz = sub {
                my $fn = shift;
                my ( $rv ) = ( $fn =~ /^(.+)\.fasta$/ );
                #print "stripped: $rv\n";
                return( $rv );
              }; # End sub $baz
          my $bar = sub  {
            my $fileprefix = shift;
            #print $fileprefix, "\n";
            if( ( -e "${output_path_dir_for_input_fasta_file}/${fileprefix}_cluster0.fasta" ) || ( -e "${output_path_dir_for_input_fasta_file}/${fileprefix}_cluster1.fasta" ) ) {
              #print "YES\n";
              my @cluster_files = grep { &$foo( $fileprefix, $_ ) } @all_files;
              # Strip off the .fasta suffix.
              #print "ok\n";
              my @stripped_cluster_files = map { &$baz($_) } @cluster_files;
              #print "still\n";
              #print( "\@stripped_cluster_files: ( " . join( ", ", @stripped_cluster_files ) . " )\n" );
              return( @stripped_cluster_files );
            } else {
              #print "NO\n";
              my @lst;
              push @lst, $fileprefix;
              return( @lst );
            }
          }; # End sub $bar
          my @new_files_for_regions;
          map { my @lst = &$bar($_); push @new_files_for_regions, @lst; $_ } @files_for_regions;
          ## TODO: REMOVE
          #print "\nUSING: ". join( ", ", @new_files_for_regions ). "\n";
          my @file_suffixes = map { ( $_ ) = ( $_ =~ /^$fasta_file_very_short(.+)$/ ); $_ } @new_files_for_regions;
          my $suffix_pattern = join( "\.fasta|", @file_suffixes ) . "\.fasta";
          # print "\$fasta_file_very_short: $fasta_file_very_short\n";
          #print "\$suffix_pattern: $suffix_pattern\n";
          $R_output = `export runMultiFounderPoissonFitter_inputFilenamePrefix='${output_path_dir_for_input_fasta_file}/${fasta_file_very_short}'; export runMultiFounderPoissonFitter_suffixPattern='$suffix_pattern'; export runMultiFounderPoissonFitter_outputDir='$output_path_dir_for_input_fasta_file'; export runMultiFounderPoissonFitter_runDSPFitter='TRUE'; export runMultiFounderPoissonFitter_maskOutNonsynonymousCodons="$fasta_file_mask_out_nonsynonymous_codons_in_PFitter"; R -f runMultiFounderPoissonFitter.R --vanilla --slave`;
          if( $VERBOSE ) {
            print $R_output;
            print( "done.\n" );
          }
          my $multi_region_PFitter_fitter_stats_raw =
            `cat ${output_path_dir_for_input_fasta_file}/${fasta_file_very_short}_${maybe_masked}MultiRegionPoissonFitterDir/LOG_LIKELIHOOD.results.txt`;
          my ( $multi_region_PFitter_lambda, $multi_region_PFitter_se, $multi_region_PFitter_nseq, $multi_region_PFitter_nbases, $multi_region_PFitter_mean_hd, $multi_region_PFitter_max_hd, $multi_region_PFitter_days_est_and_ci, $multi_region_PFitter_chi_sq_stat, $multi_region_PFitter_chi_sq_df, $multi_region_PFitter_chi_sq_p_value  ) =
            (
             $multi_region_PFitter_fitter_stats_raw =~ /\n[^\t]+\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(\S+)\s*$/ );
         my ( $multi_region_PFitter_days_est, $multi_region_PFitter_days_ci_low, $multi_region_PFitter_days_ci_high ) =
           ( $multi_region_PFitter_days_est_and_ci =~ /(\S+) \((\S+), (\S+)\)/ );
          my $multi_region_is_poisson = ( defined( $multi_region_PFitter_chi_sq_p_value ) && ( $multi_region_PFitter_chi_sq_p_value > 0.05 ) ) || 0;
          ## NOTE THAT the convolution is not set up to handle multi-region data because the convolution should be done within each region; so for now we just exclude these results.  TODO: implement multi-region version of the convolution.
          # my $multi_region_starlike_raw =
          #   `cat ${output_path_dir_for_input_fasta_file}/${fasta_file_short_nosuffix}_${maybe_masked}MultiRegionPoissonFitterDir/CONVOLUTION.results.txt`;
          # if( $DEBUG ) {
          #   ## TODO: REMOVE
          #   #print "MULTI-REGION PoissonFitter RAW: $multi_region_starlike_raw\n";
          # }
          # my ( $multi_region_starlike_text ) = ( $multi_region_starlike_raw =~ m/(FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY/ );
          # my $multi_region_is_starlike = ( $multi_region_starlike_text eq "FOLLOWS" );

         # DS results
         my $multi_region_DSPFitter_fitter_stats_raw =
           `cat ${output_path_dir_for_input_fasta_file}/${fasta_file_very_short}_${maybe_masked}MultiRegionPoissonFitterDir/${fasta_file_very_short}_${maybe_masked}DSPFitter.out`;
         my ( $multi_region_Bayesian_PFitter_lambda_est, $multi_region_Bayesian_PFitter_lambda_ci_low, $multi_region_Bayesian_PFitter_lambda_ci_high ) =
            ( $multi_region_DSPFitter_fitter_stats_raw =~ /Bayesian PFitter Estimated Lambda is (\S+) \(95% CI (\S+) to (\S+)\)/ );
         my ( $multi_region_Bayesian_PFitter_days_est, $multi_region_Bayesian_PFitter_days_ci_low, $multi_region_Bayesian_PFitter_days_ci_high ) =
            ( $multi_region_DSPFitter_fitter_stats_raw =~ /Bayesian PFitter Estimated Days: (\S+) \((\S+), (\S+)\)/ );
         my ( $multi_region_DS_PFitter_lambda_est, $multi_region_DS_PFitter_lambda_ci_low, $multi_region_DS_PFitter_lambda_ci_high ) =
            ( $multi_region_DSPFitter_fitter_stats_raw =~ /DS PFitter Estimated Lambda is (\S+) \(95% CI (\S+) to (\S+)\)/ );
         my ( $multi_region_DS_PFitter_days_est, $multi_region_DS_PFitter_days_ci_low, $multi_region_DS_PFitter_days_ci_high ) =
            ( $multi_region_DSPFitter_fitter_stats_raw =~ /DS PFitter Estimated Days: (\S+) \((\S+), (\S+)\)/ );
         my ( $multi_region_DS_PFitter_distance_mean, $multi_region_DS_PFitter_distance_ci_low, $multi_region_DS_PFitter_distance_ci_high ) =
            ( $multi_region_DSPFitter_fitter_stats_raw =~ /It seems that the CDF of the closest Poisson distribution is roughly (\S+)% away from the pepr-sampled empirical CDFs \(middle 95% (\S+) to (\S+)\)./ );
         my ( $multi_region_DS_PFitter_fitstext ) =
            ( $multi_region_DSPFitter_fitter_stats_raw =~ /^DSPFitter test that intersequence rate = 2 x seq-consensus rate: (BAD|OK)$/m );
         my $multi_region_DS_PFitter_fits =
            ( ( $multi_region_DS_PFitter_fitstext =~ /^OK$/ ) ? "1" : "0" );
         my ( $multi_region_DS_PFitter_assertion_low, $multi_region_DS_PFitter_assertion_high, $multi_region_DS_PFitter_R ) =
           ( $multi_region_DSPFitter_fitter_stats_raw =~ /There is .*evidence against the assertion that the Poisson rate between sequences is between (\S+) and (\S+) times the rate of sequences to the consensus \(R = (\S+)\)/ );

          print "\nInput fasta file: ${fasta_file_very_short}${fasta_file_suffix}\n";
          
          my ( $multi_region_DS_starlike_text, $multi_region_DS_PFitter_starlike_pvalue ) = ( $multi_region_DSPFitter_fitter_stats_raw =~ /^DSPFitter convolution test: (FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY\s*\(P (?:= )?([^\(]+)\)\s*$/m );
        my $multi_region_DS_PFitter_is_starlike = "0";
        if( $multi_region_DS_starlike_text eq "FOLLOWS" ) {
          $multi_region_DS_PFitter_is_starlike = "1";
        }
          my ( $multi_region_DS_lower_starlike_text, $multi_region_DS_PFitter_lower_starlike_pvalue ) = ( $multi_region_DSPFitter_fitter_stats_raw =~ /^DSPFitter lower convolution test: (FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY\s*\(P (?:= )?([^\(]+)\)\s*$/m );
        my $multi_region_DS_PFitter_lower_is_starlike = "0";
        if( $multi_region_DS_lower_starlike_text eq "FOLLOWS" ) {
          $multi_region_DS_PFitter_lower_is_starlike = "1";
        }
          my ( $multi_region_DS_upper_starlike_text, $multi_region_DS_PFitter_upper_starlike_pvalue ) = ( $multi_region_DSPFitter_fitter_stats_raw =~ /^DSPFitter upper convolution test: (FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY\s*\(P (?:= )?([^\(]+)\)\s*$/m );
        my $multi_region_DS_PFitter_upper_is_starlike = "0";
        if( $multi_region_DS_upper_starlike_text eq "FOLLOWS" ) {
          $multi_region_DS_PFitter_upper_is_starlike = "1";
        }
          
#         print "Multi-Region PoissonFitter Determination: ";
#         if( $is_starlike ) {
#           print "Star-Like Phylogeny\n";
#         } else {
#           print "Non-Star-Like Phylogeny\n";
#         }
        print "Multi-Region PoissonFitter Poisson Fit: ";
        if( $is_poisson ) {
          print "OK\n";
        } else {
          print "BAD (p = $PFitter_chi_sq_p_value)\n";
        }
        print "Multi-Region average distance to nearest Poisson CDF (2.5%, 97.5% quantiles): $multi_region_DS_PFitter_distance_mean ($multi_region_DS_PFitter_distance_ci_low, $multi_region_DS_PFitter_distance_ci_high)\n";
        print "Multi-Region DS PoissonFitter Determination: ";
        if( $multi_region_DS_PFitter_is_starlike ) {
          print "Star-Like Phylogeny (p = $multi_region_DS_PFitter_starlike_pvalue)\n";
        } else {
          print "Non-Star-Like Phylogeny (p = $multi_region_DS_PFitter_starlike_pvalue)\n";
        }
        print "Multi-Region DS lower PoissonFitter Determination (using 2.5th percentile of epsilon): ";
        if( $multi_region_DS_PFitter_lower_is_starlike ) {
          print "Star-Like Phylogeny (p = $multi_region_DS_PFitter_lower_starlike_pvalue)\n";
        } else {
          print "Non-Star-Like Phylogeny (p = $multi_region_DS_PFitter_lower_starlike_pvalue)\n";
        }
        print "Multi-Region DS upper PoissonFitter Determination (using 97.5th percentile of epsilon): ";
        if( $multi_region_DS_PFitter_upper_is_starlike ) {
          print "Star-Like Phylogeny (p = $multi_region_DS_PFitter_upper_starlike_pvalue)\n";
        } else {
          print "Non-Star-Like Phylogeny (p = $multi_region_DS_PFitter_upper_starlike_pvalue)\n";
        }

          print "Multi-Region DS Poisson Fit: $multi_region_DS_PFitter_fitstext (R=$multi_region_DS_PFitter_R).\n";
          print "Multi-Region PFitter Poisson time estimate (95\% CI): $multi_region_PFitter_days_est ($multi_region_PFitter_days_ci_low, $multi_region_PFitter_days_ci_high)\n";
          print "Multi-Region DS Poisson time estimate (95\% CI): $multi_region_DS_PFitter_days_est ($multi_region_DS_PFitter_days_ci_low, $multi_region_DS_PFitter_days_ci_high)\n";
          print "Multi-Region Bayesian Poisson time estimate (95\% CI): $multi_region_Bayesian_PFitter_days_est ($multi_region_Bayesian_PFitter_days_ci_low, $multi_region_Bayesian_PFitter_days_ci_high)\n";
          #print "\n$multi_region_PFitter_fitter_stats_raw\n";
          
          ## I'll call it more than one cluster if it doesn't conform to the model by one of the two DS versions of the PFitter tests.
          my $multi_region_paul_calls_one_cluster = 1;
          if( ( $multi_region_DS_PFitter_fits eq "0" ) || ( $multi_region_DS_PFitter_is_starlike eq "0" ) ) {
            $multi_region_paul_calls_one_cluster = 0;
          }
          if( $multi_region_paul_calls_one_cluster ) {
            print "Multi-Region Number of founders estimated using p-values from DSPfitter is: 1\n";
          } else {
            print "Multi-Region Number of founders estimated using p-values from DSPfitter is: greater than 1\n";
          }

          print OUTPUT_TABLE_FH "${fasta_file_very_short}${fasta_file_suffix}";
          if( $run_Hypermut ) {
            print OUTPUT_TABLE_FH "\t", "NA"; # fixed-hypermut or removed-hypermut
          }
          if( $run_RAP ) {
            print OUTPUT_TABLE_FH "\t", "NA"; # removed-recomb
          }
          print OUTPUT_TABLE_FH "\t", "NA"; # file
          if( $run_InSites_and_PhyML ) {
            print OUTPUT_TABLE_FH "\t", "NA"; # num.seqs
            print OUTPUT_TABLE_FH "\t", "NA"; # num.phyml.seqs
            print OUTPUT_TABLE_FH "\t", "NA"; # diversity
            print OUTPUT_TABLE_FH "\t", "NA"; # inf.to.priv.ratio
            print OUTPUT_TABLE_FH "\t", "NA"; # exceeds.diversity.threshold
            print OUTPUT_TABLE_FH "\t", "NA"; # exceeds.ratio.threshold
            print OUTPUT_TABLE_FH "\t", "NA"; # ratio.threshold
            print OUTPUT_TABLE_FH "\t", "NA"; # is.one.founder
          }
          if( $run_entropy ) {
            print OUTPUT_TABLE_FH "\t", "NA"; # mean.entropy
            print OUTPUT_TABLE_FH "\t", "NA"; # sd.entropy
          }
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_lambda;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_se;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_nseq;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_nbases;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_mean_hd;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_max_hd;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_days_est;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_days_ci_low;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_days_ci_high;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_chi_sq_stat;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_chi_sq_df;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_chi_sq_p_value;
          print OUTPUT_TABLE_FH "\t", $multi_region_is_poisson;
          print OUTPUT_TABLE_FH "\t", "NA"; # is.starlike
          print OUTPUT_TABLE_FH "\t", $multi_region_Bayesian_PFitter_lambda_est;
          print OUTPUT_TABLE_FH "\t", $multi_region_Bayesian_PFitter_lambda_ci_low;
          print OUTPUT_TABLE_FH "\t", $multi_region_Bayesian_PFitter_lambda_ci_high;
          print OUTPUT_TABLE_FH "\t", $multi_region_Bayesian_PFitter_days_est;
          print OUTPUT_TABLE_FH "\t", $multi_region_Bayesian_PFitter_days_ci_low;
          print OUTPUT_TABLE_FH "\t", $multi_region_Bayesian_PFitter_days_ci_high;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_lambda_est;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_lambda_ci_low;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_lambda_ci_high;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_days_est;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_days_ci_low;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_days_ci_high;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_distance_mean;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_distance_ci_low;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_distance_ci_high;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_assertion_low;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_assertion_high;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_fits;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_R;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_is_starlike;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_starlike_pvalue;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_lower_is_starlike;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_lower_starlike_pvalue;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_upper_is_starlike;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_upper_starlike_pvalue;
          print OUTPUT_TABLE_FH "\t", $multi_region_paul_calls_one_cluster; # "is.one.founder.alt"
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_lambda;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_se;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_nseq;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_nbases;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_mean_hd;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_max_hd;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_days_est;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_days_ci_low;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_days_ci_high;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_chi_sq_stat;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_chi_sq_df;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_chi_sq_p_value;
          print OUTPUT_TABLE_FH "\t", $multi_region_is_poisson;
          print OUTPUT_TABLE_FH "\t", $multi_region_Bayesian_PFitter_lambda_est;
          print OUTPUT_TABLE_FH "\t", $multi_region_Bayesian_PFitter_lambda_ci_low;
          print OUTPUT_TABLE_FH "\t", $multi_region_Bayesian_PFitter_lambda_ci_high;
          print OUTPUT_TABLE_FH "\t", $multi_region_Bayesian_PFitter_days_est;
          print OUTPUT_TABLE_FH "\t", $multi_region_Bayesian_PFitter_days_ci_low;
          print OUTPUT_TABLE_FH "\t", $multi_region_Bayesian_PFitter_days_ci_high;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_lambda_est;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_lambda_ci_low;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_lambda_ci_high;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_days_est;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_days_ci_low;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_days_ci_high;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_distance_mean;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_distance_ci_low;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_distance_ci_high;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_assertion_low;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_assertion_high;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_fits;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_R;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_is_starlike;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_starlike_pvalue;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_lower_is_starlike;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_lower_starlike_pvalue;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_upper_is_starlike;
          print OUTPUT_TABLE_FH "\t", $multi_region_DS_PFitter_upper_starlike_pvalue;
          if( $run_profillic ) {
            print OUTPUT_TABLE_FH "\t", "NA"; # "profillic.clusters"
            print OUTPUT_TABLE_FH "\t", "NA"; # "profillic.founder.call"
          } # End if $run_profillic
          if( $run_InSites_and_PhyML ) {
            print OUTPUT_TABLE_FH "\t", "NA"; # founder.call
          }
          if( $run_PFitter ) {
            print OUTPUT_TABLE_FH "\t", "NA"; # founder.call.alt
          }
          print OUTPUT_TABLE_FH "\n";

      } # End if this is the last region, consider combining for MultiRegionPoissonFitter. .. else ..

    } # End if $run_PFitter

    if( $recurse_on_clusters && ( $num_clusters > 1 ) ) {
      # Add them.
      opendir OUTPUT_DIR, $output_path_dir_for_input_fasta_file;
      my @cluster_files = grep { /^${fasta_file_short_nosuffix}_cluster\d+\.fasta$/ } readdir( OUTPUT_DIR );
      closedir OUTPUT_DIR;
      # print "CLUSTER FILES: ", join( ", ", @cluster_files ), "\n";
      my @new_input_fasta_files = map { $output_path_dir_for_input_fasta_file . "/" . $_ } @cluster_files;
      # print "NEW INPUT FILES: ", join( ", ", @new_input_fasta_files ), "\n";
      unshift @fasta_files_remaining_to_process, @new_input_fasta_files;
    } # End if( $recurse_on_clusters )
    
  } # End foreach $fasta_file @fasta_files_remaining_to_process
  if( $VERBOSE ) { print ".done.\n"; }

  if( $VERBOSE && $output_table_file ) { print "Closing file \"$output_table_file\".."; }
  close OUTPUT_TABLE_FH;
  if( $VERBOSE && $output_table_file ) { print ".done.\n"; }

  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  return 0;
} # identify_founders(..)

identify_founders( @ARGV );

1;


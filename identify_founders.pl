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
##      R packages you'll need:
##      ade4
##      seqinr (from bioconductor)
##      dynamicTreeCut
##
##      Try: rm -r rv217_1W_gold_standard-hiv-founder-id_-fs_resultDir/; mkdir rv217_1W_gold_standard-hiv-founder-id_-fs_resultDir/; perl -w ./identify_founders.pl -sf -O rv217_1W_gold_standard-hiv-founder-id_-fs_resultDir/ ~/src/from-git/projects/tholzman/MorgansFounderIDMethod/rv217_1W_gold_standard.list > rv217_1W_gold_standard-hiv-founder-id_-fs_resultDir/identify-founders.out
##      This next one skips RAP because it seems to be exceptionally slow with these large numbers of sequences.  Unsure how to proceed - maybe iterately evaluate subsets? Or use a different program.
##      Or: mkdir caprisa002_1W_gold_standard-hiv-founder-id_resultDir/; perl ./identify_founders.pl -V -R -P -O caprisa002_1W_gold_standard-hiv-founder-id_resultDir/ caprisa002_1W_gold_standard.list  > caprisa002_1W_gold_standard-hiv-founder-id_resultDir/identify-founders.out
##      Or: mkdir CAPRISA002_ft_seqs-hiv-founder-id_resultDir/; perl ./identify_founders.pl -O CAPRISA002_ft_seqs-hiv-founder-id_resultDir/ ~/src/from-git/projects/tholzman/MorgansFounderIDMethod/CAPRISA002_ft_seqs.txt  > CAPRISA002_ft_seqs-hiv-founder-id_resultDir/identify-founders.out
##      Or: mkdir Abrahams-2009aa-hiv-founder-id_resultDir/; perl ./identify_founders.pl -V -O Abrahams-2009aa-hiv-founder-id_resultDir/ Abrahams-2009aa/preparedFor_hiv-identify-founders.list > Abrahams-2009aa-hiv-founder-id_resultDir/identify-founders.out 
##      Whynot: mkdir new-Abrahams-2009aa-hiv-founder-id_resultDir/; perl ./identify_founders.pl -V -O new-Abrahams-2009aa-hiv-founder-id_resultDir/ Abrahams-2009aa/preparedFor_hiv-identify-founders.list > new-Abrahams-2009aa-hiv-founder-id_resultDir/new-identify-founders.out 
###******************************************************************************

use Getopt::Std; # for getopts
use File::Path qw( make_path );
use Path::Tiny;
require Sort::Fields; # for ??

use strict;
use vars qw( $opt_D $opt_V $opt_o $opt_O $opt_P $opt_R $opt_F $opt_H $opt_f $opt_I $opt_s );
use vars qw( $VERBOSE $DEBUG );

## NOTE: If there is only one sequence in a group, we don't bother writing out the file.
sub splitFastaFileOnHeaderPatterns {
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
    print "\tidentify_founders [-DV] [-PRFHfIs] [-(o|O) <output_dir>] <input_fasta_file1 or file_listing_input_files> [<input_fasta_file2> ...] \n";
    exit;
  }

  # This means -D, -o, -O, -V, ... are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_o is an optional directory to put the output; default depends on input filename.
  # opt_O is just like opt_o.  Same thing.
  # opt_V means be verbose.
  # opt_P means skip (do not run) profillic
  # opt_R means skip (do not run) RAP (alignment-internal recombination detection program)
  # opt_F means skip (do not run) PFitter
  # opt_H means skip (do not run) HYPERMUT 2.0 hypermutation detection
  # opt_f means fix hypermutated sequences, instead of removing them.
  # opt_I means run inSites online (instead of offline)
  # opt_s means be slow, ie run everything even when the diversity and insites thresholds are not exceeded.
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_o, $opt_O, $opt_P, $opt_R, $opt_F, $opt_H, $opt_f, $opt_I, $opt_s ) = ();
  if( not getopts('DVo:O:PRFHfIs') ) {
    identify_founders_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;

  my $run_profillic = !$opt_P;
  my $run_RAP = !$opt_R;
  my $run_PFitter = !$opt_F;
  my $run_Hypermut = !$opt_H;
  my $fix_hypermutated_sequences = $opt_f || 0;
  my $runInSites_online = $opt_I || 0;
  my $be_slow = $opt_s || 0;

  ## TODO: make into an arg
  my $run_DSPFitter = 1;
  
  ## TODO: DEHACKIFY MAGIC #s
  my $mean_diversity_threshold = 0.001;
  my $in_sites_ratio_threshold = 0.85; # For Abrahams and RV217
  #my $in_sites_ratio_threshold = 0.33; # For caprisa002
  my $RAP_pValueThreshold = 0.0007; # Appears to be the suggestion from the output file "(summaryTable)"'s column header, which reads "Pvalues<0.0007".
  my $hypermut2_pValueThreshold = 0.1; # Matches Abrahams 2009

  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }
  
  my @input_fasta_files = @ARGV;
  unless( scalar @input_fasta_files ) { identify_founders_usage(); }

  # Special case: if there is only one file, it might be a textfile
  # containing names of files.
  if( scalar( @input_fasta_files ) == 1 ) {
    #print $input_fasta_files[ 0 ], "\n";
    unless( ( $input_fasta_files[ 0 ] ) =~ /\.fast?a?$/ ) {
      # Try opening it, see if it's a list of files.
      if( $VERBOSE ) {
        print "Reading file names from file \"", $input_fasta_files[ 0 ], "\"..";
      }
      my $input_fasta_file_contents = path( $input_fasta_files[ 0 ] )->slurp();
      if( $VERBOSE ) {
        print ".done\n";
      }
      if( $DEBUG ) {
        print $input_fasta_file_contents;
      }
      @input_fasta_files = split( "\n", $input_fasta_file_contents );
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

  my $default_force_one_cluster = 1;
  if( $be_slow ) {
    $default_force_one_cluster = 0;
  }

  my $output_table_file = $output_path_dir . '/' . "identify_founders.tab";

  if( $VERBOSE ) { print "Opening file \"$output_table_file\" for writing.."; }
  unless( open OUTPUT_TABLE_FH, ">$output_table_file" ) {
    warn "Unable to open output file \"$output_table_file\": $!";
    return 1;
  }
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
   "file", "num.seqs", "diversity", "inf.to.priv.ratio",
   "exceeds.diversity.threshold", "exceeds.ratio.threshold",
   "is.one.founder"
   );

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
     "DS.PFitter.R"                     
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
     "multifounder.DS.PFitter.R"                     
    );
  } # End if $run_PFitter

  push @table_column_headers, "cluster.call";
  
  if( $run_profillic ) {
    push @table_column_headers, "profillic.clusters";
  } # End if $run_profillic

  my $table_header = join( "\t", @table_column_headers );
  print OUTPUT_TABLE_FH $table_header, "\n";

  # After removing/fixing hypermutated sequences and removing recombined sequences, run InSites (online) and get the informative:private site ratio stat; also run PhyML (locally) to get the mean pairwise diversity statistic, and also cluster the informative sites subalignments and compute the full consensus sequences for each cluster.
  my $id_string = "";
  my $output_path_dir_for_input_fasta_file;
  my $R_output;
  my $input_fasta_file;
  my $original_short_name_nosuffix_stripped;
  my %final_input_fasta_file_short_names_by_original_short_name_stripped;
  while ( scalar( @input_fasta_files ) > 0 ) {
    $input_fasta_file = shift @input_fasta_files;
    if( $VERBOSE ) {
      print $input_fasta_file, "\n";
    }
    my ( $input_fasta_file_path, $input_fasta_file_short ) =
      ( $input_fasta_file =~ /^(.*?)\/([^\/]+)$/ );
    unless( $input_fasta_file_short ) {
      $input_fasta_file_short = $input_fasta_file;
      $input_fasta_file_path = ".";
    }
    my ( $input_fasta_file_short_nosuffix, $input_fasta_file_suffix ) =
      ( $input_fasta_file_short =~ /^([^\.]+)(\..+)?$/ );

    if( defined $output_path_dir ) {
      $output_path_dir_for_input_fasta_file = $output_path_dir;
    } else {
      $output_path_dir_for_input_fasta_file =
        $input_fasta_file_path . "/" . $input_fasta_file_short_nosuffix . "_hiv-founder-id_resultsDir";
    }
    # remove trailing "/"
    ( $output_path_dir_for_input_fasta_file ) = ( $output_path_dir_for_input_fasta_file =~ /^(.*[^\/])\/*$/ );

    if( $input_fasta_file_short_nosuffix =~ /^(.+)(?:(?:_[LR]H)|(?:_NFLG))$/ ) {
      ( $original_short_name_nosuffix_stripped ) =
        ( $input_fasta_file_short_nosuffix =~ /^(.+)(?:(?:_[LR]H)|(?:_NFLG))$/ );
    } else {
      $original_short_name_nosuffix_stripped =
        $input_fasta_file_short_nosuffix;
    }

    # Duplicate the input file, possibly modifying it.
    my $input_fasta_file_contents = path( $input_fasta_file )->slurp();
    ## HACK: make sure there are no bangs in the input file (since sometimes there are, right now).
    if( $input_fasta_file_contents =~ /\!/ ) {
      if( $VERBOSE ) {
        print( "Input file \"$input_fasta_file\" contains illegal characters \"!\"; changing them to gaps, saving to output directory.\n" );
        ## TODOL REMOVE
        #print( $output_path_dir_for_input_fasta_file );
      }
      $input_fasta_file_contents =~ s/\!/-/g;
    } else {
      if( $VERBOSE ) {
        print( "Copying input file \"$input_fasta_file\" to output directory.\n" );
      }
    }
    # Now write it out to a temporary location in the output dir.
    $input_fasta_file_path = $output_path_dir_for_input_fasta_file;
    $input_fasta_file = "${input_fasta_file_path}/${input_fasta_file_short}";
    if( $VERBOSE ) {
      if( $input_fasta_file_contents =~ /\!/ ) {
        print( "Writing out fixed input file \"$input_fasta_file\".." );
      } else {
        print( "Writing out copy of input file to \"$input_fasta_file\".." );
      }
    }
    if( $VERBOSE ) { print "Opening file \"$input_fasta_file\" for writing..\n"; }
    unless( open input_fasta_fileFH, ">$input_fasta_file" ) {
      warn "Unable to open output file \"$input_fasta_file\": $!\n";
      return 1;
    }
    print input_fasta_fileFH $input_fasta_file_contents;
    close( input_fasta_fileFH );
    
    if( ( $input_fasta_file_contents =~ /RH\|/ ) && ( $input_fasta_file_contents =~ /LH\|/ ) ) {
      ## SPECIAL CASE: If the input file contains left-half and right-half genomes, separate them into separate input files and do both.
      if( $VERBOSE ) {
        print( "Input file contains left- and right- halves of the genome; separating them.\n" );
      }
      my @new_input_fasta_files =
        splitFastaFileOnHeaderPatterns( $output_path_dir_for_input_fasta_file, $input_fasta_file, , "NFLG\|", "LH\|", "RH\|" ); # NOTE THAT THE ORDER MATTERS.  RH last is assumed below.
      if( $VERBOSE ) {
        print( "Queueing the new files ( ", join( ", ", @new_input_fasta_files ), " ).\n" );
      }
      unshift @input_fasta_files, @new_input_fasta_files;
      next; # That's all we do with the original input file.
    }

    # This is the original unaltered fasta file (path removed), printed to the table out:
    print OUTPUT_TABLE_FH $input_fasta_file_short;
  
    print "\nInput Fasta file: $input_fasta_file_short\n";
    # First fix/remove hypermutated sequences, using an implementation of the HYPERMUT 2.0 algorithm.
    if( $run_Hypermut ) {
      if( $VERBOSE ) {
          if( $fix_hypermutated_sequences ) {
            print "Calling R to fix hypermutated sequences..";
          } else {
            print "Calling R to remove hypermutated sequences..";
          }
      }
      $R_output = `export removeHypermutatedSequences_fixInsteadOfRemove="$fix_hypermutated_sequences"; export removeHypermutatedSequences_pValueThreshold="$hypermut2_pValueThreshold"; export removeHypermutatedSequences_inputFilename="$input_fasta_file"; export removeHypermutatedSequences_outputDir="$output_path_dir_for_input_fasta_file"; R -f removeHypermutatedSequences.R --vanilla --slave`;
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
        $input_fasta_file_path = $output_path_dir_for_input_fasta_file;
        if( $fix_hypermutated_sequences ) {
          $input_fasta_file_short = "${input_fasta_file_short_nosuffix}_fixHypermutatedSequences${input_fasta_file_suffix}";
        } else {
          $input_fasta_file_short = "${input_fasta_file_short_nosuffix}_removeHypermutatedSequences${input_fasta_file_suffix}";
        }
        ( $input_fasta_file_short_nosuffix, $input_fasta_file_suffix ) =
          ( $input_fasta_file_short =~ /^([^\.]+)(\..+)?$/ );
        $input_fasta_file = "${input_fasta_file_path}/${input_fasta_file_short}";
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
      my $RAP_result_stdout = `perl runRAPOnline.pl $extra_flags $input_fasta_file $output_path_dir_for_input_fasta_file`;
      if( $RAP_result_stdout =~ /Recombinants identified/ ) {
        my ( $RAP_output_file ) = ( $RAP_result_stdout =~ /Recombinants identified \((.+)\)\s*$/ );
        $R_output = `export removeRecombinedSequences_pValueThreshold="$RAP_pValueThreshold"; export removeRecombinedSequences_RAPOutputFile="$RAP_output_file"; export removeRecombinedSequences_inputFilename="$input_fasta_file"; export removeRecombinedSequences_outputDir="$output_path_dir_for_input_fasta_file"; R -f removeRecombinedSequences.R --vanilla --slave`;
        ## extract the number fixed/removed from the output
        ( $removed_recombined_sequences ) = ( $R_output =~ /^.*\[1\]\s*(\d+)\s*$/ );
        if( $VERBOSE ) {
          print( ".done." );
        }
        # Now use the output from that..
        $input_fasta_file_path = $output_path_dir_for_input_fasta_file;
        $input_fasta_file_short = "${input_fasta_file_short_nosuffix}_removeRecombinedSequences${input_fasta_file_suffix}";
        ( $input_fasta_file_short_nosuffix, $input_fasta_file_suffix ) =
          ( $input_fasta_file_short =~ /^([^\.]+)(\..+)?$/ );
        $input_fasta_file = "${input_fasta_file_path}/${input_fasta_file_short}";
      } else {
        if( $VERBOSE ) {
          print( ".done." );
        }
        $removed_recombined_sequences = 0;
      } # End if any recombinants were identified. .. else ..
      print( "The number of recombined sequences removed is: $removed_recombined_sequences\n" );

      print OUTPUT_TABLE_FH "\t", $removed_recombined_sequences;
    } # End if $run_RAP

    print "Fasta file for analysis: $input_fasta_file_short\n";
    # This is the file for analysis (again, path stripped); maybe altered from the original input file.
    print OUTPUT_TABLE_FH "\t", $input_fasta_file_short;

    ## TODO: REMOVE
    #print( "SETTING \$final_input_fasta_file_short_names_by_original_short_name_stripped{ $original_short_name_nosuffix_stripped } = $input_fasta_file_short_nosuffix\n" );
    if( exists( $final_input_fasta_file_short_names_by_original_short_name_stripped{ $original_short_name_nosuffix_stripped } ) ) {
      push @{ $final_input_fasta_file_short_names_by_original_short_name_stripped{ $original_short_name_nosuffix_stripped } }, $input_fasta_file_short_nosuffix;
    } else {
      $final_input_fasta_file_short_names_by_original_short_name_stripped{ $original_short_name_nosuffix_stripped } = [ $input_fasta_file_short_nosuffix ];
    }

    if( $runInSites_online ) {
      `perl runInSitesOnline.pl $extra_flags $input_fasta_file $output_path_dir_for_input_fasta_file`;
    } else {
      `perl runInSitesOffline.pl $extra_flags $input_fasta_file $output_path_dir_for_input_fasta_file`;
    }
    `perl getInSitesStat.pl $extra_flags ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_informativeSites.txt ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_privateSites.txt ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_inSitesRatioStat.txt`;
    my $in_sites_ratio = `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_inSitesRatioStat.txt`;
    ( $in_sites_ratio ) = ( $in_sites_ratio =~ /^\s*(\S+)\s*$/ ); # Strip trailing newline, and any other flanking whitespace.

    # Run phyML, get stats
    `perl runPhyML.pl $extra_flags ${input_fasta_file} ${output_path_dir_for_input_fasta_file}`;

    my $pairwise_diversity_stats = `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}.phylip_phyml_pwdiversity.txt`;
    my ( $num_seqs, $mean_diversity ) =
      ( $pairwise_diversity_stats =~ /Max\s+(\d+)\s+([e\-\.\d]+)\s+/ );
    unless( defined( $num_seqs ) ) {
      warn( "UH OH: $input_fasta_file\nGOT:\n$pairwise_diversity_stats\n" );
    }
    print "Number of sequences: $num_seqs\n";
    print "Mean pairwise diversity: $mean_diversity\n";
    print "Informative sites to private sites ratio: $in_sites_ratio\n"; # Newline is no longer on there

    print OUTPUT_TABLE_FH "\t", $num_seqs;
    print OUTPUT_TABLE_FH "\t", $mean_diversity;
    print OUTPUT_TABLE_FH "\t", $in_sites_ratio;

    # Now cluster the informative sites (only relevant if one or both of the above exceeds a threshold.
    my $force_one_cluster = $default_force_one_cluster;
    my $morgane_calls_one_cluster = 1;
    my $diversity_threshold_exceeded = 0;
    my $in_sites_ratio_threshold_exceeded = 0;
    my $in_sites_cluster_call = 1;
    if( $mean_diversity > $mean_diversity_threshold ) {
      $diversity_threshold_exceeded = 1;
      if( $in_sites_ratio > $in_sites_ratio_threshold ) {
        $in_sites_ratio_threshold_exceeded = 1;
        print( "DIVERSITY THRESHOLD EXCEEDED\n" );
        print( "RATIO THRESHOLD EXCEEDED TOO\n" );
        $morgane_calls_one_cluster = 0;
        $force_one_cluster = 0;
      } else {
        print( "diversity threshold exceeded\n" );
      }
    } elsif( $in_sites_ratio > $in_sites_ratio_threshold ) {
        $in_sites_ratio_threshold_exceeded = 1;
        print( "ratio threshold exceeded\n" );
    }
    print OUTPUT_TABLE_FH "\t", $diversity_threshold_exceeded;
    print OUTPUT_TABLE_FH "\t", $in_sites_ratio_threshold_exceeded;
    ## "is.one.cluster":
    print OUTPUT_TABLE_FH "\t", $morgane_calls_one_cluster;

    if( $morgane_calls_one_cluster ) {
      print "Number of founders estimated using informative:private sites ratio and diversity thresholding is: 1\n";
      $in_sites_cluster_call = 1; # go with Morgane's call if it is 1.
    } else {
      print "Number of founders estimated using informative:private sites ratio and diversity thresholding is: greater than 1\n";
    }

    ## Now run PoissonFitter.
    my ( $PFitter_lambda, $PFitter_se, $PFitter_nseq, $PFitter_nbases, $PFitter_mean_hd, $PFitter_max_hd, $PFitter_chi_sq_stat, $PFitter_chi_sq_df, $PFitter_chi_sq_p_value ) = 0;
    my $PFitter_time_est_and_ci = "";
    my $PFitter_time_est = 0;
    my $PFitter_time_ci_low = 0;
    my $PFitter_time_ci_high = 0;
    my $is_poisson = 1;
    my $is_starlike = 1;
    my ( $Bayesian_PFitter_lambda_est, $Bayesian_PFitter_lambda_ci_low, $Bayesian_PFitter_lambda_ci_high ) = 0;
    my ( $Bayesian_PFitter_days_est, $Bayesian_PFitter_days_ci_low, $Bayesian_PFitter_days_ci_high ) = 0;
    my ( $DS_PFitter_lambda_est, $DS_PFitter_lambda_ci_low, $DS_PFitter_lambda_ci_high ) = 0;
    my ( $DS_PFitter_days_est, $DS_PFitter_days_ci_low, $DS_PFitter_days_ci_high ) = 0;
    my ( $DS_PFitter_distance_mean, $DS_PFitter_distance_ci_low, $DS_PFitter_distance_ci_high ) = 0;
    my $DS_PFitter_fitstext = "FITS (for very acute infection).";
    my $DS_PFitter_fits = 1;
    my $DS_PFitter_assertion_low = 1.5;
    my $DS_PFitter_assertion_high = 2.5;
    my $DS_PFitter_R = "1.0";;
    if( $run_PFitter ) {
      if( $mean_diversity == 0 ) {
        if( $VERBOSE ) {
          print "\nNOT running PoissonFitter because there is no variation whatsoever, and PFitter doesn't handle that case.\n";
        }
        print "PoissonFitter Determination: Degenerate Phylogeny (all seqs are the same)\n";
        print "PoissonFitter Poisson Fit: NA\n";
        print "DS Poisson Fit: $DS_PFitter_fitstext (R=$DS_PFitter_R).\n";
        print "Average distance to nearest Poisson CDF (2.5%, 97.5% quantiles): $DS_PFitter_distance_mean ($DS_PFitter_distance_ci_low, $DS_PFitter_distance_ci_high)\n";
        print "PoissonFitter Poisson time estimate (95\% CI): $PFitter_time_est ($PFitter_time_ci_low, $PFitter_time_ci_high)\n";
        #print "\n$PFitter_fitter_stats_raw\n";
        print "DS Poisson time estimate (95\% CI): $DS_PFitter_days_est ($DS_PFitter_days_ci_low, $DS_PFitter_days_ci_high)\n";
        print "Bayesian Poisson time estimate (95\% CI): $Bayesian_PFitter_days_est ($Bayesian_PFitter_days_ci_low, $Bayesian_PFitter_days_ci_high)\n";
      } else {
        if( $VERBOSE ) {
          print "\nCalling R to run PoissonFitter..";
        }
        $R_output = `export runPoissonFitter_inputFilename="$input_fasta_file"; export runPoissonFitter_outputDir="$output_path_dir_for_input_fasta_file"; export runPoissonFitter_runDSPFitter="$run_DSPFitter"; R -f runPoissonFitter.R --vanilla --slave`;
        if( $VERBOSE ) {
          print( "\tdone.\n" );
        }
        ## TODO: REMOVE
        # print "GOT: $R_output\n";
        my $PFitter_fitter_stats_raw =
          `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_PoissonFitterDir/LOG_LIKELIHOOD.results.txt`;
        ( $PFitter_lambda, $PFitter_se, $PFitter_nseq, $PFitter_nbases, $PFitter_mean_hd, $PFitter_max_hd, $PFitter_time_est_and_ci, $PFitter_chi_sq_stat, $PFitter_chi_sq_df, $PFitter_chi_sq_p_value  ) =
          (
           $PFitter_fitter_stats_raw =~ /\n[^\t]+\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(\S+)\s*$/ );
        ( $PFitter_time_est, $PFitter_time_ci_low, $PFitter_time_ci_high ) =
          ( $PFitter_time_est_and_ci =~ /(\S+) \((\S+), (\S+)\)/ );

        $is_poisson = ( defined( $PFitter_chi_sq_p_value ) && ( $PFitter_chi_sq_p_value > 0.05 ) ) || 0;
        my $starlike_raw =
          `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_PoissonFitterDir/CONVOLUTION.results.txt`;
        if( $DEBUG ) {
          ## TODO: REMOVE
          # print "PoissonFitter RAW: $starlike_raw\n";
        }
        my ( $starlike_text ) = ( $starlike_raw =~ m/(FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY/ );
        $is_starlike = ( $starlike_text eq "FOLLOWS" );

        # DS results
        my $DSPFitter_fitter_stats_raw =
          `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_PoissonFitterDir/${input_fasta_file_short_nosuffix}_DSPFitter.out`;
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
          ( $DSPFitter_fitter_stats_raw =~ /^((?:DOES NOT )?FIT.+)$/m );
        $DS_PFitter_fits =
          ( ( $DS_PFitter_fitstext =~ /^FITS.+$/m ) ? "1" : "0" );
        ( $DS_PFitter_assertion_low, $DS_PFitter_assertion_high, $DS_PFitter_R ) =
          ( $DSPFitter_fitter_stats_raw =~ /There is .*evidence against the assertion that the Poisson rate between sequences is between (\S+) and (\S+) times the rate of sequences to the consensus \(R = (\S+)\)/ );

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
        print "DS Poisson Fit: $DS_PFitter_fitstext (R=$DS_PFitter_R).\n";
        print "Average distance to nearest Poisson CDF (2.5%, 97.5% quantiles): $DS_PFitter_distance_mean ($DS_PFitter_distance_ci_low, $DS_PFitter_distance_ci_high)\n";
        print "PoissonFitter Poisson time estimate (95\% CI): $PFitter_time_est ($PFitter_time_ci_low, $PFitter_time_ci_high)\n";
        #print "\n$PFitter_fitter_stats_raw\n";
        print "DS Poisson time estimate (95\% CI): $DS_PFitter_days_est ($DS_PFitter_days_ci_low, $DS_PFitter_days_ci_high)\n";
        print "Bayesian Poisson time estimate (95\% CI): $Bayesian_PFitter_days_est ($Bayesian_PFitter_days_ci_low, $Bayesian_PFitter_days_ci_high)\n";
      } # End if( $mean_diversity > 0 );
      print OUTPUT_TABLE_FH "\t", $PFitter_lambda;
      print OUTPUT_TABLE_FH "\t", $PFitter_se;
      print OUTPUT_TABLE_FH "\t", $PFitter_nseq;
      print OUTPUT_TABLE_FH "\t", $PFitter_nbases;
      print OUTPUT_TABLE_FH "\t", $PFitter_mean_hd;
      print OUTPUT_TABLE_FH "\t", $PFitter_max_hd;
      print OUTPUT_TABLE_FH "\t", $PFitter_time_est;
      print OUTPUT_TABLE_FH "\t", $PFitter_time_ci_low;
      print OUTPUT_TABLE_FH "\t", $PFitter_time_ci_high;
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
    } # End if( $run_PFitter )

    my $tmp_extra_flags = $extra_flags;
    if( $force_one_cluster ) {
      $tmp_extra_flags .= "-f ";
    }
    my $num_clusters = `perl clusterInformativeSites.pl $tmp_extra_flags $input_fasta_file ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_informativeSites.txt $output_path_dir_for_input_fasta_file`;
    # There might be extra text in there.  Look for the telltale "[1]" output
    #warn "GOT $num_clusters\n";
    ( $num_clusters ) = ( $num_clusters =~ /\[1\] (\d+?)\s*/ );
    if( !$morgane_calls_one_cluster ) {
      $in_sites_cluster_call = $num_clusters;
    }

    # Print out the number of clusters
    print "Number of clusters found when clustering informative sites: $num_clusters\n";
    print "Number of founders estimated by the Informative Sites method: $in_sites_cluster_call\n";

    if( $run_PFitter ) {
      if( $in_sites_cluster_call == 1 ) {
        ## Avoid NA in the table output.  Multifounder results default to single-founder results.
        print OUTPUT_TABLE_FH "\t", $PFitter_lambda;
        print OUTPUT_TABLE_FH "\t", $PFitter_se;
        print OUTPUT_TABLE_FH "\t", $PFitter_nseq;
        print OUTPUT_TABLE_FH "\t", $PFitter_nbases;
        print OUTPUT_TABLE_FH "\t", $PFitter_mean_hd;
        print OUTPUT_TABLE_FH "\t", $PFitter_max_hd;
        print OUTPUT_TABLE_FH "\t", $PFitter_time_est;
        print OUTPUT_TABLE_FH "\t", $PFitter_time_ci_low;
        print OUTPUT_TABLE_FH "\t", $PFitter_time_ci_high;
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
      } else {
       die unless( $num_clusters > 1 );
        ## Now run PoissonFitter on the clusters.
        if( $VERBOSE ) {
          print "Calling R to run MultiFounderPoissonFitter..";
        }
       $R_output = `export runMultiFounderPoissonFitter_inputFilenamePrefix="${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}"; export runMultiFounderPoissonFitter_outputDir="$output_path_dir_for_input_fasta_file"; export runMultiFounderPoissonFitter_suffixPattern=""; export runMultiFounderPoissonFitter_runDSPFitter="TRUE"; R -f runMultiFounderPoissonFitter.R --vanilla --slave`;
       if( $VERBOSE ) {
         print( "\tdone.\n" );
       }
       my $multifounder_PFitter_fitter_stats_raw =
         `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_MultiFounderPoissonFitterDir/LOG_LIKELIHOOD.results.txt`;
       my ( $multifounder_PFitter_lambda, $multifounder_PFitter_se, $multifounder_PFitter_nseq, $multifounder_PFitter_nbases, $multifounder_PFitter_mean_hd, $multifounder_PFitter_max_hd, $multifounder_PFitter_time_est_and_ci, $multifounder_PFitter_chi_sq_stat, $multifounder_PFitter_chi_sq_df, $multifounder_PFitter_chi_sq_p_value  ) =
         (
          $multifounder_PFitter_fitter_stats_raw =~ /\n[^\t]+\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(\S+)\s*$/ );
       my ( $multifounder_PFitter_time_est, $multifounder_PFitter_time_ci_low, $multifounder_PFitter_time_ci_high ) =
         ( $multifounder_PFitter_time_est_and_ci =~ /(\S+) \((\S+), (\S+)\)/ );
       my $multifounder_is_poisson = ( defined( $multifounder_PFitter_chi_sq_p_value ) && ( $multifounder_PFitter_chi_sq_p_value > 0.05 ) ) || 0;
        ## NOTE THAT the convolution is not set up to handle multi-founder data because the convolution should be done within each founder; so for now we just exclude these results.  TODO: implement multi-founder version of the convolution.
        # my $multifounder_starlike_raw =
        #   `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_MultiFounderPoissonFitterDir/CONVOLUTION.results.txt`;
        # if( $DEBUG ) {
        #   ## TODO: REMOVE
        #   #print "MULTI-FOUNDER PoissonFitter RAW: $multifounder_starlike_raw\n";
        # }
        # my ( $multifounder_starlike_text ) = ( $multifounder_starlike_raw =~ m/(FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY/ );
        # my $multifounder_is_starlike = ( $multifounder_starlike_text eq "FOLLOWS" );

       # DS results
       my $multifounder_DSPFitter_fitter_stats_raw =
         `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_MultiFounderPoissonFitterDir/${input_fasta_file_short_nosuffix}_DSPFitter.out`;
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
          ( $multifounder_DSPFitter_fitter_stats_raw =~ /^((?:DOES NOT )?FIT.+)$/m );
       my $multifounder_DS_PFitter_fits = "0";
       if( $multifounder_DS_PFitter_fitstext =~ /^FITS.+$/m ) {
         $multifounder_DS_PFitter_fits = "1";
       }
       my ( $multifounder_DS_PFitter_assertion_low, $multifounder_DS_PFitter_assertion_high, $multifounder_DS_PFitter_R ) =
          ( $multifounder_DSPFitter_fitter_stats_raw =~ /There is .*evidence against the assertion that the Poisson rate between sequences is between (\S+) and (\S+) times the rate of sequences to the consensus \(R = (\S+)\)/ );

        # print "Multi-Founder PoissonFitter Determination: ";
        # if( $multifounder_is_starlike ) {
        #   print "Star-Like Phylogenies within clusters";
        # } else {
        #   print "Non-Star-Like Phylogenies within clusters";
       # }
        print "Multi-Founder PoissonFitter Poisson Fit: ";
        if( $multifounder_is_poisson ) {
          print "OK\n";
        } else {
          print "BAD (p = $multifounder_PFitter_chi_sq_p_value)\n";
        }
        print "Multi-Founder DS Poisson Fit: $multifounder_DS_PFitter_fitstext (R=$multifounder_DS_PFitter_R).\n";
        print "Multi-Founder Average distance to nearest Poisson CDF (2.5%, 97.5% quantiles): $multifounder_DS_PFitter_distance_mean ($multifounder_DS_PFitter_distance_ci_low, $multifounder_DS_PFitter_distance_ci_high)\n";
        print "Multi-Founder PoissonFitter Poisson time estimate (95\% CI): $multifounder_PFitter_time_est ($multifounder_PFitter_time_ci_low, $multifounder_PFitter_time_ci_high)\n";
        print "Multi-Founder DS Poisson time estimate (95\% CI): $multifounder_DS_PFitter_days_est ($multifounder_DS_PFitter_days_ci_low, $multifounder_DS_PFitter_days_ci_high)\n";
        print "Multi-Founder Bayesian Poisson time estimate (95\% CI): $multifounder_Bayesian_PFitter_days_est ($multifounder_Bayesian_PFitter_days_ci_low, $multifounder_Bayesian_PFitter_days_ci_high)\n";
        #print "\n$multifounder_PFitter_fitter_stats_raw\n";
      
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_lambda;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_se;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_nseq;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_nbases;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_mean_hd;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_max_hd;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_time_est;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_time_ci_low;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_time_ci_high;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_chi_sq_stat;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_chi_sq_df;
        print OUTPUT_TABLE_FH "\t", $multifounder_PFitter_chi_sq_p_value;
        print OUTPUT_TABLE_FH "\t", $multifounder_is_poisson;
#        print OUTPUT_TABLE_FH "\t", $multifounder_is_starlike;
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
      } # End if $num_clusters > 1
    } # End if $run_PFitter

    print OUTPUT_TABLE_FH "\t", $in_sites_cluster_call;

    ## Now try it the more profillic way.  This is a hybrid approach that makes profiles only of the informative sites.
    if( $run_profillic ) {
      if( $force_one_cluster ) {
        # Avoid NA in the table output.  Just print that the estimated number of founders is 1.  Which it is.
        print OUTPUT_TABLE_FH "\t", 1;
      } else {
        my $alignment_profiles_output_file = "${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_profileToAlignmentProfile.alignmentprofs";
        if( $VERBOSE ) {
          print "Running Profillic..\n";
        }
        my $alignment_profiles_output_files_list_file = "${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_profillic_AlignmentProfilesList.txt";
        `perl runProfillic.pl $extra_flags $input_fasta_file $alignment_profiles_output_files_list_file $output_path_dir_for_input_fasta_file`;

        if( $VERBOSE ) {
          print "Clustering..\n";
        }
        $R_output = `perl clusterProfillicAlignmentProfiles.pl $extra_flags $input_fasta_file $alignment_profiles_output_files_list_file $output_path_dir_for_input_fasta_file`;
        my ( $num_profillic_clusters ) = ( $R_output =~ /^.*\[1\]\s*(\d+)\s*$/ );
        # Print out the number of clusters
        print "Number of clusters found using profillic: $num_profillic_clusters\n";
        if( $morgane_calls_one_cluster ) {
          print "Number of founders estimated by the Informative Sites Profillic method: 1\n";
        } else {
          print "Number of founders estimated by the Informative Sites Profillic method: $num_profillic_clusters\n";
        }
        print OUTPUT_TABLE_FH "\t", $num_profillic_clusters;
      }
    } # End if $run_profillic

    print OUTPUT_TABLE_FH "\n";
    
    # If this is a half-genome dataset, should we call
    # runMultiFounderPoissonFitter to put together the two datasets
    # for better Poisson estimation?
    if( $run_PFitter ) {
      my $did_it = 0;
      if( $input_fasta_file_short =~ /_RH/ ) {
        my $nflg_version = $input_fasta_file;
        $nflg_version =~ s/_RH/_NFLG/;
        my $lh_version = $input_fasta_file;
        $lh_version =~ s/_RH/_LH/;
        # Do it only if there are at least the two; do it after "RH".
        # NOTE THAT WE ASSUME RH IS LAST. SEE ABOVE.
        if( ( -e $lh_version ) || ( -e $nflg_version ) ) {
          # OK, do it.
          $did_it = 1;
          ## Now run PoissonFitter on the regions.
          if( $VERBOSE ) {
            print "Calling R to run MultiRegionPoissonFitter..";
          }
          # NOTE that this strips off any modifiers like remove/fix HypermutatedSequences and removeRecombinedSequences: (this is dealt with in runMultiFounderPoissonFitter.R)
          my ( $input_fasta_file_very_short ) =
            ( $input_fasta_file_short =~ /^(.+)_RH/ );
          #print "\$input_fasta_file_very_short: $input_fasta_file_very_short\n"; 
          if( !exists ( $final_input_fasta_file_short_names_by_original_short_name_stripped{ $input_fasta_file_very_short } ) ) {
            die( "\$final_input_fasta_file_short_names_by_original_short_name_stripped{ $input_fasta_file_very_short } doesn't exist!" );
          }
          if( $VERBOSE ) {
            print "[Combining these sequences: ", join( ",", @{ $final_input_fasta_file_short_names_by_original_short_name_stripped{ $input_fasta_file_very_short } } ), "]..";
          }
          
          my @file_suffixes = map { ( $_ ) = ( $_ =~ /^$input_fasta_file_very_short(.+)$/ ); $_ } @{ $final_input_fasta_file_short_names_by_original_short_name_stripped{ $input_fasta_file_very_short } };
          my $suffix_pattern = join( "\.fasta|", @file_suffixes ) . "\.fasta";
          # print "\$input_fasta_file_very_short: $input_fasta_file_very_short\n";
          # print "\$suffix_pattern: $suffix_pattern\n";
          $R_output = `export runMultiFounderPoissonFitter_inputFilenamePrefix='${output_path_dir_for_input_fasta_file}/${input_fasta_file_very_short}'; export runMultiFounderPoissonFitter_suffixPattern='$suffix_pattern'; export runMultiFounderPoissonFitter_outputDir='$output_path_dir_for_input_fasta_file'; export runMultiFounderPoissonFitter_runDSPFitter='TRUE'; R -f runMultiFounderPoissonFitter.R --vanilla --slave`;
          # print $R_output;
          if( $VERBOSE ) {
            print( "\tdone.\n" );
          }
          my $multi_region_PFitter_fitter_stats_raw =
            `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_very_short}_MultiRegionPoissonFitterDir/LOG_LIKELIHOOD.results.txt`;
          my ( $multi_region_PFitter_lambda, $multi_region_PFitter_se, $multi_region_PFitter_nseq, $multi_region_PFitter_nbases, $multi_region_PFitter_mean_hd, $multi_region_PFitter_max_hd, $multi_region_PFitter_time_est_and_ci, $multi_region_PFitter_chi_sq_stat, $multi_region_PFitter_chi_sq_df, $multi_region_PFitter_chi_sq_p_value  ) =
            (
             $multi_region_PFitter_fitter_stats_raw =~ /\n[^\t]+\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(\S+)\s*$/ );
         my ( $multi_region_PFitter_time_est, $multi_region_PFitter_time_ci_low, $multi_region_PFitter_time_ci_high ) =
           ( $multi_region_PFitter_time_est_and_ci =~ /(\S+) \((\S+), (\S+)\)/ );
          my $multi_region_is_poisson = defined( $multi_region_PFitter_chi_sq_p_value ) && ( $multi_region_PFitter_chi_sq_p_value > 0.05 );
          ## NOTE THAT the convolution is not set up to handle multi-region data because the convolution should be done within each region; so for now we just exclude these results.  TODO: implement multi-region version of the convolution.
          # my $multi_region_starlike_raw =
          #   `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_MultiRegionPoissonFitterDir/CONVOLUTION.results.txt`;
          # if( $DEBUG ) {
          #   ## TODO: REMOVE
          #   #print "MULTI-REGION PoissonFitter RAW: $multi_region_starlike_raw\n";
          # }
          # my ( $multi_region_starlike_text ) = ( $multi_region_starlike_raw =~ m/(FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY/ );
          # my $multi_region_is_starlike = ( $multi_region_starlike_text eq "FOLLOWS" );

         # DS results
         my $multi_region_DSPFitter_fitter_stats_raw =
           `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_very_short}_MultiRegionPoissonFitterDir/${input_fasta_file_very_short}_DSPFitter.out`;
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
            ( $multi_region_DSPFitter_fitter_stats_raw =~ /^((?:DOES NOT )?FIT.+)$/m );
         my $multi_region_DS_PFitter_fits =
            ( ( $multi_region_DS_PFitter_fitstext =~ /^FITS.+$/m ) ? "1" : "0" );
         my ( $multi_region_DS_PFitter_assertion_low, $multi_region_DS_PFitter_assertion_high, $multi_region_DS_PFitter_R ) =
            ( $multi_region_DSPFitter_fitter_stats_raw =~ /There is .*evidence against the assertion that the Poisson rate between sequences is between (\S+) and (\S+) times the rate of sequences to the consensus \(R = (\S+)\)/ );
          
          # print "Multi-Region PoissonFitter Determination: ";
          # if( $multi_region_is_starlike ) {
          #   print "Star-Like Phylogenies within clusters";
          # } else {
          #   print "Non-Star-Like Phylogenies within clusters";
          # }
          print "\nInput fasta file: ${input_fasta_file_very_short}${input_fasta_file_suffix}\n";
          print "Multi-Region Poisson Fit: ";
          if( $multi_region_is_poisson ) {
            print "OK\n";
          } else {
            print "BAD (p = $multi_region_PFitter_chi_sq_p_value)\n";
          }
        print "Multi-Region DS Poisson Fit: $multi_region_DS_PFitter_fitstext (R=$multi_region_DS_PFitter_R).\n";
        print "Multi-Region Average distance to nearest Poisson CDF (2.5%, 97.5% quantiles): $multi_region_DS_PFitter_distance_mean ($multi_region_DS_PFitter_distance_ci_low, $multi_region_DS_PFitter_distance_ci_high)\n";
        print "Multi-Region PFitter Poisson time estimate (95\% CI): $multi_region_PFitter_time_est ($multi_region_PFitter_time_ci_low, $multi_region_PFitter_time_ci_high)\n";
        print "Multi-Region DS Poisson time estimate (95\% CI): $multi_region_DS_PFitter_days_est ($multi_region_DS_PFitter_days_ci_low, $multi_region_DS_PFitter_days_ci_high)\n";
        print "Multi-Region Bayesian Poisson time estimate (95\% CI): $multi_region_Bayesian_PFitter_days_est ($multi_region_Bayesian_PFitter_days_ci_low, $multi_region_Bayesian_PFitter_days_ci_high)\n";
          #print "\n$multi_region_PFitter_fitter_stats_raw\n";
          
          print OUTPUT_TABLE_FH "${input_fasta_file_very_short}${input_fasta_file_suffix}";
          print OUTPUT_TABLE_FH "\t", "NA"; # fixed-hypermut or removed-hypermut
          print OUTPUT_TABLE_FH "\t", "NA"; # removed-recomb
          print OUTPUT_TABLE_FH "\t", "NA"; # file
          print OUTPUT_TABLE_FH "\t", "NA"; # num.seqs
          print OUTPUT_TABLE_FH "\t", "NA"; # diversity
          print OUTPUT_TABLE_FH "\t", "NA"; # inf.to.priv.ratio
          print OUTPUT_TABLE_FH "\t", "NA"; # exceeds.diversity.threshold
          print OUTPUT_TABLE_FH "\t", "NA"; # exceeds.ratio.threshold
          print OUTPUT_TABLE_FH "\t", "NA"; # is.one.founder
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_lambda;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_se;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_nseq;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_nbases;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_mean_hd;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_max_hd;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_time_est;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_time_ci_low;
          print OUTPUT_TABLE_FH "\t", $multi_region_PFitter_time_ci_high;
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
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.PFitter.lambda
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.PFitter.se
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.PFitter.nseq
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.PFitter.nbases
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.PFitter.mean.hd
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.PFitter.max.hd
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.PFitter.time.est
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.PFitter.time.ci.low
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.PFitter.time.ci.high
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.PFitter.chi.sq.stat
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.PFitter.chi.sq.df
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.PFitter.chi.sq.p.value
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.is.poisson
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.Bayesian.PFitter.lambda.est
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.Bayesian.PFitter.lambda.ci.low
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.Bayesian.PFitter.lambda.ci.high
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.Bayesian.PFitter.days.est
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.Bayesian.PFitter.days.ci.low
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.Bayesian.PFitter.days.ci.high
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.DS.PFitter.lambda.est
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.DS.PFitter.lambda.ci.low
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.DS.PFitter.lambda.ci.high
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.DS.PFitter.days.est
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.DS.PFitter.days.ci.low
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.DS.PFitter.days.ci.high
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.DS.PFitter.distance.est
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.DS.PFitter.distance.ci.low
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.DS.PFitter.distance.ci.high
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.DS.PFitter.assertion.low
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.DS.PFitter.assertion.high
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.DS.PFitter.fits
          print OUTPUT_TABLE_FH "\t", "NA"; # multifounder.DS.PFitter.R
          print OUTPUT_TABLE_FH "\t", "NA"; # cluster.call
          if( $run_profillic ) {
            print OUTPUT_TABLE_FH "\t", "NA"; # "profillic.clusters"
          } # End if $run_profillic
          print OUTPUT_TABLE_FH "\n";
        } # End if this is the RH one, and if there is also an LH one, run MultiRegionPoissonFitter.
      } else {  # if this is an RH or LH one, consider combining for MultiRegionPoissonFitter. .. else ..
        if( $input_fasta_file_short =~ /_LH/ ) {
          ## I've lazily assumed that if there are _LH and _NFLG then there are also _RH.  That's the only case we miss because we don't care about _LH alone (since we only need to merge if there are multiple, and we've already handled the cases with _RH).
          ## HERE WE CHECK THIS ASSUMPTION.
          my $nflg_version = $input_fasta_file;
          $nflg_version =~ s/_LH/_NFLG/;
          my $rh_version = $input_fasta_file;
          $rh_version =~ s/_LH/_RH/;
          if( -e $nflg_version ) {
            die unless( -e $rh_version );
          }
        }
      } # End if this is an RH or LH one, consider combining for MultiRegionPoissonFitter. .. else ..

    } # End if $run_PFitter

  } # End foreach $input_fasta_file
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


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
##      Try: mkdir rv217_1W_gold_standard-hiv-founder-id_resultDir/; perl ./identify_founders.pl -O rv217_1W_gold_standard-hiv-founder-id_resultDir/ ~/src/from-git/projects/tholzman/MorgansFounderIDMethod/rv217_1W_gold_standard.list > rv217_1W_gold_standard-hiv-founder-id_resultDir/identify-founders.out
##      Or: rm -r rv217_1W_gold_standard-hiv-founder-id_-f_resultDir/; mkdir rv217_1W_gold_standard-hiv-founder-id_-f_resultDir/; perl ./identify_founders.pl -fV -O rv217_1W_gold_standard-hiv-founder-id_-f_resultDir/ ~/src/from-git/projects/tholzman/MorgansFounderIDMethod/rv217_1W_gold_standard.list > rv217_1W_gold_standard-hiv-founder-id_-f_resultDir/identify-founders.out
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
use vars qw( $opt_D $opt_V $opt_o $opt_O $opt_P $opt_R $opt_f $opt_I $opt_s );
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
    print "\tidentify_founders [-DV] [-PRfI] [-(o|O) <output_dir>] <input_fasta_file1 or file_listing_input_files> [<input_fasta_file2> ...] \n";
    exit;
  }

  # This means -D, -o, -O, -V, ... are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_o is an optional directory to put the output; default depends on input filename.
  # opt_O is just like opt_o.  Same thing.
  # opt_V means be verbose.
  # opt_P means skip (do not run) profillic
  # opt_R means skip (do not run) RAP (alignment-internal recombination detection program)
  # opt_f means fix hypermutated sequences, instead of removing them.
  # opt_I means run inSites online (instead of offline)
  # opt_s means be slow, ie run everything even when the diversity and insites thresholds are not exceeded.
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_o, $opt_O, $opt_P, $opt_R, $opt_f, $opt_I, $opt_s ) = ();
  if( not getopts('DVo:O:PRfIs') ) {
    identify_founders_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;

  my $run_profillic = !$opt_P;
  my $run_RAP = !$opt_R;
  my $fix_hypermutated_sequences = $opt_f || 0;
  my $runInSites_online = $opt_I || 0;
  my $be_slow = $opt_s || 0;
  my $run_PFitter = 1;

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

    if( 1 ) {
      my $input_fasta_file_contents = path( $input_fasta_file )->slurp();
      ## HACK: make sure there are no bangs in the input file (since sometimes there are, right now).
      if( $input_fasta_file_contents =~ /\!/ ) {
        if( $VERBOSE ) {
          print( "Input file \"$input_fasta_file\" contains illegal characters \"!\"; changing them to gaps\n" );
          ## TODOL REMOVE
          print( $output_path_dir_for_input_fasta_file );
        }
        $input_fasta_file_contents =~ s/\!/-/g;

        # Now write it out to a temporary location in the output dir.
        $input_fasta_file_path = $output_path_dir_for_input_fasta_file;
        $input_fasta_file = "${input_fasta_file_path}/${input_fasta_file_short}";
        if( $VERBOSE ) {
          print( "Writing out fixed input file \"$input_fasta_file\".." );
        }
        if( $VERBOSE ) { print "Opening file \"$input_fasta_file\" for writing..\n"; }
        unless( open input_fasta_fileFH, ">$input_fasta_file" ) {
            warn "Unable to open output file \"$input_fasta_file\": $!\n";
            return 1;
          }
        print input_fasta_fileFH $input_fasta_file_contents;
        close( input_fasta_fileFH );
      }
      if( ( $input_fasta_file_contents =~ /RH\|/ ) && ( $input_fasta_file_contents =~ /LH\|/ ) ) {
        ## SPECIAL CASE: If the input file contains left-half and right-half genomes, separate them into separate input files and do both.
        my @new_input_fasta_files =
          splitFastaFileOnHeaderPatterns( $output_path_dir_for_input_fasta_file, $input_fasta_file, , "NFLG\|", "LH\|", "RH\|" ); # NOTE THAT THE ORDER MATTERS.  RH last is assumed below.
        if( $VERBOSE ) {
          print( "Queueing the new files ( ", join( ", ", @new_input_fasta_files ), " ).\n" );
        }
        unshift @input_fasta_files, @new_input_fasta_files;
        next; # That's all we do with the original input file.
      }
    }

    # First fix/remove hypermutated sequences, using an implementation of the HYPERMUT 2.0 algorithm.
    if( $VERBOSE ) {
        if( $fix_hypermutated_sequences ) {
          print "Calling R to fix hypermutated sequences..";
        } else {
          print "Calling R to remove hypermutated sequences..";
        }
    }
    $R_output = `export removeHypermutatedSequences_fixInsteadOfRemove="$fix_hypermutated_sequences"; export removeHypermutatedSequences_pValueThreshold="$hypermut2_pValueThreshold"; export removeHypermutatedSequences_inputFilename="$input_fasta_file"; export removeHypermutatedSequences_outputDir="$output_path_dir_for_input_fasta_file"; R -f removeHypermutatedSequences.R --vanilla --slave`;
    ## extract the number fixed/removed from the output
    #$R_output = gsub( "^.*\\[1\\]\\s*(\\d+)", "\\1", $R_output );
#    if( $VERBOSE ) {
        if( $fix_hypermutated_sequences ) {
          print( "The number of hypermutated sequences fixed is: $R_output" );
        } else {
          print( "The number of hypermutated sequences removed is: $R_output" );
        }
#    }
    # Now use the output from that..
    $input_fasta_file_path = $output_path_dir_for_input_fasta_file;
    if( $fix_hypermutated_sequences ) {
      $input_fasta_file_short = "${input_fasta_file_short_nosuffix}_fixHypermutatedSequences${input_fasta_file_suffix}";
    } else {
      $input_fasta_file_short = "${input_fasta_file_short_nosuffix}_removeHypermutatedSequences${input_fasta_file_suffix}";
    }
    ( $input_fasta_file_short_nosuffix, $input_fasta_file_suffix ) =
      ( $input_fasta_file_short =~ /^([^\.]+)(\..+)?$/ );
    $input_fasta_file = "${input_fasta_file_path}/${input_fasta_file_short}";
    
    if( $VERBOSE ) {
      print ".done.\n";
    }

    if( $run_RAP ) {
      ## Run RAP on LANL, which gives individual
      ## sequences that are recombinants of other individual sequences,
      ## allowing those to be flagged for removal just like hypermutated
      ## ones.
      if( $VERBOSE ) {
        print "Running RAP at LANL to compute recombined sequences..";
      }
      my $RAP_result_stdout = `perl runRAPOnline.pl $extra_flags $input_fasta_file $output_path_dir_for_input_fasta_file`;
      if( $RAP_result_stdout =~ /Recombinants identified/ ) {
        my ( $RAP_output_file ) = ( $RAP_result_stdout =~ /Recombinants identified \((.+)\)\s*$/ );
        $R_output = `export removeRecombinedSequences_pValueThreshold="$RAP_pValueThreshold"; export removeRecombinedSequences_RAPOutputFile="$RAP_output_file"; export removeRecombinedSequences_inputFilename="$input_fasta_file"; export removeRecombinedSequences_outputDir="$output_path_dir_for_input_fasta_file"; R -f removeRecombinedSequences.R --vanilla --slave`;
        ## extract the number fixed/removed from the output
        my ( $removed_recombined_sequences ) = ( $R_output =~ /^.*\[1\]\s*(\d+)/ );
        if( $VERBOSE ) {
          print( "\tdone. The number of recombined sequences removed is: $removed_recombined_sequences\n" );
        }
        # Now use the output from that..
        $input_fasta_file_path = $output_path_dir_for_input_fasta_file;
        $input_fasta_file_short = "${input_fasta_file_short_nosuffix}_removeRecombinedSequences${input_fasta_file_suffix}";
        ( $input_fasta_file_short_nosuffix, $input_fasta_file_suffix ) =
          ( $input_fasta_file_short =~ /^([^\.]+)(\..+)?$/ );
        $input_fasta_file = "${input_fasta_file_path}/${input_fasta_file_short}";
      } else {
        if( $VERBOSE ) {
          print( "\tdone. The number of recombined sequences removed is: 0\n" );
        }
      } # End if any recombinants were identified. .. else ..
    } # End if $run_RAP

    ## TODO: REMOVE
    print( "SETTING \$final_input_fasta_file_short_names_by_original_short_name_stripped{ $original_short_name_nosuffix_stripped } = $input_fasta_file_short_nosuffix\n" );
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
    my $in_sites_stat = `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_inSitesRatioStat.txt`;

    # Run phyML, get stats
    `perl runPhyML.pl $extra_flags ${input_fasta_file} ${output_path_dir_for_input_fasta_file}`;

    my $pairwise_diversity_stats = `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}.phylip_phyml_pwdiversity.txt`;
    my ( $num_seqs, $mean_diversity ) =
      ( $pairwise_diversity_stats =~ /Max\s+(\d+)\s+([e\-\.\d]+)\s+/ );
    unless( defined( $num_seqs ) ) {
      warn( "UH OH: $input_fasta_file\nGOT:\n$pairwise_diversity_stats\n" );
    }
    print "Input fasta file: $input_fasta_file_short\n";
    print "Number of sequences: $num_seqs\n";
    print "Mean pairwise diversity: $mean_diversity\n";
    print "Informative sites to private sites ratio: $in_sites_stat"; # Newline is already on there

    # Now cluster the informative sites (only relevant if one or both of the above exceeds a threshold.
    my $force_one_cluster = $default_force_one_cluster;
    my $morgane_calls_one_cluster = 1;
    if( $mean_diversity > $mean_diversity_threshold ) {
      if( $in_sites_stat > $in_sites_ratio_threshold ) {
        print( "DIVERSITY THRESHOLD EXCEEDED\n" );
        print( "RATIO THRESHOLD EXCEEDED TOO\n" );
        $morgane_calls_one_cluster = 0;
        $force_one_cluster = 0;
      } else {
        print( "diversity threshold exceeded\n" );
      }
    } elsif( $in_sites_stat > $in_sites_ratio_threshold ) {
        print( "ratio threshold exceeded\n" );
    }

    if( $morgane_calls_one_cluster ) {
      print "\nNumber of founders estimated using Morgane's informative:private sites ratio and diversity threshold is: 1\n";
    } else {
      print "\nNumber of founders estimated using Morgane's informative:private sites ratio and diversity threshold is: greater than 1\n";
    }

    ## Now run PoissonFitter.
    if( $run_PFitter ) {
      if( $mean_diversity == 0 ) {
        if( $VERBOSE ) {
          print "\nNOT running PoissonFitter because there is no variation whatsoever, and PFitter doesn't handle that case.\n";
        }
      } else {
        if( $VERBOSE ) {
          print "\nCalling R to run PoissonFitter..";
        }
        $R_output = `export runPoissonFitter_inputFilename="$input_fasta_file"; export runPoissonFitter_outputDir="$output_path_dir_for_input_fasta_file"; R -f runPoissonFitter.R --vanilla --slave`;
        if( $VERBOSE ) {
          print( "\tdone.\n" );
        }
        my $poisson_fitter_stats_raw =
          `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_PoissonFitterDir/LOG_LIKELIHOOD.results.txt`;
        my ( $poisson_time_est_and_ci, $poisson_fit_stat ) = ( $poisson_fitter_stats_raw =~ /\n[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t[^\t]+\t[^\t]+\t(\S+)\s*$/ );
        my $is_poisson = defined( $poisson_fit_stat ) && ( $poisson_fit_stat > 0.05 );
        my $starlike_raw =
          `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_PoissonFitterDir/CONVOLUTION.results.txt`;
        if( $DEBUG ) {
          ## TODO: REMOVE
          # print "PoissonFitter RAW: $starlike_raw\n";
        }
        my ( $starlike_text ) = ( $starlike_raw =~ m/(FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY/ );
        my $is_starlike = ( $starlike_text eq "FOLLOWS" );
        print "PoissonFitter Determination: ";
        if( $is_starlike ) {
          print "Star-Like Phylogeny";
        } else {
          print "Non-Star-Like Phylogeny";
        }
        print "\nPoisson Fit: ";
        if( $is_poisson ) {
          print "OK";
        } else {
          print "BAD (p = $poisson_fit_stat)";
        }
        print "\nPoisson time estimate (95\% CI): $poisson_time_est_and_ci\n";
        #print "\n$poisson_fitter_stats_raw\n";
      } # End if( $mean_diversity > 0 );
    } # End if( $run_PFitter )

    # Should we call runMultiFounderPoissonFitter to put together the two datasets for better Poisson estimation?
    if( $input_fasta_file_short =~ /_RH/ ) {
      my $nflg_version = $input_fasta_file;
      $nflg_version =~ s/_RH/_NFLG/;
      my $lh_version = $input_fasta_file;
      $lh_version =~ s/_RH/_LH/;
      # Do it only if there are at least the two; do it after "RH".
      # NOTE THAT WE ASSUME RH IS LAST. SEE ABOVE.
      if( ( $input_fasta_file_short =~ /_RH/ ) && ( ( -e $lh_version ) || ( -e $nflg_version ) ) ) {
        # OK, do it.

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
        print "Combining these sequences: ", join( ",", @{ $final_input_fasta_file_short_names_by_original_short_name_stripped{ $input_fasta_file_very_short } } ), "\n";
        
        my @file_suffixes = map { ( $_ ) = ( $_ =~ /^$input_fasta_file_very_short(.+)$/ ); $_ } @{ $final_input_fasta_file_short_names_by_original_short_name_stripped{ $input_fasta_file_very_short } };
        my $suffix_pattern = join( "\$\.fasta|", @file_suffixes ) . "\$\.fasta";
        print "\$input_fasta_file_very_short: $input_fasta_file_very_short\n";
        print "\$suffix_pattern: $suffix_pattern\n";
        $R_output = `export runMultiFounderPoissonFitter_inputFilenamePrefix='$input_fasta_file_very_short'; export runMultiFounderPoissonFitter_suffixPattern='$suffix_pattern'; export runMultiFounderPoissonFitter_outputDir='$output_path_dir_for_input_fasta_file'; R -f runMultiFounderPoissonFitter.R --vanilla --slave`;
        if( $VERBOSE ) {
          print( "\tdone.\n" );
        }
        my $multi_region_poisson_fitter_stats_raw =
          `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_MultiRegionPoissonFitterDir/LOG_LIKELIHOOD.results.txt`;
        my ( $multi_region_poisson_time_est_and_ci, $multi_region_poisson_fit_stat ) =
          ( $multi_region_poisson_fitter_stats_raw =~ /\n[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t[^\t]+\t[^\t]+\t(\S+)\s*$/ );
        my $multi_region_is_poisson = defined( $multi_region_poisson_fit_stat ) && ( $multi_region_poisson_fit_stat > 0.05 );
        ## NOTE THAT the convolution is not set up to handle multi-region data because the convolution should be done within each region; so for now we just exclude these results.  TODO: implement multi-region version of the convolution.
        # my $multi_region_starlike_raw =
        #   `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_MultiRegionPoissonFitterDir/CONVOLUTION.results.txt`;
        # if( $DEBUG ) {
        #   ## TODO: REMOVE
        #   #print "MULTI-REGION PoissonFitter RAW: $multi_region_starlike_raw\n";
        # }
        # my ( $multi_region_starlike_text ) = ( $multi_region_starlike_raw =~ m/(FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY/ );
        # my $multi_region_is_starlike = ( $multi_region_starlike_text eq "FOLLOWS" );
        # print "Multi-Region PoissonFitter Determination: ";
        # if( $multi_region_is_starlike ) {
        #   print "Star-Like Phylogenies within clusters";
        # } else {
        #   print "Non-Star-Like Phylogenies within clusters";
        # }
        print "\nMulti-Region Poisson Fit: ";
        if( $multi_region_is_poisson ) {
          print "OK";
        } else {
          print "BAD (p = $multi_region_poisson_fit_stat)";
        }
        print "\nMulti-Region Poisson time estimate (95\% CI): $multi_region_poisson_time_est_and_ci\n";
        #print "\n$multi_region_poisson_fitter_stats_raw\n";
      } # End if this is the RH one, and if there is also an LH one, run MultiRegionPoissonFitter.
    } # End if this is an RH or LH one, consider combining for PoissonFitter.
    elsif( $input_fasta_file_short =~ /_LH/ ) {
      ## I've lazily assumed that if there are _LH and _NFLG then there are also _RH.  That's the only case we miss because we don't care about _LH alone (since we only need to merge if there are multiple, and we've already handled the cases with _RH).
      ## HERE WE CHECK THIS ASSUMPTION.
      my $nflg_version = $input_fasta_file;
      $nflg_version =~ s/_LH/_NFLG/;
      my $rh_version = $input_fasta_file;
      $rh_version =~ s/_LH/_RH/;
      if( -e $nflg_version ) {
        stopifnot( -e $rh_version );
      }
    }

    my $tmp_extra_flags = $extra_flags;
    if( $force_one_cluster ) {
      $tmp_extra_flags .= "-f ";
    }
    my $num_clusters = `perl clusterInformativeSites.pl $tmp_extra_flags $input_fasta_file ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_informativeSites.txt $output_path_dir_for_input_fasta_file`;
    # There might be extra text in there.  Look for the telltale "[1]" output
    #warn "GOT $num_clusters\n";
    ( $num_clusters ) = ( $num_clusters =~ /\[1\] (\d+?)\s*/ );
    
    # Print out the number of clusters
    print "\nNumber of founders estimated by clustering informative sites: $num_clusters\n\n";

    if( $num_clusters > 1 ) {
      ## Now run PoissonFitter on the clusters.
      if( $VERBOSE ) {
        print "Calling R to run MultiFounderPoissonFitter..";
      }
      $R_output = `export runMultiFounderPoissonFitter_inputFilenamePrefix="$input_fasta_file"; export runMultiFounderPoissonFitter_outputDir="$output_path_dir_for_input_fasta_file"; R -f runMultiFounderPoissonFitter.R --vanilla --slave`;
      if( $VERBOSE ) {
        print( "\tdone.\n" );
      }
      my $multi_founder_poisson_fitter_stats_raw =
        `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_MultiFounderPoissonFitterDir/LOG_LIKELIHOOD.results.txt`;
      my ( $multi_founder_poisson_time_est_and_ci, $multi_founder_poisson_fit_stat ) =
        ( $multi_founder_poisson_fitter_stats_raw =~ /\n[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t[^\t]+\t[^\t]+\t(\S+)\s*$/ );
      my $multi_founder_is_poisson = defined( $multi_founder_poisson_fit_stat ) && ( $multi_founder_poisson_fit_stat > 0.05 );
      ## NOTE THAT the convolution is not set up to handle multi-founder data because the convolution should be done within each founder; so for now we just exclude these results.  TODO: implement multi-founder version of the convolution.
      # my $multi_founder_starlike_raw =
      #   `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_MultiFounderPoissonFitterDir/CONVOLUTION.results.txt`;
      # if( $DEBUG ) {
      #   ## TODO: REMOVE
      #   #print "MULTI-FOUNDER PoissonFitter RAW: $multi_founder_starlike_raw\n";
      # }
      # my ( $multi_founder_starlike_text ) = ( $multi_founder_starlike_raw =~ m/(FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY/ );
      # my $multi_founder_is_starlike = ( $multi_founder_starlike_text eq "FOLLOWS" );
      # print "Multi-Founder PoissonFitter Determination: ";
      # if( $multi_founder_is_starlike ) {
      #   print "Star-Like Phylogenies within clusters";
      # } else {
      #   print "Non-Star-Like Phylogenies within clusters";
      # }
      print "\nMulti-Founder Poisson Fit: ";
      if( $multi_founder_is_poisson ) {
        print "OK";
      } else {
        print "BAD (p = $multi_founder_poisson_fit_stat)";
      }
      print "\nMulti-Founder Poisson time estimate (95\% CI): $multi_founder_poisson_time_est_and_ci\n";
      #print "\n$multi_founder_poisson_fitter_stats_raw\n";
    } # End if $num_clusters > 1

    ## Now try it the more profillic way.
    if( $run_profillic ) {
      if( !$force_one_cluster ) {
        my $alignment_profiles_output_file = "${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_profileToAlignmentProfile.alignmentprofs";
        print "Running Profillic..\n";
        my $alignment_profiles_output_files_list_file = "${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_profillic_AlignmentProfilesList.txt";
        `perl runProfillic.pl $extra_flags $input_fasta_file $alignment_profiles_output_files_list_file $output_path_dir_for_input_fasta_file`;
        
        print "Clustering..\n";
        my $num_profillic_clusters = `perl clusterProfillicAlignmentProfiles.pl $extra_flags $input_fasta_file $alignment_profiles_output_files_list_file $output_path_dir_for_input_fasta_file`;
        # Print out the number of clusters
        print "Number of founders estimated using Profillic: $num_profillic_clusters\n";
      }
    } # End if $run_profillic

  } # End foreach $input_fasta_file

  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  return 0;
} # identify_founders(..)

identify_founders( @ARGV );

1;


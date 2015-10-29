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
##      Try: mkdir rv217_1W_gold_standard-hiv-founder-id_resultDir/; perl ./identify_founders.pl -O rv217_1W_gold_standard-hiv-founder-id_resultDir/ ~/src/from-git/projects/tholzman/MorgansFounderIDMethod/rv217_1W_gold_standard.list > rv217_1W_gold_standard-hiv-founder-id_resultDir/identify-founders.out
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
use vars qw( $opt_D $opt_V $opt_o $opt_O $opt_P $opt_R $opt_f $opt_I );
use vars qw( $VERBOSE $DEBUG );

sub identify_founders {
  @ARGV = @_;

  sub identify_founders_usage {
    print "\tidentify_founders [-DV] [-PRfI] [-(o|O) <output_dir>] <input_fasta_file1 or file_listing_input_files> [<input_fasta_file2> ...] \n";
    exit;
  }

  # This means -D, -o, -O, -V, -P, -R are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_o is an optional directory to put the output; default depends on input filename.
  # opt_O is just like opt_o.  Same thing.
  # opt_V means be verbose.
  # opt_P means skip (do not run) profillic
  # opt_R means skip (do not run) RAP (alignment-internal recombination detection program)
  # opt_f means fix hypermutated sequences, instead of removing them.
  # opt_I means run inSites online (instead of offline)
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_o, $opt_O, $opt_P, $opt_R, $opt_f, $opt_I ) = ();
  if( not getopts('DVo:O:PRfI') ) {
    identify_founders_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;

  my $run_profillic = !$opt_P;
  my $run_RAP = !$opt_R;
  my $fix_hypermutated_sequences = $opt_f || 0;
  my $runInSites_online = $opt_I || 0;
  
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

  # After removing/fixing hypermutated sequences and removing recombined sequences, run InSites (online) and get the informative:private site ratio stat; also run PhyML (locally) to get the mean pairwise diversity statistic, and also cluster the informative sites subalignments and compute the full consensus sequences for each cluster.
  my $id_string = "";
  my $output_path_dir_for_input_fasta_file;
  my $R_output;
  foreach my $input_fasta_file ( @input_fasta_files ) {
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

    ## HACK: make sure there are no bangs in the input file (since there are, right now).
    if( 1 ) {
      my $input_fasta_file_contents = path( $input_fasta_file )->slurp();
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
    }

    # First fix/remove hypermutated sequences, using an implementation of the HYPERMUT 2.0 algorithm.
    if( $VERBOSE ) {
        if( $fix_hypermutated_sequences ) {
          print "Calling R to fix hypermutated sequences..";
        } else {
          print "Calling R to remove hypermutated sequences..";
        }
    }
    my $hypermut2_pValueThreshold = 0.1; # Matches Abrahams 2009
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
      my $RAP_pValueThreshold = 0.0007; # Appears to be the suggestion from the output file "(summaryTable)"'s column header, which reads "Pvalues<0.0007".
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
    my $mean_diversity_threshold = 0.001;
    #my $in_sites_ratio_threshold = 0.85; # For Abrahams and RV217
    my $in_sites_ratio_threshold = 0.33; # For caprisa002
    my $force_one_cluster = 1;
    if( $mean_diversity > $mean_diversity_threshold ) {
      if( $in_sites_stat > $in_sites_ratio_threshold ) {
        print( "DIVERSITY THRESHOLD EXCEEDED\n" );
        print( "RATIO THRESHOLD EXCEEDED TOO\n" );
        $force_one_cluster = 0;
      } else {
        print( "diversity threshold exceeded\n" );
      }
    } elsif( $in_sites_stat > $in_sites_ratio_threshold ) {
        print( "ratio threshold exceeded\n" );
    }

    ## Now run PoissonFitter.
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
      my $multi_founder_starlike_raw =
        `cat ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_MultiFounderPoissonFitterDir/CONVOLUTION.results.txt`;
      if( $DEBUG ) {
        ## TODO: REMOVE
        #print "MULTI-FOUNDER PoissonFitter RAW: $multi_founder_starlike_raw\n";
      }
      my ( $multi_founder_starlike_text ) = ( $multi_founder_starlike_raw =~ m/(FOLLOWS|DOES NOT FOLLOW) A STAR-PHYLOGENY/ );
      my $multi_founder_is_starlike = ( $multi_founder_starlike_text eq "FOLLOWS" );
      print "Multi-Founder PoissonFitter Determination: ";
      if( $multi_founder_is_starlike ) {
        print "Star-Like Phylogenies within clusters";
      } else {
        print "Non-Star-Like Phylogenies within clusters";
      }
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


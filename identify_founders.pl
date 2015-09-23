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
##      Try: perl ./identify_founders.pl -O rv217_1W_gold_standard-hiv-founder-id_resultDir/ -V rv217_1W_gold_standard.list > rv217_1W_gold_standard-hiv-founder-id_resultDir/identify-founders.out
##      
###******************************************************************************

use Getopt::Std; # for getopts
use File::Path qw( make_path );
use Path::Tiny;

use strict;
use vars qw( $opt_D $opt_V $opt_o $opt_O );
use vars qw( $VERBOSE $DEBUG );

sub identify_founders {
  @ARGV = @_;

  sub identify_founders_usage {
    print "\tidentify_founders [-DV] [-(o|O) <output_dir>] <input_fasta_file1 or file_listing_input_files> [<input_fasta_file2> ...] \n";
    exit;
  }

  # This means -D, -o, -O, -V are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_o is an optional directory to put the output; default depends on input filename.
  # opt_O is just like opt_o.  Same thing.
  # opt_V means be verbose.
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_o, $opt_O ) = ();
  if( not getopts('DVo:O:') ) {
    identify_founders_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;

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
  }

  my $extra_flags = "";
  if( $DEBUG ) {
    $extra_flags .= "-D ";
  }
  if( $VERBOSE ) {
    $extra_flags .= "-V ";
  }
  # Run InSites (online) and get the informative:private site ratio stat; also run PhyML (locally) to get the mean pairwise diversity statistic, and also cluster the informative sites subalignments and compute the full consensus sequences for each cluster.
  my $id_string = "";
  my $output_path_dir_for_input_fasta_file;
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

    `perl runInSitesOnline.pl $extra_flags $input_fasta_file $output_path_dir_for_input_fasta_file`;
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
    my $in_sites_ratio_threshold = 0.20;
    my $force_one_cluster = 0;
    if( $mean_diversity > $mean_diversity_threshold ) {
      print( "DIVERSITY THRESHOLD EXCEEDED\n" );
      $force_one_cluster = 1;
    }
    if( $in_sites_stat > $in_sites_ratio_threshold ) {
      print( "RATIO THRESHOLD EXCEEDED\n" );
      $force_one_cluster = 1;
    }
    my $tmp_extra_flags = $extra_flags;
    if( $force_one_cluster ) {
      $tmp_extra_flags .= "-f ";
    }
    my $num_clusters = `perl clusterInformativeSites.pl $tmp_extra_flags $input_fasta_file ${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_informativeSites.txt $output_path_dir_for_input_fasta_file`;
    # Print out the number of clusters
    print "Number of founders: $num_clusters\n\n";
  } # End foreach $input_fasta_file

  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  return 0;
} # identify_founders(..)

identify_founders( @ARGV );

1;


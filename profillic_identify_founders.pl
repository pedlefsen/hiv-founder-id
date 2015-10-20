#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) profillic_identify_founders
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##  Description:
##
##      Script for using Profillic (https://github.com/galosh/profillic) to
##      estimate the number of founders and perform sequence
##      reconstruction of the founders.
##
##      Note that this creates output files in subdirectories named
##      after the input fasta file name (unless you specify an output dir).
##
##      Or: mkdir Abrahams-2009aa-profillic_identify_founders_resultDir/; perl ./profillic_identify_founders.pl -V -O Abrahams-2009aa-profillic_identify_founders_resultDir/ Abrahams-2009aa/preparedFor_hiv-identify-founders.list > Abrahams-2009aa-profillic_identify_founders_resultDir/identify-founders.out 
###******************************************************************************

use Getopt::Std; # for getopts
use File::Path qw( make_path );
use Path::Tiny;

use strict;
use vars qw( $opt_D $opt_V $opt_o $opt_O $opt_f );
use vars qw( $VERBOSE $DEBUG );

sub profillic_identify_founders {
  @ARGV = @_;

  sub profillic_identify_founders_usage {
    print "\tprofillic_identify_founders [-DV] [-f] [-(o|O) <output_dir>] <input_fasta_file1 or file_listing_input_files> [<input_fasta_file2> ...] \n";
    exit;
  }

  # This means -D, -o, -O, -V, -P, -R are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_o is an optional directory to put the output; default depends on input filename.
  # opt_O is just like opt_o.  Same thing.
  # opt_V means be verbose.
  # opt_f means fix hypermutated sequences, instead of removing them.
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_o, $opt_O, $opt_f ) = ();
  if( not getopts('DVo:O:f') ) {
    profillic_identify_founders_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;

  my $fix_hypermutated_sequences = $opt_f || 0;
  
  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }
  
  my @input_fasta_files = @ARGV;
  unless( scalar @input_fasta_files ) { profillic_identify_founders_usage(); }

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

  # First remove/fix hypermutated sequences and remove recombined sequences
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
        $input_fasta_file_path . "/" . $input_fasta_file_short_nosuffix . "_profillic_identify_founders_resultsDir";
    }
    # remove trailing "/"
    ( $output_path_dir_for_input_fasta_file ) = ( $output_path_dir_for_input_fasta_file =~ /^(.*[^\/])\/*$/ );

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
    $R_output = gsub( "^.*\\[1\\]\\s*(\\d+)", "\\1", $R_output );
#    if( $VERBOSE ) {
        if( $fix_hypermutated_sequences ) {
          print( "The number of hypermutated sequences fixed is: $R_output" );
        } else {
          print( "The number of hypermutated sequences removed is: $R_output" );
        }
#    }
    # Now use the output from that..
    $input_fasta_file_path = $output_path_dir_for_input_fasta_file;
    $input_fasta_file_short = "${input_fasta_file_short_nosuffix}_removeHypermutatedSequences${input_fasta_file_suffix}";
    ( $input_fasta_file_short_nosuffix, $input_fasta_file_suffix ) =
      ( $input_fasta_file_short =~ /^([^\.]+)(\..+)?$/ );
    $input_fasta_file = "${input_fasta_file_path}/${input_fasta_file_short}";
    
    if( $VERBOSE ) {
      print ".done.\n";
    }

    ## Run RAP on LANL, which gives individual
    ## sequences that are recombinants of other individual sequences,
    ## allowing those to be flagged for removal just like hypermutated
    ## ones.
    my $RAP_pValueThreshold = 0.0007; # Appears to be the suggestion from the output file "(summaryTable)"'s column header, which reads "Pvalues<0.0007".
    if( $VERBOSE ) {
      print "Running RAP at LANL to compute recombined sequences..";
    }
    my $RAP_result_stdout = `perl runRAPOnline.pl $extra_flags $input_fasta_file $output_path_dir_for_input_fasta_file`;
    if( $VERBOSE ) {
      print "\tdone. Got $RAP_result_stdout\n";
    }
    if( $RAP_result_stdout =~ /Recombinants identified \(/ ) {
      my $RAP_output_file = $output_path_dir . "/" . $input_fasta_file_short_nosuffix . "_RAP.txt";
      if( $VERBOSE ) {
        print "Calling R to remove recombined sequences..";
      }
      $R_output = `export removeRecombinedSequences_pValueThreshold="$RAP_pValueThreshold"; export removeRecombinedSequences_RAPOutputFile="$RAP_output_file"; export removeRecombinedSequences_inputFilename="$input_fasta_file"; export removeRecombinedSequences_outputDir="$output_path_dir_for_input_fasta_file"; R -f removeRecombinedSequences.R --vanilla --slave`;
      ## extract the number fixed/removed from the output
      $R_output = gsub( "^.*\\[1\\]\\s*(\\d+)", "\\1", $R_output );
      if( $VERBOSE ) {
        print( "\tdone. The number of recombined sequences removed is: $R_output" );
      }
      # Now use the output from that..
      $input_fasta_file_path = $output_path_dir_for_input_fasta_file;
      $input_fasta_file_short = "${input_fasta_file_short_nosuffix}_removeRecombinedSequences${input_fasta_file_suffix}";
      ( $input_fasta_file_short_nosuffix, $input_fasta_file_suffix ) =
        ( $input_fasta_file_short =~ /^([^\.]+)(\..+)?$/ );
      $input_fasta_file = "${input_fasta_file_path}/${input_fasta_file_short}";
    } # End if any recombinants were identified.
   
    ## Now try it the more profillic way.
    if( $VERBOSE ) {
      print "Running Profillic..\n";
    }
    my $alignment_profiles_output_files_list_file = "${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_profillic_AlignmentProfilesList.txt";
    my $runProfillic_output =
      `perl runProfillicFromScratch.pl $tmp_extra_flags $input_fasta_file $alignment_profiles_output_files_list_file $output_path_dir_for_input_fasta_file`;

    my ( $self_entropy ) = ( $runProfillic_output =~ /Self Entropy: (\S+)/ );
    if( $VERBOSE ) {
      print( "\tdone. The self entropy of the profile is: $self_entropy" );
    }

    #print "Self Entropy: $self_entropy\n";
    #my $self_entropy_threshold = 
    my $force_one_cluster = 1;#( $self_entropy > $self_entropy_threshold );

    if( !$force_one_cluster ) {
        my $alignment_profiles_output_file = "${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_profileToAlignmentProfile.alignmentprofs";
        print "Clustering..\n";
        my $num_profillic_clusters = `perl clusterProfillicAlignmentProfiles.pl $tmp_extra_flags $input_fasta_file $alignment_profiles_output_files_list_file $output_path_dir_for_input_fasta_file`;
        # Print out the number of clusters
        print "Number of founders estimated using Profillic: $num_profillic_clusters\n";
    } # End if !$force_one_cluster

  } # End foreach $input_fasta_file

  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  return 0;
} # profillic_identify_founders(..)

profillic_identify_founders( @ARGV );

1;


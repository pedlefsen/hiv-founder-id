#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) clusterProfillicAlignmentProfiles
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##  Description:
##      Script for using Profile HMMs to cluster sequences.
##      By default it writes the output to
##      stdout; use -o to write to a file.  Note that this calls runProfillic.pl
##      to create a file of alignment profiles.  This file is then broken into
##      input-sequence-specfiic alignment profile files. And then this
##      calls out to the R script clusterProfillicAlignmentProfiles.R.
##      
###******************************************************************************

use Getopt::Std; # for getopts
use Text::Wrap; # for wrap
$Text::Wrap::columns = 72;# TODO: DEHACKIFY MAGIC #

use strict;
use vars qw( $opt_D $opt_V $opt_f );
use vars qw( $VERBOSE $DEBUG );


sub clusterProfillicAlignmentProfiles {
  @ARGV = @_;

  sub clusterProfillicAlignmentProfiles_usage {
    print "\tclusterProfillicAlignmentProfiles [-DVf] <input_fasta_filename> <profillic_alignment_profile_files_list_file> [<output dir>]\n";
    exit;
  }

  # This means -D, -V, and -f are ok, but nothin' else.
  # opt_f means don't actually cluster them; instead put them all into one cluster ("cluster 0").
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_f ) = ();
  if( not getopts('DVf') ) {
    clusterProfillicAlignmentProfiles_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;

  my $force_one_cluster = $opt_f;
  
  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }

  my $input_fasta_file = shift @ARGV || clusterProfillicAlignmentProfiles_usage();
  my ( $input_fasta_file_path, $input_fasta_file_short ) =
    ( $input_fasta_file =~ /^(.*?)\/([^\/]+)$/ );
  unless( $input_fasta_file_short ) {
    $input_fasta_file_short = $input_fasta_file;
    $input_fasta_file_path = ".";
  }

  my $alignment_profile_files_list_file = shift @ARGV || clusterProfillicAlignmentProfiles_usage();
  
  my $output_dir = shift @ARGV || $input_fasta_file_path;
  # Remove the trailing "/" if any
  if( defined( $output_dir ) ) {
    ( $output_dir ) = ( $output_dir =~ /^(.*[^\/])\/*$/ );
  }

  if( $VERBOSE ) { print "Output will be written in directory \"$output_dir\".."; }
  if( !-e $output_dir ) {
    `mkdir $output_dir`;
  }

  if( $VERBOSE ) {
    print "Calling R to cluster the alignment profiles..";
    if( $force_one_cluster ) {
      print( "Forcing one cluster..\n" );
    } else {
      print( "Clustering..\n" );
    }
  }
  my $R_output = `export clusterProfillicAlignmentProfiles_forceOneCluster="$force_one_cluster"; export clusterProfillicAlignmentProfiles_alignmentProfileFilesListFilename="$alignment_profile_files_list_file"; export clusterProfillicAlignmentProfiles_fastaFilename="$input_fasta_file"; export clusterProfillicAlignmentProfiles_outputDir="$output_dir"; R -f clusterProfillicAlignmentProfiles.R --vanilla --slave`;
  # The output has the number of clusters.
  print( $R_output );

  if( $VERBOSE ) {
    print "\t.done.\n";
  }
    
  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  if( $VERBOSE ) {
    print ".Done.\n";
  }
  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  return 0;
} # clusterProfillicAlignmentProfiles(..)

clusterProfillicAlignmentProfiles( @ARGV );

1;


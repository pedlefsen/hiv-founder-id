#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) evaluateFounderCalls
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##  Description:
##      Script for evaluating the founder sequence predictions.
##      
###******************************************************************************

use Getopt::Std; # for getopts
use File::Path qw( make_path );
use Path::Tiny;
use Statistics::Descriptive;
use Readonly;

use strict;
use vars qw( $opt_D $opt_V );
use vars qw( $VERBOSE $DEBUG );

sub evaluateFounderCalls {
  @ARGV = @_;

  sub evaluateFounderCalls_usage {
    print "\tevaluateFounderCalls [-DV] <predictions_ungapped_fasta_file> <ideal_ungapped_fasta_file> [<output_dir>]\n";
    exit;
  }

  # This means -D and -V are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # But first reset the opt vars.
  ( $opt_D, $opt_V ) = ();
  if( not getopts('DV') ) {
    evaluateFounderCalls_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;
  
  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }
  
  my $predictions_fasta_file = shift @ARGV || evaluateFounderCalls_usage();
  my ( $predictions_fasta_file_path, $predictions_fasta_file_short ) =
    ( $predictions_fasta_file =~ /^(.*?)\/([^\/]+)$/ );
  unless( $predictions_fasta_file_short ) {
    $predictions_fasta_file_short = $predictions_fasta_file;
    $predictions_fasta_file_path = ".";
  }
  my ( $predictions_fasta_file_short_nosuffix, $predictions_fasta_file_suffix ) =
    ( $predictions_fasta_file_short =~ /^([^\.]+)(\..+)?$/ );

  my $ideal_fasta_file = shift @ARGV || evaluateFounderCalls_usage();
  my ( $ideal_fasta_file_path, $ideal_fasta_file_short ) =
    ( $ideal_fasta_file =~ /^(.*?)\/([^\/]+)$/ );
  unless( $ideal_fasta_file_short ) {
    $ideal_fasta_file_short = $ideal_fasta_file;
    $ideal_fasta_file_path = ".";
  }
  my ( $ideal_fasta_file_short_nosuffix, $ideal_fasta_file_suffix ) =
    ( $ideal_fasta_file_short =~ /^([^\.]+)(\..+)?$/ );

  my $output_path_dir = shift @ARGV ||
    $predictions_fasta_file_path . "/" . $predictions_fasta_file_short_nosuffix . "_hiv-founder-id_resultsDir";
  # Remove the trailing "/" if any
  if( defined( $output_path_dir ) ) {
    ( $output_path_dir ) = ( $output_path_dir =~ /^(.*[^\/])\/*$/ );
  }
  #make_path( $output_path_dir );

  my $extra_flags = "";
  if( $DEBUG ) {
    $extra_flags .= "-D ";
  }
  if( $VERBOSE ) {
    $extra_flags .= "-V ";
  }
  
  my $combined_fasta_file = "${output_path_dir}/${predictions_fasta_file_short_nosuffix}_combinedWith_${ideal_fasta_file_short_nosuffix}${ideal_fasta_file_suffix}";
  `cat $ideal_fasta_file > $combined_fasta_file`;
  `cat $predictions_fasta_file >> $combined_fasta_file`;
  my $nucs_genecutter_path = "${output_path_dir}/${predictions_fasta_file_short_nosuffix}_combinedWith_${ideal_fasta_file_short_nosuffix}_allnucs";
  my $proteins_genecutter_path = "${output_path_dir}/${predictions_fasta_file_short_nosuffix}_combinedWith_${ideal_fasta_file_short_nosuffix}_allproteins";
  if( !( -e $nucs_genecutter_path ) || !( -e $proteins_genecutter_path ) ) {
    my $GeneCutter_result_stdout = `perl runGeneCutterOnline.pl $extra_flags $combined_fasta_file $output_path_dir`;
    if( $VERBOSE ) {
      print $GeneCutter_result_stdout;
    }
    # OK now unzip the files.
    `mkdir $nucs_genecutter_path`;
    `unzip ${output_path_dir}/${predictions_fasta_file_short_nosuffix}_combinedWith_${ideal_fasta_file_short_nosuffix}_allnucs.zip -d $nucs_genecutter_path`;
    
    `mkdir $proteins_genecutter_path`;
    `unzip ${output_path_dir}/${predictions_fasta_file_short_nosuffix}_combinedWith_${ideal_fasta_file_short_nosuffix}_allproteins.zip -d $proteins_genecutter_path`;
  }

  ## ERE I AM.  NOW CALL AN R PROGRAM WITH EACH OF THE OUTPUT FILES IN THOSE DIRS.
  my $R_output;
  my @proteins_genecutter_files = <"${proteins_genecutter_path}/*">;
  foreach my $protein_genecutter_file ( @proteins_genecutter_files ) {
    if( `grep Error $protein_genecutter_file` eq "" ) {
      if( $VERBOSE ) {
        print "Calling R to evaluate $protein_genecutter_file..";
      }
      $R_output = `export evaluateFounderCalls_isProteins="TRUE"; export evaluateFounderCalls_inputFilename="$protein_genecutter_file"; export evaluateFounderCalls_idealFilename="$ideal_fasta_file"; export evaluateFounderCalls_outputDir="$output_path_dir"; R -f evaluateFounderCalls.R --vanilla --slave`;
      if( $VERBOSE ) {
        print( $R_output );
      }
      ## get the output filename.
      my ( $statistics_file ) = ( $R_output =~ /^\[1\]\s*(\d+)\s*$/m );
      
    } else {
      # NO RESULT.
      warn( "No result in $protein_genecutter_file" );
    }
  }
  if( $VERBOSE ) {
    print( ".done\n" );
  }

  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  return 0;
} # evaluateFounderCalls(..)

######################################################################################

evaluateFounderCalls( @ARGV );

1;


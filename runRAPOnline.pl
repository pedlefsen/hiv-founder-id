#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) runRAPOnline
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##  Description:
##      Script for running RAP using the online LANL tools to determine which
##      sequences are recombinants of other sequences in a set.
##      
###******************************************************************************

use Getopt::Std; # for getopts
use File::Path qw( make_path );
use Path::Tiny;

# For screenscraping
use WWW::Mechanize;
use LWP::Simple;
#use HTTP::Request::Common qw(POST);
#use LWP::UserAgent;

use strict;
use vars qw( $opt_D $opt_V );
use vars qw( $VERBOSE $DEBUG );

sub runRAPOnline {
  @ARGV = @_;

  sub runRAPOnline_usage {
    print "\trunRAPOnline [-DV] <input_fasta_file> [<output_dir>]\n";
    exit;
  }

  # This means -D and -V are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # But first reset the opt vars.
  ( $opt_D, $opt_V ) = ();
  if( not getopts('DV') ) {
    runRAPOnline_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;
  
  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }
  
  my $input_fasta_file = shift @ARGV || runRAPOnline_usage();
  my ( $input_fasta_file_path, $input_fasta_file_short ) =
    ( $input_fasta_file =~ /^(.*?)\/([^\/]+)$/ );
  unless( $input_fasta_file_short ) {
    $input_fasta_file_short = $input_fasta_file;
    $input_fasta_file_path = ".";
  }
  my ( $input_fasta_file_short_nosuffix, $input_fasta_file_suffix ) =
    ( $input_fasta_file_short =~ /^([^\.]+)(\..+)?$/ );

  my $output_path_dir = shift @ARGV ||
    $input_fasta_file_path . "/" . $input_fasta_file_short_nosuffix . "_hiv-founder-id_resultsDir";
  # Remove the trailing "/" if any
  if( defined( $output_path_dir ) ) {
    ( $output_path_dir ) = ( $output_path_dir =~ /^(.*[^\/])\/*$/ );
  }

  my $R_output;

  # RAP seems to not be able to handle duplicate sequences.
  ## STEP 1: remove duplicate sequences.
  if( $VERBOSE ) {
    print "Calling R to create a version of the fasta file in which duplicate sequences are removed..";
  }
  $R_output = `export removeDuplicateSequencesFromAlignedFasta_inputFilename="$input_fasta_file"; export removeDuplicateSequencesFromAlignedFasta_outputDir="$output_path_dir"; R -f removeDuplicateSequencesFromAlignedFasta.R --vanilla --slave`;
  # The output has the file name of the consensus file.
  if( $DEBUG ) {
    print( "GOT: $R_output\n" );
  }
  if( $VERBOSE ) {
    print ".done.\n";
  }
  # Parse it to get the output fasta filename.
  my ( $fasta_file_no_duplicates ) = ( $R_output =~ /\"([^\"]+)\"/ );
  my ( $table_file_no_duplicates ) =
    ( $fasta_file_no_duplicates =~ /^(.+)$input_fasta_file_suffix$/ ) . ".tbl";
  if( $DEBUG ) {
    print "Fasta file with duplicates removed: $fasta_file_no_duplicates\n";
    print "Table of duplicates removed: $table_file_no_duplicates\n";
  }

  # RAP has a problem with certain characters in fasta headers. 
  ## STEP 2: Rename the seqs to just their numbers.
  if( $VERBOSE ) {
    print "Calling R to create a version of the fasta file in which seqs are numbered instead of named..";
  }
  $R_output = `export computeConsensusSequenceFromAlignedFasta_inputFilename="$fasta_file_no_duplicates"; export computeConsensusSequenceFromAlignedFasta_outputDir="$output_path_dir"; export computeConsensusSequenceFromAlignedFasta_includeFullAlignment="TRUE"; export  computeConsensusSequenceFromAlignedFasta_includeConsensus="FALSE"; export computeConsensusSequenceFromAlignedFasta_useSeqeunceNumbersAsNames="TRUE"; R -f computeConsensusSequenceFromAlignedFasta.R --vanilla --slave`;
  # The output has the file name of the consensus file.
  if( $DEBUG ) {
    print( "GOT: $R_output\n" );
  }
  if( $VERBOSE ) {
    print ".done.\n";
  }
  # Parse it to get the filename.
  my ( $fasta_file_readyForRAP ) = ( $R_output =~ /\"([^\"]+)\"/ );
  if( $DEBUG ) {
    print "Fasta file with consensus: $fasta_file_readyForRAP\n";
  }

  ## Set up the table mapping sequence names to numbers.
      if( $VERBOSE ) {
        print "Reading sequence names from file \"", $fasta_file_no_duplicates, "\"..";
      }
      my $fasta_file_no_duplicates_contents =
           path( $fasta_file_no_duplicates )->slurp();
      if( $VERBOSE ) {
        print ".done\n";
      }
      if( $DEBUG ) {
        #print $fasta_file_no_duplicates_contents;
      }
      my ( @seq_names ) =  ( $fasta_file_no_duplicates_contents =~ /\n?>[ \t]*(.+) *\n/g );

  my $mech = WWW::Mechanize->new( autocheck => 1 );
  $mech->get( "http://www.hiv.lanl.gov/content/sequence/RAP/RAP.html" );

  my $result = $mech->submit_form(
                                  form_name => 'input',
                                  fields    => { alignmentFile => $fasta_file_readyForRAP }
                       );
  
  my $content = $mech->content();
  if( $DEBUG ) {
    print "OK\n";#, \$content is $content\n";
  }
  my ( $no_recombinants ) = ( $content =~ /No recombinants found/ );
  my ( $format_error ) = ( $content =~ /'Format Conversion' probably did not recognize the format of your input alignment/ );
  if( $no_recombinants || $format_error ) {
    if( $VERBOSE ) {
      if( $format_error ) {
        print "ERROR: File format of $fasta_file_readyForRAP not recognized by RAP!\n";
      }
      if( $no_recombinants ) {
        print "No recombinants identified.\n";
      }
      select STDOUT;
      $| = $old_autoflush;
    }
    return 0;
  }

  # TODO: Save the highlighter plots?
  
  my ( $RAP_id ) = ( $content =~ /\/(\d+)\/summaryTable/ );
  unless( defined( $RAP_id ) ) {
    warn "No RAP_id in:\n$content\n";
  }
  my $RAP_output_file = $output_path_dir . "/" . $input_fasta_file_short_nosuffix . "_RAP.txt";
  my $RAP_output_file_contents = get "http://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/tmp/RAP/${RAP_id}/summaryTable";

  ## TODO: REMOVE
  #print "RAP ID: $RAP_id\n";

    ## Need to map the numbers back to names.
  my @RAP_output_file_lines = split( "\n", $RAP_output_file_contents );
  my ( $firstpart, $recombinant, $parents_before_split, $rest ) = undef;
  my @parents;
  for( my $line_i = 2; $line_i < scalar( @RAP_output_file_lines ); $line_i++ ) {
## TODO: REMOVE!
##    print "LINE $line_i: $RAP_output_file_lines[ $line_i ]\n";
    ( $firstpart, $recombinant, $parents_before_split, $rest ) =
      ( $RAP_output_file_lines[ $line_i ] =~ /^(Set \d+) (\d+) ([,\d]+) (.+)\s*$/ );
    @parents = split( ",", $parents_before_split );
    if( $DEBUG ) {
      print "PARENTS before split: $parents_before_split\n";
      print "PARENTS: ", join( ", ", @parents ), "\n";
    }
    $RAP_output_file_lines[ $line_i ] = $firstpart . ' ' . $seq_names[ $recombinant - 1 ] . ' ' . join( ",", @seq_names[ map { $_ - 1 } @parents ] ) . ' ' . $rest;
  }
  $RAP_output_file_contents = join( "\n", @RAP_output_file_lines ) . "\n";

        if( $VERBOSE ) {
          print( "Writing out RAP output file \"$RAP_output_file\".." );
        }
        if( $VERBOSE ) { print "Opening file \"$RAP_output_file\" for writing..\n"; }
        unless( open RAP_output_fileFH, ">$RAP_output_file" ) {
            warn "Unable to open output file \"$RAP_output_file\": $!\n";
            return 1;
          }
        print RAP_output_fileFH $RAP_output_file_contents;
        close( RAP_output_fileFH );

  if( $VERBOSE ) {
    print "Recombinants identified ($RAP_output_file)\n";
    select STDOUT;
    $| = $old_autoflush;
  }

  return 0;
} # runRAPOnline(..)

runRAPOnline( @ARGV );

1;


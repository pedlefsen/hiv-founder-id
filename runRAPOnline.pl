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

  my $mech = WWW::Mechanize->new( autocheck => 1 );
  $mech->get( "http://www.hiv.lanl.gov/content/sequence/RAP/RAP.html" );

  my $result = $mech->submit_form(
                                  form_name => 'input',
                                  fields    => { alignmentFile => $input_fasta_file }
                       );
  
  my $content = $mech->content();
  if( $DEBUG ) {
    print "OK\n";#, \$content is $content\n";
  }
  my ( $no_recombinants ) = ( $content =~ /No recombinants found/ );

  if( $no_recombinants ) {
    if( $VERBOSE ) {
      print "No recombinants identified.\n";
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


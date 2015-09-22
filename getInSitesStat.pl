#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) getInSitesStat
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##  Description:
##      Script for parsing output from the InformativeSites DiveIn
##      tool (at the Mullins lab web site) to calculate the stat used
##      to estimate the number of founders and perform sequence
##      reconstruction of the founders.  Takes two input files (first can be missing):
##      getInSitesStat 123456789012.txt 123456789012_priv.txt
##      By default it writes the output to stdout; use -o to write to a file.
##      If the first file does not exist (no informative sites) use "-" for first argument.
##      
###******************************************************************************

use Getopt::Std; # for getopts

use strict;
use vars qw( $opt_D $opt_V $opt_o $opt_O );
use vars qw( $VERBOSE $DEBUG );

sub getInSitesStat {
  @ARGV = @_;

  sub getInSitesStat_usage {
    print "\tgetInSitesStat [-DV] [-(o|O) <Output File>] <jobnum.txt> <jobnum_priv.txt> [<Output File]\n";
    exit;
  }

  # This means -D, -o, -O, -V are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_o is an optional filename to put the output.  Otherwise, STDOUT.
  # opt_O is just like opt_o except it'll overwrite the file if it exists.
  # opt_V means be verbose.
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_o, $opt_O ) = ();
  if( not getopts('DVo:O:') ) {
    getInSitesStat_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;
  
  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }

  ## Sometimes this file is missing, because there are no informative sites.  We take "-" to indicate this condition.
  my $informative_sites_file = shift @ARGV || getInSitesStat_usage();
  my ( $informative_sites_file_path, $informative_sites_file_short ) =
    ( $informative_sites_file =~ /^(.*?)\/([^\/]+)$/ );
  unless( $informative_sites_file_short ) {
    $informative_sites_file_short = $informative_sites_file;
    $informative_sites_file_path = ".";
  }
  my ( $informative_sites_file_short_nosuffix, $informative_sites_file_suffix ) =
    ( $informative_sites_file_short =~ /^([^\.]+)(\..+)?$/ );

  my $private_sites_file = shift @ARGV || getInSitesStat_usage();
  my ( $private_sites_file_path, $private_sites_file_short ) =
    ( $private_sites_file =~ /^(.*?)\/([^\/]+)$/ );
  unless( $private_sites_file_short ) {
    $private_sites_file_short = $private_sites_file;
    $private_sites_file_path = ".";
  }
  my ( $private_sites_file_short_nosuffix, $private_sites_file_suffix ) =
    ( $private_sites_file_short =~ /^([^\.]+)(\..+)?$/ );

  my $informative_no_gaps = 0; # default for when informative sites file is empty
  if( $informative_sites_file eq '-' ) {
    $informative_no_gaps = 0;
  } else {
    if( $VERBOSE ) { print "Opening file \"$informative_sites_file\".."; }
    
    unless( open( INFORMATIVE_SITES_FILE_FH, $informative_sites_file ) ) {
      warn "Unable to open informative_sites_file file \"$informative_sites_file\": $!";
      return 1;
    }
    
    if( $VERBOSE ) { print ".done.\n"; }
    
    if( $VERBOSE ) { print "Reading informative sites file.."; }
    my ( $line );
    while( $line = <INFORMATIVE_SITES_FILE_FH> ) {
    
      ( $line ) = ( $line =~ /^(.+?)\s*$/ );  # Chop and Chomp won't remove ^Ms
  
      if( $VERBOSE ) { print "."; }
      next unless $line =~ /^Alignment/;
      ( $informative_no_gaps ) = ( $line =~ /Alignment\s\d+\s(\d+)\s\d+/ );
      last;
    }
    if( $VERBOSE ) { print ".done.\n"; }
  
    if( $VERBOSE ) { print "Closing file \"$informative_sites_file\".."; }
    close INFORMATIVE_SITES_FILE_FH;
    if( $VERBOSE ) { print ".done.\n"; }
  } # End if $informative_sites_file eq '-' .. else ..  
  
  my $private_no_gaps;
  if( $VERBOSE ) { print "Opening file \"$private_sites_file\".."; }
  unless( open( PRIVATE_SITES_FILE_FH, $private_sites_file ) ) {
    warn "Unable to open private_sites_file file \"$private_sites_file\": $!";
    return 1;
  }
  if( $VERBOSE ) { print ".done.\n"; }
  
  if( $VERBOSE ) { print "Reading private sites file.."; }
  my ( $line );
  while( $line = <PRIVATE_SITES_FILE_FH> ) {
  
    ( $line ) = ( $line =~ /^(.+?)\s*$/ );  # Chop and Chomp won't remove ^Ms
  
    if( $VERBOSE ) { print "."; }
    next unless $line =~ /^Alignment/;
    ( $private_no_gaps ) = ( $line =~ /Alignment\s\d+\s(\d+)\s\d+/ );
    last;
  }
  if( $VERBOSE ) { print ".done.\n"; }
  
  if( $VERBOSE ) { print "Closing file \"$private_sites_file\".."; }
  close PRIVATE_SITES_FILE_FH;
  if( $VERBOSE ) { print ".done.\n"; }

  if( $VERBOSE ) {
    print "\$informative_no_gaps is $informative_no_gaps\n";
    print "\$private_no_gaps is $private_no_gaps\n";
  }
  # Named for Morgane Rolland, whose advice lead us to consider this statistic.
  my $morganes_stat = $informative_no_gaps / $private_no_gaps;
  if( $VERBOSE ) {
    print "\$morganes_stat is $morganes_stat\n";
  }

  my $output_file = $opt_o || $opt_O || shift @ARGV;
  if( $output_file ) {
  
    if( !$opt_O && -e $output_file ) {
      print "The file \"$output_file\" exists.  Use -O to overwrite.\n";
      exit;
    }
  
    if( $VERBOSE ) { print "Opening file \"$output_file\" for writing.."; }
  
    unless( open OUTPUT_FH, ">$output_file" ) {
      warn "Unable to open output file \"$output_file\": $!";
      return 1;
    }
  
    if( $VERBOSE ) { print ".done.\n"; }
  
  } else {
  
    unless( open OUTPUT_FH, ">&STDOUT" ) {
      warn "Unable to rename STDOUT: $!";
      return 1;
    }
  
  }

  if( $VERBOSE ) {
    print "Printing output..";
  }
  
  print OUTPUT_FH "$morganes_stat\n";
  
  if( $VERBOSE ) { print ".done.\n"; }
  
  if( $VERBOSE && $output_file ) { print "Closing file \"$output_file\".."; }
  close OUTPUT_FH;
  if( $VERBOSE && $output_file ) { print ".done.\n"; }

  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  return 0;
} # getInSitesStat(..)

getInSitesStat( @ARGV );

1;


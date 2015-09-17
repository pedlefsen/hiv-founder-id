#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) clusterInformativeSites
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##  Description:
##      Script for parsing the informativeSites output from the
##      InformativeSites DiveIn tool (at the Mullins lab web site) to
##      cluster the sequences.  By default it writes the output to
##      stdout; use -o to write to a file.  Note that this creates a
##      multiple alignment of just the informative sites and then
##      calls out to the R script clusterInformativeSites.R.
##      
###******************************************************************************

use Getopt::Std; # for getopts
use Text::Wrap; # for wrap
$Text::Wrap::columns = 72;# TODO: DEHACKIFY MAGIC #

use strict;
use vars qw( $opt_D $opt_V $opt_o $opt_O );
use vars qw( $VERBOSE $DEBUG );

sub clusterInformativeSites {
  @ARGV = @_;

  sub clusterInformativeSites_usage {
    print "\tclusterInformativeSites [-DV] [-(o|O) <Output File>] <informativeSites.txt> [<Output File>]\n";
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
    clusterInformativeSites_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;
  
  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }

  my $informative_sites_file = shift @ARGV || clusterInformativeSites_usage();
  my ( $informative_sites_file_path, $informative_sites_file_short ) =
    ( $informative_sites_file =~ /^(.*?)\/([^\/]+)$/ );
  unless( $informative_sites_file_short ) {
    $informative_sites_file_short = $informative_sites_file;
    $informative_sites_file_path = ".";
  }
  my ( $informative_sites_file_short_nosuffix, $informative_sites_file_suffix ) =
    ( $informative_sites_file_short =~ /^([^\.]+)(\..+)?$/ );

  if( $VERBOSE ) { print "Opening file \"$informative_sites_file\".."; }
  
  unless( open( INFORMATIVE_SITES_FILE_FH, $informative_sites_file ) ) {
    warn "Unable to open informative_sites_file file \"$informative_sites_file\": $!";
    return 1;
  }
  
  if( $VERBOSE ) { print ".done.\n"; }
  
  if( $VERBOSE ) { print "Reading informative sites file.."; }
  my ( $line );
  my %seqs;
  my @fields;
  my @seqorder;
  my @seq_as_list;
  my $seq;
  while( $line = <INFORMATIVE_SITES_FILE_FH> ) {
  
    ( $line ) = ( $line =~ /^(.+?)\s*$/ );  # Chop and Chomp won't remove ^Ms
  
    next unless $line;
    next if $line =~ /^\s/;
    last if $line =~ /^Alignment/;
  
    if( $VERBOSE ) { print "."; }
  
    @fields = split "\t", $line;
  
    if( $DEBUG ) {
      print $fields[ 0 ], ": ";
      print join( "\t", @fields[ 4..$#fields ] ), "\n";
    }
    push @seqorder, $fields[ 0 ];
    $seqs{ $fields[ 0 ] } = [ @fields[ 4..$#fields ] ];
  }
  if( $VERBOSE ) { print ".done.\n"; }
  
  if( $VERBOSE ) { print "Closing file \"$informative_sites_file\".."; }
  close INFORMATIVE_SITES_FILE_FH;
  if( $VERBOSE ) { print ".done.\n"; }

  # Ok, great.  Now, for all non-consensus entries, replace '.' with the consensus value.
  die unless( $seqorder[ 0 ] eq 'Consensus' );
  shift @seqorder; # Remove `Consensus`
  
  if( $DEBUG ) {
    print "Sequences: ", join( ", ", @seqorder ), "\n";
  }

  my @consensus_as_list = @{ $seqs{ "Consensus" } };
  my $alignment_length = scalar( @consensus_as_list );
  delete $seqs{ "Consensus" };

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
    print "Printing alignment of informative sites..";
  }
  my $seq_wrapped;
  foreach my $seq_name ( @seqorder ) {
    if( $DEBUG ) {
      print $seq_name, "\n";
    }
    @seq_as_list = @{ $seqs{ $seq_name } };

    die unless( scalar( @seq_as_list ) == $alignment_length );
    for( my $i = 0; $i < scalar( @seq_as_list ); $i++ ) {
      if( $seq_as_list[ $i ] eq '.' ) {
        $seq_as_list[ $i ] = $consensus_as_list[ $i ];
      }
    }
    $seq = join "", @seq_as_list;
    
    # add newlines / wraparound to $seq
    $seq_wrapped = wrap( '', '',  $seq );
    print OUTPUT_FH "> $seq_name\n$seq_wrapped\n";
  } # End foreac $seq_name
  
  if( $VERBOSE ) { print ".done.\n"; }
  
  if( $VERBOSE && $output_file ) { print "Closing file \"$output_file\".."; }
  close OUTPUT_FH;
  if( $VERBOSE && $output_file ) { print ".done.\n"; }

  if( $VERBOSE ) {
    print "Calling R to cluster informative sites..";
  }
  my $R_output = `export clusterInformativeSites_inputFilename="$output_file"; echo \$clusterInformativeSites_inputFilename; R -f clusterInformativeSites.R --vanilla --slave`;

  if( $VERBOSE ) {
    print ".done.\n";
  }
  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  return 0;
} # clusterInformativeSites(..)

clusterInformativeSites( @ARGV );

1;


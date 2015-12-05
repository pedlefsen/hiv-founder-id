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
use vars qw( $opt_D $opt_V $opt_f );
use vars qw( $VERBOSE $DEBUG );

sub clusterInformativeSites {
  @ARGV = @_;

  sub clusterInformativeSites_usage {
    print "\tclusterInformativeSites [-DVf] <input_fasta_filename> <informativeSites.txt> [<output dir>]\n";
    exit;
  }

  # This means -D, -V, and -f are ok, but nothin' else.
  # opt_f means don't actually cluster them; instead put them all into one cluster ("cluster 0").
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_f ) = ();
  if( not getopts('DVf') ) {
    clusterInformativeSites_usage();
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

  my $original_fasta_filename = shift @ARGV || clusterInformativeSites_usage();

  my $informative_sites_file = shift @ARGV || clusterInformativeSites_usage();
  my ( $informative_sites_file_path, $informative_sites_file_short ) =
    ( $informative_sites_file =~ /^(.*?)\/([^\/]+)$/ );
  unless( $informative_sites_file_short ) {
    $informative_sites_file_short = $informative_sites_file;
    $informative_sites_file_path = ".";
  }

  my $output_dir = shift @ARGV || $informative_sites_file_path;
  # Remove the trailing "/" if any
  if( defined( $output_dir ) ) {
    ( $output_dir ) = ( $output_dir =~ /^(.*[^\/])\/*$/ );
  }

  if( $VERBOSE ) { print "Output will be written in directory \"$output_dir\".."; }

  my ( $line );
  my %seqs;
  my @fields;
  my @seqorder;
  my @seq_as_list;
  my $seq;
  if( -e $informative_sites_file ) {
    if( $VERBOSE ) { print "Opening file \"$informative_sites_file\".."; }
    
    unless( open( INFORMATIVE_SITES_FILE_FH, $informative_sites_file ) ) {
      warn "Unable to open informative_sites_file file \"$informative_sites_file\": $!";
      return 1;
    }
    
    if( $VERBOSE ) { print ".done.\n"; }
    
    if( $VERBOSE ) { print "Reading informative sites file.."; }
    while( $line = <INFORMATIVE_SITES_FILE_FH> ) {
    
      #( $line ) = ( $line =~ /^(.+?)\s*$/ );  # Chop and Chomp won't remove ^Ms
      chomp( $line ); # Actually we don't want to remove trailing tabs...
      
      if( $DEBUG ) {
        print "LINE: $line\n";
      }
    
      next unless $line;
      next if $line =~ /^\s/;
      last if $line =~ /^Alignment/;
    
      if( $VERBOSE ) { print "."; }
    
      @fields = split( "\t", $line, -1 ); # The -1 means no limit, which forces inclusion of flanking empty fields.
    
      ## NOTE that there is a bug (as of September, 2015) in the inSites
      ## code online, in which flanking gaps are not printed in the
      ## output table.  THIS MUST BE FIXED WHEN READING/USING THE FILE.
      @fields = map { if( $_ =~ /^\s*$/ ) { '-' } else { $_ } } @fields;
      
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
  } # End if there is a file to read.

  # Maybe there are no informative sites.  We are just done, then.
  my $insites_fasta_file = "";
  if( !scalar( @seqorder ) ) {
    if( $VERBOSE ) {
      print( "No informative sites..\n" );
    }
    $force_one_cluster = 1;
    $insites_fasta_file = $original_fasta_filename;
  } else {
    # Ok, great.  Now, for all entries, replace '.' with the consensus value.
    unless( $seqorder[ 0 ] eq 'Consensus' ) {
      die( "First sequence is not 'Consensus' it is instead '$seqorder[ 0 ]'" );
    }
  
    shift @seqorder; # Remove `Consensus`
    
    if( $DEBUG ) {
      print "Sequences: ", join( ", ", @seqorder ), "\n";
    }
  
    my @consensus_as_list = @{ $seqs{ "Consensus" } };
    my $alignment_length = scalar( @consensus_as_list );
    delete $seqs{ "Consensus" };
  
    $insites_fasta_file =
      "${output_dir}/${informative_sites_file_short}.fasta";
    if( $VERBOSE ) { print "Opening file \"$insites_fasta_file\" for writing.."; }
    unless( open INSITES_FASTA_FH, ">$insites_fasta_file" ) {
      warn "Unable to open insites_fasta file \"$insites_fasta_file\": $!";
      return 1;
    }
    if( $VERBOSE ) { print ".done.\n"; }
  
    if( $VERBOSE ) {
      print "Printing alignment of informative sites..";
    }
    my $seq_wrapped;
    foreach my $seq_name ( @seqorder ) {
      if( $DEBUG ) {
        print $seq_name, "\n";
      }
      @seq_as_list = @{ $seqs{ $seq_name } };
  
      unless( scalar( @seq_as_list ) == $alignment_length ) {
        print( "Alignment length error in inSites file $insites_fasta_file: got " . scalar( @seq_as_list ) . ", expected $alignment_length" );
        print( ">Consensus\n", join( "", @consensus_as_list ), "\n" );
        print( ">$seq_name\n", join( "", @seq_as_list ), "\n" );
        die( "Alignment length error in inSites file $insites_fasta_file: got " . scalar( @seq_as_list ) . ", expected $alignment_length" );
      }
      for( my $i = 0; $i < scalar( @seq_as_list ); $i++ ) {
        if( $seq_as_list[ $i ] eq '.' ) {
          $seq_as_list[ $i ] = $consensus_as_list[ $i ];
        }
      }
      $seq = join "", @seq_as_list;
      
      # add newlines / wraparound to $seq
      $seq_wrapped = wrap( '', '',  $seq );
      print INSITES_FASTA_FH "> $seq_name\n$seq_wrapped\n";
    } # End foreac $seq_name
    
    if( $VERBOSE ) { print ".done.\n"; }
    
    if( $VERBOSE ) { print "Closing file \"$insites_fasta_file\".."; }
    close INSITES_FASTA_FH;
    if( $VERBOSE ) { print ".done.\n"; }
  } # End if there are no informative sites .. else ..

  if( $VERBOSE ) {
    print "Calling R to cluster informative sites..";
    if( $force_one_cluster ) {
      print( "Forcing one cluster..\n" );
    } else {
      print( "Clustering..\n" );
    }
  }
  my $R_output = `export clusterInformativeSites_forceOneCluster="$force_one_cluster"; export clusterInformativeSites_inputFilename="$insites_fasta_file"; export clusterInformativeSites_originalFastaFilename="$original_fasta_filename"; export clusterInformativeSites_outputDir="$output_dir"; R -f clusterInformativeSites.R --vanilla --slave`;
  print( $R_output );

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


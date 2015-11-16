#! /usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) translate_ambiguous
##  Author:
##      Paul T. Edlefsen   pedlefse@scharp.org
##  Description:
##      Script for translating DNA to AA, with "---" translated to "-"
##      and any other gap-combination translated to "X" unless there
##      is an unambiguous translation available (eg if the last
##      position is the gap char and the others uniquely identify an amino acid).
##
##
#******************************************************************************
#*  Copyright (C) 2015    Paul T. Edlefsen                                    *
#*  All rights reserved.                                                      *
#******************************************************************************

use Getopt::Std;
use Bio::Tools::CodonTable;
use Bio::Perl; # for revcom(..)

use strict;

use vars qw( $DEBUG $VERBOSE
             $opt_D $opt_V $opt_o $opt_O $opt_c );

sub translate_ambiguous {
  @ARGV = @_;

  sub translate_ambiguous_usage {
    print "\ttranslate_ambiguous [-DV] [-(o|O) <Output Fasta File>] <Input Fasta File> <frame|1> <reverse|0> <complement|0>\n";
    return 1;
  }

  # This means -D, -o, -O, -V are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_o is an optional filename to put the output.  Otherwise, STDOUT.
  # opt_O is just like opt_o except it'll overwrite the file if it exists.
  # opt_V means be verbose.
  # opt_c is the (1-indexed) column to insert newlines at (default 40)
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_o, $opt_O, $opt_c ) = ();
  if( not getopts('DVo:O:c:') ) {
    return translate_ambiguous_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;
  
  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }

  my $newline_column = $opt_c || 40;

  my $input_fasta_file = shift @ARGV;
  unless( $input_fasta_file ) { return translate_ambiguous_usage() };

  my $frame = shift @ARGV || 1;
  my $reverse = shift @ARGV || 0;
  my $complement = shift @ARGV || 0;

  my $output_file = $opt_o || $opt_O;
  if( defined( $output_file ) && ( $output_file ne '-' ) ) {
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

  if( $VERBOSE ) { print "Opening input Fasta file \"$input_fasta_file\".."; }
  open( INPUT_FASTA_FH, $input_fasta_file ) or
    die "Unable to open fasta file \"$input_fasta_file\": $!";
  if( $VERBOSE ) { print ".done.\n"; }

  if( $VERBOSE ) { print "Parsing Fasta file and writing output.."; }

  my $myCodonTable   = Bio::Tools::CodonTable->new();

  sub translate_ambiguous_codon {
    my $codon = shift;

    ## TODO: REMOVE
    #print "translate_ambiguous_codon( '$codon' )\n";
    unless( length( $codon ) ) { return '' };

    if( $codon =~ /^-+$/ ) {
      return '-';
    }

    # Here's the tricksy bit.  If it's not all-gaps, treat gaps as Ns
    # to get the most generous translation possible.
    $codon =~ s/-/N/g;

    return $myCodonTable->translate( $codon );
  } # translate_ambiguous_codon (..)

  sub translate_sequence {
    my $sequence = shift;
    my $frame = shift || 1;
    my $reverse = shift || 0;
    my $complement = shift || 0;

    if( $complement ) {
      my $seqobj    = Bio::PrimarySeq->new(-seq => $sequence );
      $sequence = revcom( $seqobj )->seq();
      unless( $reverse ) {
        $sequence = reverse( $sequence );
      }
    } elsif( $reverse ) {
      $sequence = reverse( $sequence );
    }

    if( $frame > 1 ) {
      if( $frame >3 ) {
        die( "\$frame is $frame, but it should be 1, 2, or 3." );
      } 
      $sequence = substr $sequence, ( $frame - 1 );
    } 

    #print "SEQUENCE:\n$sequence\n";
    my @sequence_codons = ( $sequence =~ /(\S\S?\S?)/g );
    #print "CODONS:\n", join( "\n", @sequence_codons ), "\n";

    #print "FIRST CODON, $sequence_codons[ 0 ], translates to ", translate_ambiguous_codon( $sequence_codons[ 0 ] ), "\n";

    my @translated_codons = map { translate_ambiguous_codon( $_ ) } @sequence_codons;

    return( join( '', @translated_codons ) );
  } # translate_sequence(..)

  sub print_sequence {
    my $sequence = shift;

    my $seq_len = length( $sequence );
    for( my $seq_pos = 0; $seq_pos < $seq_len; $seq_pos += $newline_column ) {
      if( ( $seq_pos + $newline_column ) > $seq_len ) {
        print OUTPUT_FH substr( $sequence, $seq_pos, ( $seq_len - $seq_pos ) ), "\n\n";
      } else {
        print OUTPUT_FH substr( $sequence, $seq_pos, $newline_column ), "\n";
      }
    }
  } # print_sequence (..)

  my $sequence;
  my $sequence_header;
  my $have_a_match = 0;
  while( <INPUT_FASTA_FH> ) {
  
    if( /^>/ ) {
      if( defined( $sequence ) ) {

        # translate it first.
        $sequence = translate_sequence( $sequence, $frame, $reverse, $complement );

        print OUTPUT_FH $sequence_header;
        print_sequence( $sequence );
        undef $sequence;
      }
      $sequence_header = $_;
      $have_a_match = 1;
      if( $DEBUG ) {
        print "Found sequence header: $sequence_header";
      }
    } elsif( $have_a_match ) {
      # Get rid of all preexisting whitespace
      $_ =~ s/\s//g;
  
      $sequence .= $_;
    }
  }
  if( defined( $sequence ) ) {

    # translate it first.
    $sequence = translate_sequence( $sequence, $frame, $reverse, $complement );

    print OUTPUT_FH $sequence_header;
    print_sequence( $sequence );
    undef $sequence;
  }
  if( $VERBOSE ) { print ".done.\n"; }


  if( $VERBOSE  ) { print "Closing file \"$input_fasta_file\".."; }
  close INPUT_FASTA_FH;
  if( $VERBOSE ) { print ".done.\n"; }

  if( $VERBOSE && $output_file ) { print "Closing file \"$output_file\".."; }
  close OUTPUT_FH;
  if( $VERBOSE && $output_file ) { print ".done.\n"; }

  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  return 0;
} # translate_ambiguous(..)

translate_ambiguous( @ARGV );

1;


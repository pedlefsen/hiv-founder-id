#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) runGeneCutterOnline
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##  Description:
##      Script for running GeneCutter using the online LANL tools to determine which
##      sequences are recombinants of other sequences in a set.
##      
###******************************************************************************

use Getopt::Std; # for getopts
use File::Path qw( make_path );
use Path::Tiny;
use Try::Tiny;

# For screenscraping
use WWW::Mechanize;
use LWP::Simple;

use strict;
use vars qw( $opt_D $opt_V $opt_e $opt_p );
use vars qw( $VERBOSE $DEBUG );

sub runGeneCutterOnline {
  @ARGV = @_;

  my $DEFAULT_EMAIL_ADDRESS = "pedlefsen\@gmail.com";

  sub runGeneCutterOnline_usage {
    print "\trunGeneCutterOnline [-DVp] [-e <youremail\@ddr.ess>] <input_fasta_file> [<output_dir>]\n";
    exit;
  }

  # This means -D and -V are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # opt_p means that the input files are pre-aligned.
  # opt_e means use a different email address than the default.
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_e, $opt_p ) = ();
  if( not getopts('DVep') ) {
    runGeneCutterOnline_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;

  my $is_prealigned = $opt_p || 0;
  my $emailAddress = $opt_e || $DEFAULT_EMAIL_ADDRESS;
  
  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }
  
  my $input_fasta_file = shift @ARGV || runGeneCutterOnline_usage();
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

  my $fasta_file_readyForGeneCutter = $input_fasta_file;
  
  my $mech = WWW::Mechanize->new( autocheck => 1 );
  $mech->get( "http://www.hiv.lanl.gov/content/sequence/GENE_CUTTER/cutter.html" );

  my $result = $mech->submit_form(  form_number => 2,
                                  fields    => {
                                                ORGANISM => "HIV-1",
                                                UPLOAD => $fasta_file_readyForGeneCutter,
                                                PREALIGNED => ( $is_prealigned ? "YES" : "NO" ),
                                                SEG => "ALL",
                                                INSERTSTDSEQ => "YES",
                                                REMOVESTDSEQ => "NO",
                                                ALIGN => "YES", # TODO: TRY NO? THIS IS FOR CODON-ALIGNING THE REGION; I THINK IT'S JUST TWEAKING IT.
                                                PROTEIN => "YES3", # THIS TOLERATES SOME AMBIGUITY I GUESS
                                                ACTION => "DOWNPC",
                                                VIEW => "YES"
                                               }
                       );
  
  my $content = $mech->content();
  if( $DEBUG ) {
    print "OK\n"#\$content is $content\n";
  }
  $mech->submit_form( form_number => 2,
                                  fields    => {
                                                titleFromUser => $input_fasta_file_short_nosuffix,
                                                EMAIL => $emailAddress,
                                                EMAIL2 => $emailAddress
                                               }
 );
  my $content2 = $mech->content();
  if( $DEBUG ) {
    print "OK2\n";# \$content2 is $content2\n";
  }
  
  my ( $jobTitle, $jobID ) = ( $content2 =~ /Your job has been submitted .+The job title is \<b\>(.+)\<\/b\> and your reference number is \<b\>(\d+)\<\/b\>/ );
  unless( defined $jobTitle ) {
    die( "Error running GeneCutter online: GOT $content2" );
  }
  if( $VERBOSE ) {
    print "JOB ID IS $jobID\n";
  }
  ## Results page
  my $results_url = "http://www.hiv.lanl.gov/tmp/download/GENE_CUTTER/$jobID/FRAMESET_PRO.html";
  ## Actual results are here:
  my $pro_prepresults_url = "http://www.hiv.lanl.gov/tmp/download/GENE_CUTTER/$jobID/printpro.html";
  my $pro_results_url = "http://www.hiv.lanl.gov/tmp/download/GENE_CUTTER/$jobID/CONTROLS_PRO.html";
  my $nuc_prepresults_url = "http://www.hiv.lanl.gov/tmp/download/GENE_CUTTER/$jobID/printnucs.html";
  my $nuc_results_url = "http://www.hiv.lanl.gov/tmp/download/GENE_CUTTER/$jobID/CONTROLS_NUCS.html";
  my $resultsContent = undef;
  try {
    $mech->get( $results_url );
    $resultsContent = $mech->content();
    if( $DEBUG ) {
      print( "got $resultsContent" );
    }
  } catch {
    ## DO NOTHING.
    #warn "caught error: $_";
  };
  while( !defined( $resultsContent ) || ( $resultsContent =~ /No such file/ ) ) {
    if( $VERBOSE || $DEBUG ) {
      print( "trying again to get the gene cutter results in 1 second..\n" );
    }
      sleep( 1 );
      try {
        $mech->get( $results_url );
        $resultsContent = $mech->content();
      } catch {
        # do nothing.
      }
    }
  if( $DEBUG ) {
    print( "finally, got $resultsContent" );
  }

  $mech->get( $nuc_prepresults_url );
  $mech->get( $nuc_results_url );
  #my $nuc_results_content = $mech->content();
  $mech->submit_form( 
                                  fields    => {
                                                REGION => "ALL",
                                                OUTFORMAT => "FASTA",
                                                PROTEINLIST => "-GAG-POL-VIF-VPR-TAT-REV-VPU-ENV-NEF-FULL_SEQUENCE",
                                                PROTEIN_FLAG => "NO",
                                                DIR => "/tmp/download/GENE_CUTTER/$jobID"
                                               }
 );
  $mech->save_content( "${output_path_dir}/${input_fasta_file_short_nosuffix}_allnucs.zip" );
  #my $content3 = $mech->content();
  if( $DEBUG ) {
    print "OK3\n";# \$content3 is $content3\n";
  }

  $mech->get( $pro_prepresults_url );
  $mech->get( $pro_results_url );
  my $pro_results_content = $mech->content();

  $mech->submit_form( 
                                  fields    => {
                                                REGION => "ALL",
                                                OUTFORMAT => "FASTA",
                                                PROTEINLIST => "-GAG-POL-VIF-VPR-TAT-REV-VPU-ENV-NEF-FULL_SEQUENCE",
                                                PROTEIN_FLAG => "YES",
                                                DIR => "/tmp/download/GENE_CUTTER/$jobID"
                                               }
 );
  $mech->save_content( "${output_path_dir}/${input_fasta_file_short_nosuffix}_allproteins.zip" );
  #my $content4 = $mech->content();
  if( $DEBUG ) {
    print "OK4\n";# \$content4 is $content4\n";
  }

  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }

  return 0;
} # runGeneCutterOnline(..)

runGeneCutterOnline( @ARGV );

1;


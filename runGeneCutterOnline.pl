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
use File::Temp qw / :POSIX /;
use Path::Tiny;
use Try::Tiny;

# For screenscraping
use WWW::Mechanize;
use LWP::Simple;

use strict;
use vars qw( $opt_D $opt_V $opt_e $opt_p $opt_P $opt_R $opt_x );
use vars qw( $VERBOSE $DEBUG );

sub runGeneCutterOnline {
  @ARGV = @_;

  my $DEFAULT_EMAIL_ADDRESS = "pedlefsen\@gmail.com";

  sub runGeneCutterOnline_usage {
    print "\trunGeneCutterOnline [-DVp] [-P <proteins_list>] [-R <region>] [-e <youremail\@ddr.ess>] [-x <output_zipfile_prefix>] <input_fasta_file> [<output_dir>]\n";
    exit;
  }

  # This means -D and -V, -p, -e, -P, -R, -x are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # opt_p means that the input files are pre-aligned.
  # opt_e means use a different email address than the default.
  # opt_P changes the proteins to evaluate string (default: "-GAG-POL-VIF-VPR-TAT-REV-VPU-ENV-NEF")a
  # opt_R changes the region of the genome that is being evaluated (default: "ALL")
  # opt_x is an optional prefix to add before the "_allnucs.zip" and "_allproteins.zip"
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_p, $opt_e, $opt_P, $opt_R, $opt_x ) = ();
  if( not getopts('DVpe:P:R:x:') ) {
    runGeneCutterOnline_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;

  my $is_prealigned = $opt_p || 0;
  my $emailAddress = $opt_e || $DEFAULT_EMAIL_ADDRESS;
  my $proteins_to_evaluate = $opt_P || "-GAG-POL-VIF-VPR-TAT-REV-VPU-ENV-NEF";
  # NOTE: -FULL_SEQUENCE seems to not ever give work without giving the in-fasta-file-output message eg "Error: I can't open file /tmp/download/GENE_CUTTER/30488/FULL_SEQUENCE.NA.RAW.: No such file or directory"
  my $region = $opt_R || "ALL";
  my $output_zip_prefix = $opt_x || "";
  if( length( $output_zip_prefix ) > 0 ) {
    unless( $output_zip_prefix =~ /^_/ ) {
      $output_zip_prefix = '_' . $output_zip_prefix;
    }
  }
#  if( $proteins_to_evaluate !~ /^-/ ) {
#    stop( "The -P argument must be formatted as '-GENE1-GENE2-GENE3' eg '-GAG-POL-VIF-VPR-TAT-REV-VPU-ENV-NEF'" );
#  }
  
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
  
  # Genecutter now rejects filenames containing dots apparently, or something (ie even in the directory name), so we use a temporary symlink in the current dir and then strip the dots.
  $fasta_file_readyForGeneCutter = File::Temp::tempnam( ".", "tmp" ) . ".fasta"; # this doesn't actually create the file.
  `cp $input_fasta_file $fasta_file_readyForGeneCutter`;
  # Ok and now annoyingly we have to strip the "./" because the latest version of genecutter doesn't like it.
  ( $fasta_file_readyForGeneCutter ) = ( $fasta_file_readyForGeneCutter =~ /^\.\/(.+)$/ );
  if( $VERBOSE ) {
    print "Created temporary copy of the input fasta file at $fasta_file_readyForGeneCutter\n";
  }

  my $mech = WWW::Mechanize->new( autocheck => 1 );
  $mech->get( "http://www.hiv.lanl.gov/content/sequence/GENE_CUTTER/cutter.html" );

  my                                   %fields    = (
                                                ORGANISM => "HIV-1",
                                                UPLOAD => $fasta_file_readyForGeneCutter,
                                                PREALIGNED => ( $is_prealigned ? "YES" : "NO" ),
                                                SEG => $region,
                                                INSERTSTDSEQ => "YES",
                                                REMOVESTDSEQ => "NO",
                                                ALIGN => "YES", # THIS IS FOR CODON-ALIGNING THE REGION; I THINK IT'S JUST TWEAKING IT IF IT'S PREALIGNED.
                                                PROTEIN => "YES3", # THIS TOLERATES SOME AMBIGUITY I GUESS
                                                ACTION => "DOWNPC",
                                                VIEW => "YES"
                                               );
  if( $DEBUG ) {
    foreach my $key ( keys %fields ) {
      print $key .  " => " . $fields{ $key } . "\n";
    }
  }

  my $result = $mech->submit_form(  form_number => 2,
                                  fields    => \%fields
                       );
  
  my $content = $mech->content();
  if( $DEBUG ) {
    print "OK\n";#\$content is $content\n";
  }
  $mech->submit_form( form_number => 2,
                                  fields    => {
                                                titleFromUser => "notitle", #$input_fasta_file_short_nosuffix
                                                EMAIL => $emailAddress,
                                                EMAIL2 => $emailAddress
                                               }
 );
  my $content2 = $mech->content();
  if( $DEBUG ) {
    print "OK2\n";# \$content2 is $content2\n";
  }
  
  if( !defined( $content2 ) || ( $content2 =~ /Request Rejected/ ) ) {
    ## Update: one reason this was happening (at least) was that the run name was too long or something, so see above where I changed the name to "notitle".
    die( "RAN INTO THE DREADED BUG AT GENECUTTER THAT WE DON'T UNDERSTAND.\n" );
  }
  if( $DEBUG ) {
    print( "got $content2" );
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
                                                REGION => $region,
                                                OUTFORMAT => "FASTA",
                                                PROTEINLIST => $proteins_to_evaluate,
                                                PROTEIN_FLAG => "NO",
                                                DIR => "/tmp/download/GENE_CUTTER/$jobID"
                                               }
 );
  $mech->save_content( "${output_path_dir}/${input_fasta_file_short_nosuffix}${output_zip_prefix}_allnucs.zip" );
  #my $content3 = $mech->content();
  if( $DEBUG ) {
    print "OK3\n";# \$content3 is $content3\n";
  }

  $mech->get( $pro_prepresults_url );
  $mech->get( $pro_results_url );
  my $pro_results_content = $mech->content();

  $mech->submit_form( 
                                  fields    => {
                                                REGION => $region,
                                                OUTFORMAT => "FASTA",
                                                PROTEINLIST => $proteins_to_evaluate,
                                                PROTEIN_FLAG => "YES",
                                                DIR => "/tmp/download/GENE_CUTTER/$jobID"
                                               }
 );
  $mech->save_content( "${output_path_dir}/${input_fasta_file_short_nosuffix}${output_zip_prefix}_allproteins.zip" );
  #my $content4 = $mech->content();
  if( $DEBUG ) {
    print "OK4\n";# \$content4 is $content4\n";
  }

  ## See above where we created a tmpfile for this.
  `rm $fasta_file_readyForGeneCutter`;
  if( $VERBOSE ) {
    print "Removed temporary file $fasta_file_readyForGeneCutter\n";
  }

  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }

  return 0;
} # runGeneCutterOnline(..)

runGeneCutterOnline( @ARGV );

1;


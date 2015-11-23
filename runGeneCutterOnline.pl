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
#use HTTP::Request::Common qw(POST);
#use LWP::UserAgent;

use strict;
use vars qw( $opt_D $opt_V $opt_e );
use vars qw( $VERBOSE $DEBUG );

sub runGeneCutterOnline {
  @ARGV = @_;

  my $DEFAULT_EMAIL_ADDRESS = "pedlefsen\@gmail.com";

  sub runGeneCutterOnline_usage {
    print "\trunGeneCutterOnline [-DV] [-e <youremail\@ddr.ess>] <input_fasta_file> [<output_dir>]\n";
    exit;
  }

  # This means -D and -V are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # opt_e means use a different email address than the default.
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_e ) = ();
  if( not getopts('DVe') ) {
    runGeneCutterOnline_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;

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

  ## Ok, for some reason GeneCutter has strict rules about header lines...
  ## Reformat the file.
  my $formattedForGeneCutter = 0;
    if( 0 ) {
      my $input_fasta_file_contents = path( $input_fasta_file )->slurp();
      if( $input_fasta_file_contents =~ /[\|\/]/ ) {
        if( $VERBOSE ) {
          print( "Input file \"$input_fasta_file\" contains illegal characters \"|\"; changing them to \"-x-BAR-x-\" or \"-x-SLASH-x-\"\n" );
          ## TODOL REMOVE
          #print( $output_path_dir );
        }
        $formattedForGeneCutter = 1;
        $input_fasta_file_contents =~ s/\|/-x-BAR-x-/g;
        $input_fasta_file_contents =~ s/\//-x-SLASH-x-/g;
        $input_fasta_file_contents =~ s/\\/-x-BACKSLASH-x-/g;
        $input_fasta_file_contents =~ s/\.\./-x-DOTDOT-x-/g;
        $input_fasta_file_contents =~ s/\./-x-DOT-x-/g;
        # Now write it out to a temporary location in the output dir.
        $input_fasta_file_path = $output_path_dir;
        $input_fasta_file_short_nosuffix = "${input_fasta_file_short_nosuffix}_formattedForGeneCutter";
        $input_fasta_file_short = "${input_fasta_file_short_nosuffix}${input_fasta_file_suffix}";
        $input_fasta_file = "${input_fasta_file_path}/${input_fasta_file_short}";
        
        if( $VERBOSE ) {
          print( "Writing out fixed input file \"$input_fasta_file\".." );
        }
        if( $VERBOSE ) { print "Opening file \"$input_fasta_file\" for writing..\n"; }
        unless( open input_fasta_fileFH, ">$input_fasta_file" ) {
            warn "Unable to open output file \"$input_fasta_file\": $!\n";
            return 1;
          }
        print input_fasta_fileFH $input_fasta_file_contents;
        close( input_fasta_fileFH );
      }
    }  

  my $fasta_file_readyForGeneCutter = $input_fasta_file;
  
  my $R_output;

  # GeneCutter seems to not be able to handle duplicate sequences.
  ## STEP 1: remove duplicate sequences.
#   if( $VERBOSE ) {
#     print "Calling R to create a version of the fasta file in which duplicate sequences are removed..";
#   }
#   $R_output = `export removeDuplicateSequencesFromAlignedFasta_inputFilename="$input_fasta_file"; export removeDuplicateSequencesFromAlignedFasta_outputDir="$output_path_dir"; R -f removeDuplicateSequencesFromAlignedFasta.R --vanilla --slave`;
#   # The output has the file name of the consensus file.
#   if( $DEBUG ) {
#     print( "GOT: $R_output\n" );
#   }
#   if( $VERBOSE ) {
#     print ".done.\n";
#   }
#   # Parse it to get the output fasta filename.
#   my ( $fasta_file_no_duplicates ) = ( $R_output =~ /\"([^\"]+)\"/ );
#   print "Fasta file with duplicates removed: $fasta_file_no_duplicates\n";
# 
#   my ( $table_file_no_duplicates ) =
#     ( $fasta_file_no_duplicates =~ /^(.+)$input_fasta_file_suffix$/ );
#   $table_file_no_duplicates .= ".tbl";
#   
#   if( -e $table_file_no_duplicates ) {
#     print "Table of duplicates removed: $table_file_no_duplicates\n";
#   }
# 
#   # GeneCutter has a problem with certain characters in fasta headers. 
#   ## STEP 2: Rename the seqs to just their numbers.
#   if( $VERBOSE ) {
#     print "Calling R to create a version of the fasta file in which seqs are numbered instead of named..";
#   }
#   my ( $fasta_file_no_duplicates_numbered ) =
#     ( $fasta_file_no_duplicates =~ /^(.+)$input_fasta_file_suffix$/ );
#   $fasta_file_no_duplicates_numbered .= "_numbered$input_fasta_file_suffix";
#   $R_output = `export computeConsensusSequenceFromAlignedFasta_inputFilename="$fasta_file_no_duplicates"; export computeConsensusSequenceFromAlignedFasta_outputFilename="$fasta_file_no_duplicates_numbered"; export computeConsensusSequenceFromAlignedFasta_includeFullAlignment="TRUE"; export  computeConsensusSequenceFromAlignedFasta_includeConsensus="FALSE"; export computeConsensusSequenceFromAlignedFasta_useSeqeunceNumbersAsNames="TRUE"; R -f computeConsensusSequenceFromAlignedFasta.R --vanilla --slave`;
#   # The output has the file name of the "consensus file" which is not in fact the consensus but the fasta with renamed seqs.
#   if( $DEBUG ) {
#     print( "GOT: $R_output\n" );
#   }
#   if( $VERBOSE ) {
#     print ".done.\n";
#   }
#   # Parse it to get the filename.
#   my ( $fasta_file_readyForGeneCutter ) = ( $R_output =~ /\"([^\"]+)\"/ );
#   if( $DEBUG ) {
#     print "Fasta file ready for GeneCutter: $fasta_file_readyForGeneCutter\n";
#   }

  ## Set up the table mapping sequence names to numbers.
#       if( $VERBOSE ) {
#         print "Reading sequence names from file \"", $fasta_file_no_duplicates, "\"..";
#       }
#       my $fasta_file_no_duplicates_contents =
#            path( $fasta_file_no_duplicates )->slurp();
#       if( $VERBOSE ) {
#         print ".done\n";
#       }
#       if( $DEBUG ) {
#         #print $fasta_file_no_duplicates_contents;
#       }
#       my ( @seq_names ) =  ( $fasta_file_no_duplicates_contents =~ /\n?>[ \t]*(.+) *\n/g );

  my $mech = WWW::Mechanize->new( autocheck => 1 );
  $mech->get( "http://www.hiv.lanl.gov/content/sequence/GENE_CUTTER/cutter.html" );

  my $result = $mech->submit_form(  form_number => 2,
                                  fields    => {
                                                ORGANISM => "HIV-1",
                                                UPLOAD => $fasta_file_readyForGeneCutter,
                                                PREALIGNED => "YES",
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
  ## TODO: REMOVE?
  print "JOB ID IS $jobID\n";
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
    print( "got $resultsContent" );
  } catch {
    ## DO NOTHING.
    #warn "caught error: $_";
  };
    while( !defined( $resultsContent ) || ( $resultsContent =~ /No such file/ ) ) {
      print( "trying again to get the gene cutter results.\n" );
      sleep( 1 );
      try {
        $mech->get( $results_url );
        $resultsContent = $mech->content();
      } catch {
        # do nothing.
      }
    }
    print( "finally, got $resultsContent" );

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
  $mech->save_content( "allnucs.zip" );
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
  $mech->save_content( "allproteins.zip" );
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


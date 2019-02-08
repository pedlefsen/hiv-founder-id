#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) runGeneCutterOnline
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##  Description:
##      Script for running GeneCutter using the online LANL tools to split
##      nflg or half-genome sequence alignments into gene-specific subalignments.
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
use vars qw( $opt_D $opt_V $opt_e $opt_p $opt_P $opt_R $opt_x $opt_n );
use vars qw( $VERBOSE $DEBUG );

sub runGeneCutterOnline {
  @ARGV = @_;

  my $DEFAULT_EMAIL_ADDRESS = "pedlefsen\@gmail.com";

  sub runGeneCutterOnline_usage {
    print "\trunGeneCutterOnline [-DVp] [-P <proteins_list>] [-R <region>] [-e <youremail\@ddr.ess>] [-x <output_zipfile_prefix>] <input_fasta_file> [<output_dir>]\n";
    exit;
  }

  # This means -D and -V, -p, -e, -P, -R, -x, -n are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # opt_p means that the input files are pre-aligned.
  # opt_e means use a different email address than the default.
  # opt_P changes the proteins to evaluate string (default: "-GAG-POL-VIF-VPR-TAT-REV-VPU-ENV-NEF")a
  # opt_R changes the region of the genome that is being evaluated (default: "ALL")
  # opt_x is an optional prefix to add before the "_allnucs.zip" and "_allproteins.zip"
  # opt_n means get nucleotides only, not also translated into proteins.
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_p, $opt_e, $opt_P, $opt_R, $opt_x, $opt_n ) = ();
  if( not getopts('DVpe:P:R:x:n') ) {
    runGeneCutterOnline_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;

  my $is_unaligned = !$opt_p;
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
  my $get_amino_acids_too = !$opt_n;
  
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

  if ($DEBUG){
    print "\nChecking the input file names and output path:";
    print "\n\$input_fasta_file = $input_fasta_file";
    print "\n\$input_fasta_file_path = $input_fasta_file_path";
    print "\n\$input_fasta_file_short = $input_fasta_file_short";
    print "\n\$output_path_dir = $output_path_dir\n";
  }

  # GeneCutter has a problem with certain characters in fasta headers. 
  ## Step 1: Set up the table mapping sequence names to numbers.
      if( $VERBOSE ) {
        print "Reading sequence names from file \"", $input_fasta_file, "\"..";
      }
      my $input_fasta_file_contents =
           path( $input_fasta_file )->slurp();
      if( $VERBOSE ) {
        print ".done\n";
      }
      if( $DEBUG ) {
        #print $input_fasta_file_contents;
      }
      my ( @seq_names ) =  ( $input_fasta_file_contents =~ /\n?>[ \t]*(.+) *\n/g );

  ## STEP 2: Rename the seqs to just their numbers.
  my ( $input_fasta_file_numbered ) =
    ( $input_fasta_file =~ /^(.+)$input_fasta_file_suffix$/ );
  $input_fasta_file_numbered .= "_numbered$input_fasta_file_suffix";

  if( !( -e $input_fasta_file_numbered ) ) {
    if( $VERBOSE ) {
      print "Creating a version of the fasta file in which seqs are numbered instead of named ($input_fasta_file_numbered)..";
    }
    if( $VERBOSE ) { print "Opening file \"$input_fasta_file\".."; }
    unless( open( INPUT_FASTA_FH, $input_fasta_file ) ) {
      warn "Unable to open input_fasta file \"$input_fasta_file\": $!";
      return 1;
    }
    if( $VERBOSE ) { print ".done.\n"; }
    
    if( $VERBOSE ) { print "Opening file \"$input_fasta_file_numbered\" for writing.."; }
    unless( open OUTPUT_FH, ">$input_fasta_file_numbered" ) {
      warn "Unable to open output file \"$input_fasta_file_numbered\": $!";
      return 1;
    }
    if( $VERBOSE ) { print ".done.\n"; }
  
    for( my $counter = 1; <INPUT_FASTA_FH>; ) {
      if( $_ =~ /^>/ ) { print OUTPUT_FH ">$counter\n"; $counter = $counter + 1 } else { print OUTPUT_FH $_ }
    }
  
    if( $VERBOSE ) { print ".done.\n"; }
    
    if( $VERBOSE && $input_fasta_file_numbered ) { print "Closing file \"$input_fasta_file_numbered\".."; }
    close OUTPUT_FH;
    if( $VERBOSE && $input_fasta_file_numbered ) { print ".done.\n"; }
    if( $VERBOSE ) { print "Closing file \"$input_fasta_file\".."; }
    close INPUT_FASTA_FH;
    if( $VERBOSE ) { print ".done.\n"; }
  } # End unless it already exists, create the numbered version of the input fasta file.

  # Also fix the filename.
  my $fasta_file_readyForLANL = File::Temp::tempnam( ".", "tmp" ) . ".fasta";
  # remove leading dots and slashes from filename:
  ( $fasta_file_readyForLANL ) = ( $fasta_file_readyForLANL =~ /^\.\/(.+)$/ );
  # remove any underscores:
  $fasta_file_readyForLANL =~ s/\_//;

  `cp $input_fasta_file_numbered $fasta_file_readyForLANL`;

  if( $VERBOSE ) {
    print "Created temporary copy of the input fasta file at $fasta_file_readyForLANL\n";
  }

  my $mech = WWW::Mechanize->new( autocheck => 1 );
  
  $mech->get( "http://www.hiv.lanl.gov/content/sequence/GENE_CUTTER/cutter.html" );

  if( $DEBUG ){
    print "\n===========================================================\n";
    print "Initial Form requested:\n";
    #    print $mech->content();
    print "\nEnd of initial form\n";
    print "============================================================\n";
  }

  my                                   %fields    = (
                                                ORGANISM => "HIV-1",
                                                UPLOAD => $fasta_file_readyForLANL,
                                                PREALIGNED => ( $is_unaligned ? "YES" : "NO" ), # "PREALIGNED" actually means "unaligned" from page text
                                                SEG => $region,
                                                INSERTSTDSEQ => "YES",
                                                REMOVESTDSEQ => "NO",
                                                     ALIGN => ( $is_unaligned ? "YES" : "NO" ), #"YES", # THIS IS FOR CODON-ALIGNING THE REGION; I THINK IT'S JUST TWEAKING IT IF IT'S ACTUALLY prealigned (but note "PREALIGNED" in this form actually means unaligned).
                                                compensatingNum => "15",
                                                PROTEIN => ( $get_amino_acids_too ? "YES3" : "NO" ),
                                                ACTION => "DOWNPC",
                                                VIEW => "YES"
                                               );
  if( $DEBUG ) {
    print "\n===========================================================\n";
    print "Content that will be submitted for the intial form\n";
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
  my $results_url = "http://www.hiv.lanl.gov/tmp/download/GENE_CUTTER/$jobID/FRAMESET_NUCS.html"; # PL changed this from .../FRAMESET_PRO.html
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
  for( my $num_waits = 1; !defined( $resultsContent ) || ( $resultsContent eq "" ) || ( $resultsContent =~ /No such file/ ); $num_waits++ ) {
    if( $VERBOSE || $DEBUG ) {
      print( "trying again to get the gene cutter results in $num_waits seconds..\n" );
    }
      sleep( $num_waits );
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
  $mech->submit_form( form_number => 1,
                                  fields    => {
                                                REGION => $region,
                                                OUTFORMAT => "FASTA",
                                                PROTEINLIST => $proteins_to_evaluate,
                                                PROTEIN_FLAG => "NO",
                                                DIR => "/tmp/download/GENE_CUTTER/$jobID"
                                               }
                    );

  ## TODO: ERE I AM. The problem is that these will be numbered; we need to expand it, replace numbers with names, and re-collapse it!
  $mech->save_content( "${output_path_dir}/${input_fasta_file_short_nosuffix}${output_zip_prefix}_allnucs.zip" );
  #my $content3 = $mech->content();
  if( $DEBUG ) {
    print "OK3\n";# \$content3 is $content3\n";
  }

  if( $get_amino_acids_too ) {
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
  } # End if $get_amino_acids_too

  ## See above where we created a tmpfile for this.
  `rm $fasta_file_readyForLANL`;
  if( $VERBOSE ) {
    print "Removed temporary file $fasta_file_readyForLANL\n";
  }

  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }

  return 0;
} # runGeneCutterOnline(..)

runGeneCutterOnline( @ARGV );

1;


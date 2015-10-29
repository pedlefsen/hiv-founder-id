#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) runGARDOnline
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##  Description:
##      Script for running GARD using the online datamonkey tools (at
##      the www.datamonkey.org web site). Note that this creates output files in
##      subdirectories named after the input fasta file name.
##      
###******************************************************************************

use Getopt::Std; # for getopts
use File::Path qw( make_path );
use Path::Tiny;

# For screenscraping
use WWW::Mechanize;
#use HTTP::Request::Common qw(POST);
#use LWP::UserAgent;

use strict;
use vars qw( $opt_D $opt_V );
use vars qw( $VERBOSE $DEBUG );

sub runGARDOnline {
  @ARGV = @_;

  sub runGARDOnline_usage {
    print "\trunGARDOnline [-DV] <input_fasta_file> [<output_dir>]\n";
    exit;
  }

  # This means -D and -V are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # But first reset the opt vars.
  ( $opt_D, $opt_V ) = ();
  if( not getopts('DV') ) {
    runGARDOnline_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;
  
  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }
  
  my $input_fasta_file = shift @ARGV || runGARDOnline_usage();
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
  $mech->get( "http://www.datamonkey.org/dataupload.php" );
  #print( $mech->content() );
  my $result = $mech->submit_form(
                                  form_name => 'uploadform',
                                  fields    => { datatype => '0', code => '0', upfile => $input_fasta_file }
                       );
  
  # Have to submit it again, ie my $result2 = $result->submit();
  my $content = $mech->content();
  if( $DEBUG ) {
    print "OK1\n \$content is $content\n";
  }
  my ( $datamonkey_filename ) = ( $content =~ m/name='filename' value=\'([^\']+)\'/ );
#   <FORM method='POST' enctype='multipart/form-data' action='http://www.datamonkey.org/cgi-bin/datamonkey/finishUpload.pl'>
# <DIV style = 'margin:10px; text-align:center;'><input type='Hidden' name='filename' value='upload.21231639675385.1'><input type='Hidden' name='genCodeID' value='0'>
# <input type='Submit' value='Proceed to the analysis menu' style = 'background-color:purple; color:white; font-size:18px;'> <!--<a href='http://www.datamonkey.org/help/tree.php' target = '_blank' class = 'INFO'>Help</a>--> </DIV></FORM>	
  my $result2 = $mech->submit_form(
                                   form_number => 2,
                                   fields => {
                                     genCodeID => '0',
                                     filename => $datamonkey_filename
                                             }

  );
  $content = $result2->content();

  if( $DEBUG ) {
    print "OK2\n \$content is $content\n";
  }

  my ( $datamonkey_sequences ) = ( $content =~ m/name='sequences' value=\'([^\']+)\'/ );
  my ( $datamonkey_sites ) = ( $content =~ m/name='sites' value=\'([^\']+)\'/ );
  my ( $datamonkey_partitions ) = ( $content =~ m/name='partitions' value=\'([^\']+)\'/ );
  my ( $datamonkey_modelString ) = ( $content =~ m/NAME=\"modelString\" VALUE=\"([^\"]+)\"/ );
  my $result3 = $mech->submit_form(
                                   form_name => 'modelForm',
                                   fields => {
                                              genCodeID => '0',
                                              filename => $datamonkey_filename,
                                              sequences => $datamonkey_sequences,
                                              sites => $datamonkey_sites,
                                              partitions => $datamonkey_partitions,
                                              method => 21, # 21 is for 'GARD' (20 is SBP and 22 is for ASR)
                                              modelString => "012345", # for 'REV'
                                              NamedModels => "012345", # for 'REV'
                                              rateOption => '2', # for "Beta-Gamma",
                                              rateClasses => '4', # for 4 rate classes
                                              dNdS => '1.0', # DEFAULT
                                              rOptions => '4', # 'Estimated' global dN/dS value (the default)
                                              ambChoice => '0', # DEFAULT (Handling ambiguities: "Averaged")
                                              prime_property_choice => '0', # DEFAULT
                                              sigLevel => '0.1', # DEFAULT
                                              AC => '0', # DEFAULT
                                              AT => '0', # DEFAULT
                                              CG => '0', # DEFAULT
                                              CT => '1', # DEFAULT
                                              GT => '0' # DEFAULT
                                             }

  );
  $content = $result3->content();
  if( $DEBUG ) {
    print "OK3\n \$content is $content\n";
  }
  ## Parse out the refresh url
  my ( $datamonkey_gard_joburl ) = ( $content =~ m/http-equiv=\"refresh\" content=\"0;url=([^\"]+)\"/ );
  print( "JOB URL: $datamonkey_gard_joburl" );
  $mech->get( "http://www.datamonkey.org${datamonkey_gard_joburl}" );
  $content = $mech->content();
  if( $DEBUG ) {
    print "OK4\n \$content is $content\n";
  }

  ### ERE I AM, not done.  The problem is debugging this while the server queue is full.  Have to wait and do this another time.
  exit( 1 );

  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  return 0;
} # runGARDOnline(..)

runGARDOnline( @ARGV );

1;


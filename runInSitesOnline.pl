#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) runInSitesOnline
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##  Description:
##      Script for running InSites using the online DiveIn tools (at
##      the Mullins lab web site) to get the input needed to calculate
##      the InSites informative to privite sites ratio statistic: see
##      getInSitesStat.  Note that this creates output files in
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

sub runInSitesOnline {
  @ARGV = @_;

  sub runInSitesOnline_usage {
    print "\trunInSitesOnline [-DV] <input_fasta_file> [<output_dir>]\n";
    exit;
  }

  # This means -D and -V are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # But first reset the opt vars.
  ( $opt_D, $opt_V ) = ();
  if( not getopts('DV') ) {
    runInSitesOnline_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;
  
  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }
  
  my $input_fasta_file = shift @ARGV || runInSitesOnline_usage();
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
    ( $output_path_dir ) = ( $output_path_dir =~ /^(.+)\/*$/ );
  }
  make_path( $output_path_dir );

    ## HACK: make sure there are no bangs in the input file (since there are, right now).
    if( 1 ) {
      my $input_fasta_file_contents = path( $input_fasta_file )->slurp();
      if( $input_fasta_file_contents =~ /\!/ ) {
        if( $VERBOSE ) {
          print( "Input file \"$input_fasta_file\" contains illegal characters \"!\"; changing them to gaps\n" );
        }
        $input_fasta_file_contents =~ s/\!/-/g;
        # Now write it out to a temporary location in the output dir.
        $input_fasta_file_path = $output_path_dir;
        $input_fasta_file = "$input_fasta_file_path/$input_fasta_file_short";
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

  my $mech = WWW::Mechanize->new( autocheck => 1 );
  $mech->get( "http://indra.mullins.microbiol.washington.edu/DIVEIN/insites.html" );

  my $result = $mech->submit_form(
                                  form_name => 'alignmentForm',
                                  fields    => { local => 'DIVEIN', seqFile => $input_fasta_file, seqRadio => 'fasta', datatype => 'DNA' }
                       );
  
  # Have to submit it again, ie my $result2 = $result->submit();
  my $content = $mech->content();
  if( $DEBUG ) {
    print "OK1, \$content is $content\n";
  }
  $mech->select('seqName',{n=>[1..1000]});
  my $result2 = $mech->submit_form(
     form_name => 'grpForm'
  );
  $content = $result2->content();

   #my @links = $mech->find_all_links(
   #   tag => "a", text_regex => qr/\bdownload\b/i );
#  my $ua = LWP::UserAgent->new();
#  my $req = POST 'http://indra.mullins.microbiol.washington.edu/cgi-bin/DIVEIN/insites/insites_grp.cgi',
#                [ local => 'DIVEIN', seqFile => $input_fasta_file_contents, seqRadio => 'fasta', datatype => 'DNA' ];
#  my $content = $ua->request( $req )->as_string;

  #if( $DEBUG ) {
    print "OK2\n \$content is $content\n";
  #}
  ## ERE I AM!!
## Aligned informative sites:</td><td>None
  my ( $job_id ) = ( $content =~ /Your job id is (\d+)\./ );

  # Save all the files to  $output_path_dir
  #my @links = $mech->find_all_links(
  #   tag => "a", text_regex => qr/\bdownload\b/i );

  if( $VERBOSE ) {
    print "JOB ID IS $job_id\n";
  }

  #return 1;

  my $privContent = undef;
  $mech->get( "http://indra\.mullins\.microbiol\.washington.edu/cgi-bin/DIVEIN/insites/download\.cgi?id=$job_id&ext=_priv\.txt&&local=DIVEIN");
  $privContent = $mech->content();
  #print( "got $privContent" );
  while( !defined( $privContent ) || $content =~ /No such file/ ) {
    sleep( 1 );
    print( "trying again to get the private sites file" );
    $mech->get( "http://indra\.mullins\.microbiol\.washington.edu/cgi-bin/DIVEIN/insites/download\.cgi?id=$job_id&ext=_priv\.txt&&local=DIVEIN");
    $privContent = $mech->content();
  }

  my $privFile =  $output_path_dir . "/" . $input_fasta_file_short_nosuffix . "_privateSites.txt";
  if( $VERBOSE ) { print "Opening file \"$privFile\" for writing.."; }
  unless( open privFileFH, ">$privFile" ) {
      warn "Unable to open output file \"$privFile\": $!\n";
      return 1;
    }
  print privFileFH $privContent;
  close( privFileFH );
  if( $VERBOSE ) { print ".done\n"; }

  $mech->get("http://indra\.mullins\.microbiol\.washington\.edu/cgi-bin/DIVEIN/insites/download\.cgi?id=$job_id&ext=.txt&local=DIVEIN");
  my $informativeSitesContent = $mech->content();
  print( "got $informativeSitesContent" );
  while( !defined( $informativeSitesContent ) || ( $informativeSitesContent =~ /No such file/ ) ) {
    print( "trying again to get the informative sites file." );
    sleep( 1 );
    $mech->get("http://indra\.mullins\.microbiol\.washington\.edu/cgi-bin/DIVEIN/insites/download\.cgi?id=$job_id&ext=.txt&local=DIVEIN");
    $informativeSitesContent = $mech->content();
  }

  my $informativeSitesFile =  $output_path_dir . "/" . $input_fasta_file_short_nosuffix . "_informativeSites.txt";
  if( $VERBOSE ) { print "Opening file \"$informativeSitesFile\" for writing.."; }
  unless( open informativeSitesFileFH, ">$informativeSitesFile" ) {
      warn "Unable to open output file \"$informativeSitesFile\": $!\n";
      return 1;
    }
  print informativeSitesFileFH $informativeSitesContent;
  close( informativeSitesFileFH );
  if( $VERBOSE ) { print ".done\n"; }

  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  return 0;
} # runInSitesOnline(..)

runInSitesOnline( @ARGV );

1;


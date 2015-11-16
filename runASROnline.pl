#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) runASROnline
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##  Description:
##      Script for running ASR using the online datamonkey tools (at
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

sub runASROnline {
  @ARGV = @_;

  sub runASROnline_usage {
    print "\trunASROnline [-DV] <input_fasta_file> <outgroup_fasta_file> [<output_dir>]\n";
    exit;
  }

  # This means -D and -V are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # But first reset the opt vars.
  ( $opt_D, $opt_V ) = ();
  if( not getopts('DV') ) {
    runASROnline_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;
  
  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }
  
  my $input_fasta_file = shift @ARGV || runASROnline_usage();
  my ( $input_fasta_file_path, $input_fasta_file_short ) =
    ( $input_fasta_file =~ /^(.*?)\/([^\/]+)$/ );
  unless( $input_fasta_file_short ) {
    $input_fasta_file_short = $input_fasta_file;
    $input_fasta_file_path = ".";
  }
  my ( $input_fasta_file_short_nosuffix, $input_fasta_file_suffix ) =
    ( $input_fasta_file_short =~ /^([^\.]+)(\..+)?$/ );

  my $outgroup_fasta_file = shift @ARGV || runASROnline_usage();
  my ( $outgroup_fasta_file_path, $outgroup_fasta_file_short ) =
    ( $outgroup_fasta_file =~ /^(.*?)\/([^\/]+)$/ );
  unless( $outgroup_fasta_file_short ) {
    $outgroup_fasta_file_short = $outgroup_fasta_file;
    $outgroup_fasta_file_path = ".";
  }
  my ( $outgroup_fasta_file_short_nosuffix, $outgroup_fasta_file_suffix ) =
    ( $outgroup_fasta_file_short =~ /^([^\.]+)(\..+)?$/ );

  my $output_path_dir = shift @ARGV ||
    $input_fasta_file_path . "/" . $input_fasta_file_short_nosuffix . "_hiv-founder-id_resultsDir";
  # Remove the trailing "/" if any
  if( defined( $output_path_dir ) ) {
    ( $output_path_dir ) = ( $output_path_dir =~ /^(.*[^\/])\/*$/ );
  }

  # Make one unholy hybrid.  NOTE THIS ASSUMES THE BEST OF THINGS.  IT ASSUMES THE FASTA FILES ARE IN THE SAME ALIGNMENT - IE THEY ARE ALREADY ALIGNED TO ONE ANOTHER.  ALSO THE OUTGROUP FILE HAS A HEADER WE WILL IGNORE.

  my $input_fasta_file_with_outgroup_short_nosuffix = $input_fasta_file_short_nosuffix . "_withOutgroup";
  my $input_fasta_file_with_outgroup_short = $input_fasta_file_with_outgroup_short_nosuffix . $input_fasta_file_suffix;
  my $input_fasta_file_with_outgroup =  $output_path_dir . '/' . $input_fasta_file_with_outgroup_short;
  
  if( $VERBOSE ) { print "Opening file \"$input_fasta_file_with_outgroup\" for writing.."; }
  unless( open input_fasta_file_with_outgroupFH, ">$input_fasta_file_with_outgroup" ) {
    warn "Unable to open output file \"$input_fasta_file_with_outgroup\": $!\n";
    return 1;
  }
  if( $VERBOSE ) { print "Opening file \"$outgroup_fasta_file\" for reading.."; }
  unless( open OUTGROUP_FH, $outgroup_fasta_file ) {
    warn( "Unable to open outgroup file \"$outgroup_fasta_file\': $!\n" );
    return 1;
  }
  my $line;
  while( $line = <OUTGROUP_FH> ) {
    if( $line =~ /^\s*>/ ) {
      # It's the header.  Instead, print ">Outgroup".
      print input_fasta_file_with_outgroupFH ">Outgroup\n";
      next;
    }
    # Otherwise, just copy it over.
    print input_fasta_file_with_outgroupFH $line;
  }
  print input_fasta_file_with_outgroupFH "\n";
  if( $VERBOSE ) { print ".done.\nClosing file \"$outgroup_fasta_file\".."; }
  close( OUTGROUP_FH );
  if( $VERBOSE ) { print ".done.\n"; }
  if( $VERBOSE ) { print ".done.\nClosing file \"$input_fasta_file_with_outgroup\".."; }
  close( input_fasta_file_with_outgroupFH );
  if( $VERBOSE ) { print ".done.\n"; }
  
  # Now put the rest of the alignment on there.
  `cat $input_fasta_file >> $input_fasta_file_with_outgroup`;

  my $mech = WWW::Mechanize->new( autocheck => 1 );
  $mech->get( "http://www.datamonkey.org/dataupload.php" );
  #print( $mech->content() );
  my $result = $mech->submit_form(
                                  form_name => 'uploadform',
                                  fields    => { datatype => '0', code => '0', upfile => $input_fasta_file_with_outgroup }
                       );
  
  # Have to submit it again, ie my $result2 = $result->submit();
  my $content = $mech->content();
  if( $DEBUG ) {
    print "OK1\n \$content is $content\n";
  }
  my ( $stop_codons_found ) = ( $content =~ /stop codons found/ );
  if( $stop_codons_found ) {
    die( "STOP Codons found [one or more of: (TGA, TAA, TAG)]." );
  }
  my ( $datamonkey_filename ) = ( $content =~ m/name='filename' value=\'([^\']+)\'/ );
#   <FORM method='POST' enctype='multipart/form-data' action='http://www.datamonkey.org/cgi-bin/datamonkey/finishUpload.pl'>
# <DIV style = 'margin:10px; text-align:center;'><input type='Hidden' name='filename' value='upload.21231639675385.1'><input type='Hidden' name='genCodeID' value='0'>
# <input type='Submit' value='Proceed to the analysis menu' style = 'background-color:purple; color:white; font-size:18px;'> <!--<a href='http://www.datamonkey.org/help/tree.php' target = '_blank' class = 'INFO'>Help</a>--> </DIV></FORM>	
  my $result2 = $mech->submit_form(
                                   #form_number => 2,
                                   fields => {
                                     genCodeID => '0',
                                     filename => $datamonkey_filename
                                             }

  );
  $content = $result2->content();

  if( $DEBUG ) {
    print "OK2\n \$content is $content\n";
  }
  exit( 1 );
  
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
                                              method => 22, #  (20 is for 'SBP', 21 is for 'GARD', and 22 is for 'ASR')
                                              root => "OUTGROUP", # This is the sequence we created; Its name has become uppercase.
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

  
  ## KEEP
  ## # ASR is done when the Datamonkey Job Status page changes to include the text:
  ## # 'ASR recombination analysis was run '
  ## $mech->get( "http://www.datamonkey.org/cgi-bin/datamonkey/jobStatus.pl?file=$datamonkey_filename" );
  ## $content = $mech->content();
  ## my ( $gard_results_are_done ) = ( $content =~ m/ASR recombination analysis was run / );
  ## #    /spool/upload.${datamonkey_filename}_gard.php

  # SBP is done when the Datamonkey Job Status page changes to include the text:
  # 'Single Breakpoint Recombination was run '
  $mech->get( "http://www.datamonkey.org/cgi-bin/datamonkey/jobStatus.pl?file=$datamonkey_filename" );
  $content = $mech->content();
  if( $DEBUG ) {
    print "OK5\n \$content is $content\n";
  }
  ## ERE I AM.
  exit( 1 );
  my $asr_results_are_done = ( $content =~ m/ASR recombination analysis was run / );
  while( !$asr_results_are_done ) {
    sleep( 10 );
    print( "Trying again to get the ASR results..\n" );
    $mech->get( "http://www.datamonkey.org/cgi-bin/datamonkey/jobStatus.pl?file=$datamonkey_filename" );
    $content = $mech->content();
    if( $DEBUG ) {
      print "OK5++\n";
    }
    $asr_results_are_done = ( $content =~ m/Single Breakpoint Recombination was run / );
  }
  print( "ASR is done.\n" );
  
  $mech->get( "http://www.datamonkey.org/spool/${datamonkey_filename}_asr.php" );
  $content = $mech->content();
  if( $DEBUG ) {
    print "OK6\n \$content is $content\n";
  }

#<TABLE><TR CLASS = 'TRReport1'><TD>Model</TD><TD>010010</TD></TR><TR CLASS = 'TRReport2'><TD>Potential Breakpoints</TD><TD>945</TD></TR><TR CLASS = 'TRReport1'><TD>Processed Breakpoints</TD><TD>945</TD></TR><TR CLASS = 'TRReport2'><TD>Run time</TD><TD>00:00:50</TD></TR></TABLE></DIV><DIV class = 'RepClassSM'><b>Recombination report</b><br><i>Small sample AIC (cAIC) is the recommended default criterion</i>
# <p><table><tr class = 'HeaderClassSM'><td>IC <span class = 'INFO' onmouseover = "Tip('Which information criterion is used for inference')">?</span></td><td>Recombination <span class = 'INFO' onmouseover = "Tip('Was recombination inferred using this IC?')">?</span></td><td>IC improvement <span class = 'INFO' onmouseover = "Tip('How many points did the recombinant model improve the IC over the base (single tree) model?')">?</span></td><td>Breakpoint location<span class = 'INFO' onmouseover = "Tip('Where is the most likely breakpoint located?')">?</span></td><td>Model averaged support<span class = 'INFO' onmouseover = "Tip('Model averaged confidence for having recombination in the alignment')">?</span></td></tr><tr class = 'TRReport1'><td>AIC</td><td>No</td><td>N/A</td><td>N/A</td><td> 0.00%</td></tr>
# <tr class = 'TRReport2'><td>cAIC</td><td>No</td><td>N/A</td><td>N/A</td><td> 0.00%</td></tr>
  # <tr class = 'TRReport1'><td>BIC</td><td>No</td><td>N/A</td><td>N/A</td><td> 0.00%</td></tr></table>
  # Recombination ?	IC improvement ?	Breakpoint location?	Model averaged support?
#     my ( $asr_caic_support_for_recombination, $asr_caic_IC_improvement, $asr_caic_breakpoint_location, $asr_caic_model_averaged_support ) =
#     ( $content =~ m/cAIC\<\/td\>\<td\>(.+)\<\/td\>\<td\>(.+)\<\/td\>\<td\>(.+)\<\/td\>\<t\d>(.+)\<\/td\>\<\/tr\>/ );


  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  return 0;
} # runASROnline(..)

runASROnline( @ARGV );

1;


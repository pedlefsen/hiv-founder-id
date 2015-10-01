#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) runInSitesOffline
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##      Ted Holzman              tholzman@fhcrc.org
##  Description:
##      Script based on Paul Edlefsen's script for screen-scraping
##      the Mullins lab web site for obtaining the informative vs private
##      sites ratio statistic: see  getInSitesStat.  Note that this creates 
##      output files in subdirectories named after the input fasta file name.
##
##      This version does not do screen scraping.  It uses the perl modules 
##      that drive the Mullins divein web-page directly.  There is no 
##      CGI or Mechanize in this version.
##      
##      NOTE that there is a bug (as of September, 2015) in the
##      inSites code online, in which flanking gaps are not printed in
##      the output (informative sites) table.  THIS MUST BE FIXED WHEN
##      READING/USING THE FILE.
##
###******************************************************************************

use Getopt::Std; # for getopts
use File::Path qw( make_path );
use Path::Tiny;
use DateTime;

use strict;
use vars qw( $opt_D $opt_V );
use vars qw( $VERBOSE $DEBUG );

#TAH use the same arguments and options as in Paul's code 
sub runInSitesOffline { 
  @ARGV = @_;

  sub runInSitesOffline_usage {
    print "\trunInSitesOnline [-DV] <input_fasta_file> [<output_dir>]\n";
    exit;
  }

  # This means -D and -V are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # But first reset the opt vars.
  ( $opt_D, $opt_V ) = ();
  if( not getopts('DV') ) {
    runInSitesOffline_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;

  # if verbose, use unbuffered immediate (autoflush) output
  # save the current state of the flush state ($|) of STDOUT
  # and reset to status quo ante before exiting function
  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }

  # get the input file name 
  my $input_fasta_file = shift @ARGV || runInSitesOffline_usage();

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
    ## removed:  TAH 9/15
    ## The function of this hack is moved to divein's CleanString function
    
  open(my $seqFile_handle, "<",$input_fasta_file) or die "Couldn't open $input_fasta_file for reading\n";

  ##TAH this section emulates the <SUBMIT> action from the insites.cgi script which is
  ##essentially the computational (not display) functions of insites_grp.cgi
    
  #id is date-time instead of seq. num from divein
  my $id = $dt->ymd('')."_".$dt->hms('');
  my $datatype = "nt";
  my $seqRadio = "fasta";
  my $rsRadio = '';
  my $grpRadio = 'n';
  my $sortRadio = 'n';
  my $seqGrpHashRef;
  my @nas = qw (A C G T);

  my $seqFileLines = GetFileLines ($seqFile_handle);
  # get sequence info such as sequence number, length, array of seq name and array of seq name and sequence
  my $seqInfo = GetSequences ($seqFileLines, $seqRadio);	



my $seqNum = shift @$seqInfo;
my $seqLen = shift @$seqInfo;
my $seqNamesRef = shift @$seqInfo;
my $stdNamesRef = shift @$seqInfo;











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
  my ( $no_informative_sites ) = ( $content =~ /Aligned informative sites:<\/td><td>None/ );
  my ( $job_id ) = ( $content =~ /Your job id is (\d+)\./ );

  # Save all the files to  $output_path_dir
  #my @links = $mech->find_all_links(
  #   tag => "a", text_regex => qr/\bdownload\b/i );

  if( $VERBOSE ) {
    print "JOB ID IS $job_id\n";
    if( $no_informative_sites ) {
      print "NO INFORMATIVE SITES\n";
    }
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

  my $informativeSitesContent = "";
  if( $no_informative_sites ) {
    $informativeSitesContent = "";
  } else {
    $mech->get("http://indra\.mullins\.microbiol\.washington\.edu/cgi-bin/DIVEIN/insites/download\.cgi?id=$job_id&ext=.txt&local=DIVEIN");
    $informativeSitesContent = $mech->content();
    print( "got $informativeSitesContent" );
    while( !defined( $informativeSitesContent ) || ( $informativeSitesContent =~ /No such file/ ) ) {
      print( "trying again to get the informative sites file." );
      sleep( 1 );
      $mech->get("http://indra\.mullins\.microbiol\.washington\.edu/cgi-bin/DIVEIN/insites/download\.cgi?id=$job_id&ext=.txt&local=DIVEIN");
      $informativeSitesContent = $mech->content();
    }
  }

  my $informativeSitesFile =  $output_path_dir . "/" . $input_fasta_file_short_nosuffix . "_informativeSites.txt";
  if( $VERBOSE ) { print "Opening file \"$informativeSitesFile\" for writing.."; }
  unless( open informativeSitesFileFH, ">$informativeSitesFile" ) {
      warn "Unable to open output file \"$informativeSitesFile\": $!\n";
      return 1;
    }
  ## NOTE that there is a bug (as of September, 2015) in the inSites code online, in which flanking gaps are not printed in the output table.  THIS MUST BE FIXED WHEN READING/USING THE FILE.
  print informativeSitesFileFH $informativeSitesContent;
  close( informativeSitesFileFH );
  if( $VERBOSE ) { print ".done\n"; }

  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  return 0;
} # runInSitesOffline(..)

# borrowed and modified from insites_grp.cgi
# TAH 9/15
# usage <array of file trimmed lines> = GetFileLines(<filehandle>)
sub GetFileLines {
	my $fh = shift;
	my $line = "";
	my @buffer = <$fh>;
 	foreach my $element (@buffer) {
 		$line .= $element;
 	}
 	if ($line =~ /\r\n/) {
		$line =~ s/\r//g;
	}elsif ($line =~ /\r/) {
		$line =~ s/\r/\n/g;
	}
	my @fileLines = split /\n/, $line;
	return \@fileLines;
}
# stolen and modified from Divein.pm
# added Paul's hack to to change ! into -
# TAH 9/15
sub CleanString {
	my $string = shift;
	$string =~ s/^\s+//;	# remove leading spaces
	$string =~ s/\s+$//;	# remove ending spaces
        $string =~ s/\!/-/g;    # change bangs to gaps TAH
	return $string;
}

# stolen and modified from Divein.pm
# de-CGIed, and saving only the code pertinent to fasta files
# TAH 9/15
# usage: <info-structure> = GetSequences(<raw file lines array>,<file type code>)
sub GetSequences {
   if ($seqFileRadio ne "fasta") die "This program only deals with fasta files\n";
   my ($seqFileLines, $seqFileRadio) = @_;
   my ($seqName, %seqNameStatus, @seqNames, @stdnameNseqs, @seqInfo, %nameSeq, %countName);
   my $seqNum = my $seqLen = my $seqCount = my $nexusFlag = my $phylipFlag = 
      my $seqStartFlag = my $ntaxFlag = my $ncharFlag = my $fastSeqFlag = my $count =
         0;
   my fastaFlag = 1;
   
   foreach my $line (@$seqFileLines) {
      next if $line =~ /^\s*$/;  #skip blank lines
      $line = CleanString ($line);
      if ($line !~ /^>/) die "Bad FASTA file.  First non-blank line doesn't start with a >";
      if ($line =~ /^>(\S+)/) {#first line of a new sequence
         $seqCount++;
	 $seqNum = $seqCount;
	 $seqName = $1;
	 push @seqNames, $seqName;
	 $fastSeqFlag = 0;
      }else { #sequence lines of current sequence
	 $line =~ s/\s//g;	# remove spaces that may be in sequences
	 unless ($line =~ /^[A-Za-z\-\.\?]+$/) {
	    my @nas = split //, $line;
	    foreach my $na (@nas) {
	       unless ($na =~ /[A-Za-z\-\.\?]/) {
	          die "Error: couldn't recognize character $na in sequence $seqName. Please check the sequence file.";
	       }
	    }								
	 }
	 if (!$fastSeqFlag) {
	    $nameSeq{$seqName} = "";
	    $fastSeqFlag = 1;
	 }
	 $nameSeq{$seqName} .= $line;				
      }
    }
    $seqLen = length $nameSeq{$seqNames[0]};	# set alignment length to be the length of first sequence
#   print "length: $seqLen<br>";
    if (!@seqNames) die "Error: Couldn't get sequence names from the input sequence file. Please check the sequence file.\n";
    my %seqNamesHash;
    foreach my $seqName (@seqNames) {
       my $ciName = uc $seqName;
       # name is case-insensitive, because hyphy treats lower- and upper-case same
       if (!$seqNamesHash{$ciName}) {
          $seqNamesHash{$ciName} = 1;
       }else {
          die "Error: At least two sequences have the same name of $seqName in sequence alignment file. It may be caused by the space(s) in sequence name or case-insensitve of the name.\nPlease check your input sequence file to make sure the unique sequence name or no space in sequence name.\n";
       }		
    }
    if ($seqCount != $seqNum) {
       die "Error: Number of sequences is not equal to pre-defined sequence number of $seqNum in your uploaded sequence alignment file. 
		It may caused by the duplicated sequence names in your alignment. Please check the sequence number or remove the duplicates and upload your file again.\n";
    }
    foreach my $seqName (@seqNames) {	# check each sequence length in the alignment
       if (length $nameSeq{$seqName} != $seqLen) {						
          die "Error: Lengths of sequences are not same among the alignment. It may caused by the duplicated sequence names in your alignment. \nPlease check your input sequence file to make sure that sequences are aligned or there are no duplicates.\n";
       }
    }
    my @stdSeqNames;
    foreach my $seqName (@seqNames) {
       my $seq = $nameSeq{$seqName};
       my $stdName = $seqName;
       $stdName =~ s/\W/_/g;	# replace non-letters with "_"
       my $stdnameNseq = $stdName."\t".$seq;
       push @stdnameNseqs, $stdnameNseq;
       push @stdSeqNames, $stdName;
    }
    push @seqInfo, $seqNum, $seqLen, \@seqNames, \@stdSeqNames, \@stdnameNseqs;
    return \@seqInfo;
}









runInSitesOnline( @ARGV );

1;


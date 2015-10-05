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
    
  my $startTime = time();
  my $rand = int (rand (90)) + 10;
  my $id = $startTime.$rand;
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
  my $uploadDir = $output_path_dir;
  my $seqNameNseq = shift @$seqInfo;
  unshift @$seqNameNseq, $seqNum."\t".$seqLen;

  my $seqFile = $uploadDir.'/'.$id;
  WriteFile ($seqFile, $seqNameNseq);
   
  ### Now we simulate insites.cgi, the action of the grpForm form
  my $uploadSeqFile = $seqFile;
  my $uploadseqFile = $uploadDir.'/'.$id;
  my %nameSeq = ();
  $seqLen = 0;

  open SEQ, $uploadseqFile or die "couldn't open $uploadseqFile: $!\n";
  while (my $line = <SEQ>) {
     chomp $line;
     next if $line =~ /^\s*$/;
     if ($line =~ /^\d+\s(\d+)$/) {
        $seqLen = $1;
     } else {
	my ($name, $seq) = split /\t/, $line;
	$nameSeq{$name} = $seq;
     }
  }
  close SEQ;

  $seqNum = 0;
  my @grpNames = my @grpSeqs = my %seqGrp = ();

  my @seqs = $seqNamesRef;
  push @grpSeqs, \@seqs;
  $seqNum = scalar @seqs;
  foreach my $seqName (@seqs) {
     push @grpSeqs, \@seqs;
     $seqNum = scalar @seqs;
     foreach my $seqName (@seqs) {
        $seqGrp{$seqName} = "default";
     }
   }

   my @seqNameNseqs = ();
   for (my $i=0; $i<@grpSeqs; $i++) {
      foreach my $seqName (@{$grpSeqs[$i]}) {
         my $nameNseq = $seqName."\t".$nameSeq{$seqName};
 	 push @seqNameNseqs, $nameNseq;
      }
   }
   my $refName = "Consensus";
   my $refSeq = GetConsensus (\@seqNameNseqs, $seqLen);
   $seqFile = $uploadseqFile."_seq4insites.phy";
   $seqNum += 1;
   my $refNameNseq = $refName."\t".$refSeq;
   unshift @seqNameNseqs, $refNameNseq;
   unshift @seqNameNseqs, $seqNum."\t".$seqLen;
   WriteFile ($seqFile, \@seqNameNseqs);


   return 1;
#/  my ( $no_informative_sites ) = ( $content =~ /Aligned informative sites:<\/td><td>None/ );
#/  my ( $job_id ) = ( $content =~ /Your job id is (\d+)\./ );
#/
#/  if( $VERBOSE ) {
#/    print "JOB ID IS $job_id\n";
#/    if( $no_informative_sites ) {
#/      print "NO INFORMATIVE SITES\n";
#/    }
#/  }
#/
#/
#/  my $privContent = undef;
#/  $mech->get( "http://indra\.mullins\.microbiol\.washington.edu/cgi-bin/DIVEIN/insites/download\.cgi?id=$job_id&ext=_priv\.txt&&local=DIVEIN");
#/  $privContent = $mech->content();
#/
#/  while( !defined( $privContent ) || $content =~ /No such file/ ) {
#/    sleep( 1 );
#/    print( "trying again to get the private sites file" );
#/    $mech->get( "http://indra\.mullins\.microbiol\.washington.edu/cgi-bin/DIVEIN/insites/download\.cgi?id=$job_id&ext=_priv\.txt&&local=DIVEIN");
#/    $privContent = $mech->content();
#/  }

#/  my $privFile =  $output_path_dir . "/" . $input_fasta_file_short_nosuffix . "_privateSites.txt";
#/  if( $VERBOSE ) { print "Opening file \"$privFile\" for writing.."; }
#/  unless( open privFileFH, ">$privFile" ) {
#/      warn "Unable to open output file \"$privFile\": $!\n";
#/      return 1;
#/    }
#/  print privFileFH $privContent;
#/  close( privFileFH );
#/  if( $VERBOSE ) { print ".done\n"; }
#/
#/  my $informativeSitesContent = "";
#/  if( $no_informative_sites ) {
#/    $informativeSitesContent = "";
#/  } else {
#/    $mech->get("http://indra\.mullins\.microbiol\.washington\.edu/cgi-bin/DIVEIN/insites/download\.cgi?id=$job_id&ext=.txt&local=DIVEIN");
#/    $informativeSitesContent = $mech->content();
#/    print( "got $informativeSitesContent" );
#/    while( !defined( $informativeSitesContent ) || ( $informativeSitesContent =~ /No such file/ ) ) {
#/      print( "trying again to get the informative sites file." );
#/      sleep( 1 );
#/      $mech->get("http://indra\.mullins\.microbiol\.washington\.edu/cgi-bin/DIVEIN/insites/download\.cgi?id=$job_id&ext=.txt&local=DIVEIN");
#/      $informativeSitesContent = $mech->content();
#/    }
#/  }

#/  my $informativeSitesFile =  $output_path_dir . "/" . $input_fasta_file_short_nosuffix . "_informativeSites.txt";
#/  if( $VERBOSE ) { print "Opening file \"$informativeSitesFile\" for writing.."; }
#/  unless( open informativeSitesFileFH, ">$informativeSitesFile" ) {
#/      warn "Unable to open output file \"$informativeSitesFile\": $!\n";
#/      return 1;
#/    }
  ## NOTE that there is a bug (as of September, 2015) in the inSites code online, in which flanking gaps are not printed in the output table.  THIS MUST BE FIXED WHEN READING/USING THE FILE.
#/  print informativeSitesFileFH $informativeSitesContent;
#/  close( informativeSitesFileFH );
#/  if( $VERBOSE ) { print ".done\n"; }

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
# toDo
# save only the code pertinent to fasta files
# TAH 10/15
# usage: <info-structure> = GetSequences(<raw file lines array>,<file type code>)
sub GetSequences {
	my ($seqFileLines, $seqFileRadio) = @_;
	my ($seqName, %seqNameStatus, @seqNames, @stdnameNseqs, @seqInfo, %nameSeq, %countName);
	my $seqNum = my $seqLen = my $seqCount = my $nexusFlag = my $phylipFlag = my $fastaFlag = my $seqStartFlag = my $ntaxFlag = my $ncharFlag = my $fastSeqFlag = my $count = 0;
	foreach my $line (@$seqFileLines) {
		next if $line =~ /^\s*$/;
		$line = CleanString ($line);
		if ($seqFileRadio eq "nexus") {
			if ($nexusFlag == 0 && $line =~ /^\#NEXUS$/i) {
				$nexusFlag = 1;	
			}elsif ($nexusFlag) {
				if($line =~ /NTAX=(\d+)/i) {	# handled the case that ntax and nchar are not in the same line
					$seqNum = $1;
					$ntaxFlag = 1;
				}
				if($line =~ /NCHAR=(\d+)/i) {	# handled the case that ntax and nchar are not in the same line
					$seqLen = $1;
					$ncharFlag = 1;
				}
				if ($line =~ /^MATRIX$/i) {
					$seqStartFlag = 1;
				}elsif ($line =~ /^\;$/) {
					$seqStartFlag = 0;
					if (!$ntaxFlag) {
						$seqNum = $seqCount;
					}
				}elsif ($seqStartFlag) {
					unless ($line =~ /^\[\s+/) {	# ignore the line beginning with [
						$line =~ s/\[\d+\]$// if ($line =~ /\[\d+\]$/);	# remove [\d+] on the end of line
						if ($line =~ /^(\S+)\s+(.*)$/) {	# requires no space in sequence names, there could be spaces in sequences
							
							my $seqName = $1;
							my $seq = $2;
							$seq =~ s/\s//g;	# remove spaces that may be in sequences
							unless ($seq =~ /^[A-Za-z\-\.\?]+$/) {
								my @nas = split //, $seq;
								foreach my $na (@nas) {
									unless ($na =~ /[A-Za-z\-\.\?]/) {
										print "<p>Error: couldn't recognize character $na in sequence $seqName. Please check the sequence file.</p>";
										PrintFooter();
									}
								}								
							}
											
							if (!$seqNameStatus{$seqName}) {
								push @seqNames, $seqName;
								$seqNameStatus{$seqName} = 1;
								$nameSeq{$seqName} = '';
								$seqCount++;
							}else {
								#print "<p>$seqName</p>";
							}
							
							$nameSeq{$seqName} .= $seq;
						}
					} 				
				}
			}else {
				print "<p>The upload sequence file is not a nexus file. The nexus file must start with #NEXUS. Please check the file and upload again.</p>";
				PrintFooter();
			}
		}elsif ($seqFileRadio eq "phylip" || $seqFileRadio eq "example") {	# phylip file
			if ($phylipFlag == 0 && $line =~ /^\s*(\d+)\s+(\d+)/) {
				$seqNum = $1;
				$seqLen = $2;
				$phylipFlag = 1;
			}elsif ($phylipFlag) {
				if ($count < $seqNum) {
					$line =~ /^(\S+)\s+(.*)$/;	# requires no space in sequence names, there could be spaces in sequences
					my $seqName = $1;
					my $seq = $2;
					$seq =~ s/\s//g;	# remove spaces that may be in sequences
					unless ($seq =~ /^[A-Za-z\-\.\?]+$/) {
						my @nas = split //, $seq;
						foreach my $na (@nas) {
							unless ($na =~ /[A-Za-z\-\.\?]/) {
								print "<p>Error: couldn't recognize character $na in sequence $seqName. Please check the sequence file.</p>";
								PrintFooter();
							}
						}								
					}
					push @seqNames, $seqName;
					$countName{$count} = $seqName;
					$nameSeq{$seqName} = $seq;
					$seqCount++;
				}else {	# interleaved 
					my $index = $count % $seqNum;
					$line =~ s/\s//g;	# remove all spaces, the line only contains sequence
					my $seqName = $countName{$index};
					$nameSeq{$seqName} .= $line;
				}				
				$count++;
			}else {
				print "<p>The upload sequence file is not a phylip file. Please check the file and upload again.</p>";
				PrintFooter();
			}
		}else {	# fasta file, the code can handle both interval and sequential formats
			if (!$fastaFlag) {	# check for fasta format
				if ($line !~ /^>/) {
					print "<p>The upload sequence file is not a fasta file. Please check the file and upload again.</p>";
					PrintFooter();
				}else {
					$fastaFlag = 1;
				}	
			}
			
			if ($line =~ /^>(\S+)/) {
				$seqCount++;
				$seqNum = $seqCount;
				$seqName = $1;
				push @seqNames, $seqName;
				$fastSeqFlag = 0;
			}else {
				$line =~ s/\s//g;	# remove spaces that may be in sequences
				unless ($line =~ /^[A-Za-z\-\.\?]+$/) {
					my @nas = split //, $line;
					foreach my $na (@nas) {
						unless ($na =~ /[A-Za-z\-\.\?]/) {
							print "<p>Error: couldn't recognize character $na in sequence $seqName. Please check the sequence file.</p>";
							PrintFooter();
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
	}
	
	if ($seqFileRadio eq "nexus") {
		if (!$seqNum && !$seqLen) {
			print "<p>Error: Couldn't get information of pre-defined sequnece number and length. It is probably caused by missing statement of dimensions 
			in your nexus file. Please check the sequence file.</p>";
			PrintFooter();
		}
	}elsif ($seqFileRadio eq "fasta") {
		$seqLen = length $nameSeq{$seqNames[0]};	# set alignment length to be the length of first sequence
	}
#	print "length: $seqLen<br>";
	if (!@seqNames) {
		print "<p>Error: Couldn't get sequence names from the input sequence file. Please check the sequence file.</p>";
		PrintFooter();
	}else {
		my %seqNamesHash;
		foreach my $seqName (@seqNames) {
			#if (length $seqName > 30) {	# check length of sequence name, when compile PhyML, user can set it. now set to maximum 30 characters
			#	print "<p>Error: The length of sequence name exceeds the maximum length of 30 characters.</p>";
			#	PrintFooter();
			#}
			my $ciName = uc $seqName;
			# name is case-insensitive, because hyphy treats lower- and upper-case same
			if (!$seqNamesHash{$ciName}) {
				$seqNamesHash{$ciName} = 1;
			}else {
				print "<p>Error: At least two sequences have the same name of $seqName in sequence alignment file. It may be caused by the space(s) in sequence name or case-insensitve of the name. ";
				print "Please check your input sequence file to make sure the unique sequence name or no space in sequence name.</p>";
				PrintFooter();
			}		
		}
	}
	
	if ($seqCount != $seqNum) {
		print "<p>Error: Number of sequences is not equal to pre-defined sequence number of $seqNum in your uploaded sequence alignment file. 
		It may caused by the duplicated sequence names in your alignment. Please check the sequence number or remove the duplicates and upload your file again.</p>";
		PrintFooter();
	}
	
	foreach my $seqName (@seqNames) {	# check each sequence length in the alignment
		if (length $nameSeq{$seqName} != $seqLen) {						
			print "<p>Error: Lengths of sequences are not same among the alignment. It may caused by the duplicated sequence names in your alignment. ";
			print "Please check your input sequence file to make sure that sequences are aligned or there are no duplicates.</p>";
			PrintFooter();
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
	
#	my $seqNumNLen = $seqNum."\t".$seqLen;
#	unshift @stdnameNseqs, $seqNumNLen;	
#	my $datasize = $seqNum * $seqLen;
	push @seqInfo, $seqNum, $seqLen, \@seqNames, \@stdSeqNames, \@stdnameNseqs;
	return \@seqInfo;
}


## Stolen directly from Divein.pm
## TAH 10/1/15
sub WriteFile {
   my ($uploadFile, $fileLines) = @_;
   open OUT, ">$uploadFile" or die "couldn't open $uploadFile: $!\n";
   foreach my $line (@$fileLines) {
      print OUT $line,"\n";
   }
   close OUT;
}

## Stolen directly from Divein.pm
## To Do:  remove unused portions, e.g. protein stuff.
## maybe normalize spacing.
## TAH 10/15
sub GetConsensus {
	my ($seqNameNseq, $nchar) = @_;
	my (@seqNames, $seqArr);
	my $element = scalar @$seqNameNseq;
	for (my $i = 0; $i < $element; $i++) {
		my ($seqName, $seq) = split /\t/, $seqNameNseq->[$i];
		push @seqNames, $seqName;
		$seq = uc $seq;
		my $beginFlag = my $terminalFlag = 0;
		my $nameNformatedseq = $seqName."\t";
		for (my $j = 0; $j < $nchar; $j++) {
			my $aa = substr($seq, $j, 1);				
			# deal with leading gaps
			if ($j == 0 && $aa eq "-") {
				$beginFlag = 1;
				$aa = " ";
			}elsif ($beginFlag == 1 && $aa eq "-") {
				$aa = " ";
			}elsif ($aa ne "-") {
				$beginFlag = 0;
			}
			# deal with termianl gaps
			if (substr ($seq, $j) =~ /^\-+$/) {
				$terminalFlag = 1;
			}
			if ($terminalFlag == 1) {
				$aa = " ";
			}
			
			if ($i == 0) {
				$seqArr->[$i]->[$j] = $aa;
			}else {				
				if ($aa eq ".") {
					$seqArr->[$i]->[$j] = $seqArr->[0]->[$j];
				}else {
					$seqArr->[$i]->[$j] = $aa;
				}
			}
		}
	}
	my @consAas;
	for (my $i = 0; $i < $nchar; $i++) {
		my %aaCount;
		my $blankCount = 0;
		for (my $j = 0; $j < $element; $j++) {
			my $aa = $seqArr->[$j]->[$i];
			unless ($aa eq "?") {
				if ($aa eq " ") {	# leading or ending gaps
					$blankCount++;
				}else {
					if (!$aaCount{$aa}) {
						$aaCount{$aa} = 0;
					}
					$aaCount{$aa}++;
				}			
			}
		}
		
		my $consAa;
		if ($blankCount == $element) {
			$consAa = "-";
		}else {
			my $flag = 0;		
			foreach my $aa (keys %aaCount) {			
				if (!$flag) {
					$consAa = $aa;
					$flag = 1;
				}else {
					if ($aaCount{$aa} > $aaCount{$consAa}) {
						$consAa = $aa;
					}
				}
			}
		}		
		push @consAas, $consAa;
	}
	my $cons = join("", @consAas);
	return $cons;
}

runInSitesOffline( @ARGV );

1;


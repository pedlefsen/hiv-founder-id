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
use Sort::Fields;

use strict;
use warnings;

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

  my @seqs = @$stdNamesRef;
  push @grpSeqs, \@seqs;
  $seqNum = scalar @seqs;
  foreach my $seqName (@seqs) {
     #push @grpSeqs, \@seqs;
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
   #
   # write new sequence file "..._seq4insites.phy"
   #
   WriteFile ($seqFile, \@seqNameNseqs);

   ### This section creates the inSites output files
   ### We are primarily interested in ..._privateSites.txt and
   ### ..._informativeSites.txt -- which have been renamed to Paul's
   ### preferred names.
   ### TAH 10/15

   my $alnDisplayFile = $uploadseqFile.".aln";	            # alignment display 
   my $alnCondenseDisplayFile = $uploadseqFile."_uniq.aln"; # alignment display - condense unique sequence
   my $alnStatFile = $uploadseqFile."_aln.txt";	            # tab delimited alignment summary output file
   my $varSitesAlnFile = $uploadseqFile."_var.aln";	    # aligned variable sites file
   my $varSitesTabFile = $uploadseqFile."_var.txt";	    # tab delimited variable sites file
   my $alnFile = $uploadseqFile."_info.aln";	            # aligned informative output file
   my $privateAlnFile = $uploadseqFile."_priv.aln";	    # aligned private output file
   my $tabFile = $uploadseqFile."_informativeSites.txt";    # tab delimited informative output file
   my $privateTabFile = $uploadseqFile."_privateSites.txt"; # tab delimited private output file
   my $seqArr = my $displaySeqArr = ();
   my $element = my $count = my $maxLen = 0;
   my (@informativeSites, @seqNames, @lines, @alnLines, @privateLines, @privateAlnLines,
       @infoMutant, %infoMutStatus, @privateMutant, %privateMutStatus);
   my (%ambiguityHash, @ambiguousPos);

   foreach my $line (@seqNameNseqs){
      next if $line =~ /^\d+\s+\d+$/;
      my ($seqName, $aaSeq) = split /\t/, $line;
      $aaSeq = uc $aaSeq;
      my $nameLen = length $seqName;
      if ($nameLen > $maxLen) {
         $maxLen = $nameLen;
      }
      $seqNames[$element] = $seqName;	
      my $beginFlag = my $terminalFlag = 0;
      for (my $i = 0; $i < $seqLen; $i++) {
         my $aa = substr($aaSeq, $i, 1);
	 if ($element == 0) {	# this is reference or consensus sequence
	    $seqArr->[$element]->[$i] = $aa;
	 } else {	
            # find ambiguous positions for nucleotide sequence
	    if (!$ambiguityHash{$i}) {
	       if ($datatype eq 'nt' && $aa !~ /[ACGT\-\.\?]/i) {
	          $ambiguityHash{$i} = 1;
		  push @ambiguousPos, $i;
	       }				
	    }
	
	    # deal with leading gaps
	    if ($i == 0 && $aa eq "-") {
	       $beginFlag = 1;
	       $aa = " ";
	    } elsif ($beginFlag == 1 && $aa eq "-") {
	       $aa = " ";
	    } elsif ($aa ne "-") {
	      $beginFlag = 0;
	    }
	    # deal with termianl gaps
	    if (substr ($aaSeq, $i) =~ /^\-+$/) {
	       $terminalFlag = 1;
	    }
	    if ($terminalFlag == 1) {
	       $aa = " ";
	    }
			
	    if ($aa eq $seqArr->[0]->[$i]) {
	       $seqArr->[$element]->[$i] = ".";
	    } else {
	       $seqArr->[$element]->[$i] = $aa;
	    }
	 } # end for non-consensus lines
      } # end for each nucleotide 
      $element++;	# sequence index
   } # end of for each sequence   

   if (@ambiguousPos) {
      my $ambiguousFile = $uploadseqFile."_ambi.txt";	# tab delimited ambiguous output file
      open AMBI, ">$ambiguousFile" or die "couldn't open $ambiguousFile: $!\n";
      foreach my $pos (sort {$a <=> $b} @ambiguousPos) {
         my $site = $pos + 1;
         print AMBI "\t", $site;
      }
      print AMBI "\n";
      for (my $i = 0; $i < $element; $i++) {
         print AMBI $seqNames[$i];
	 foreach my $pos (sort {$a <=> $b} @ambiguousPos) {
	    print AMBI "\t",$seqArr->[$i]->[$pos];
	 }
	 print AMBI "\n";
      }
      close AMBI;
   } # end of if there are ambiguous positions

   my (@consAAs, %infoSiteHash, %privateSiteHash, @privateSites, @uniqNas,
       %uniqNasStatus, %posTotalCount, $posNaCount, @varSites, %varSiteStatus);
   my $infoMutHash = my $privateMutHash = my $grpAaHash = ();
   my $gapOnlyInsiteCount = my $resultFlag = my $gapOnlyPriSiteCount = 0;

   for (my $i = 0; $i < $seqLen; $i++) {
      my $gapOnlyFlag = 1;
      if (!$posTotalCount{$i}) {
         $posTotalCount{$i} = 0;
      }
      for (my $j = 0; $j < $element; $j++) {
         my $aa = $seqArr->[$j]->[$i];
	 if ($j == 0) {
	    $consAAs[$i] = $aa;
	 } else {
	    my $seqName = $seqNames[$j];
	    my $grp = "";
	    if (defined $seqGrp{$seqName}) {	# there is a group file
	       $grp = $seqGrp{$seqName};
	    } else {
	       print "No defined group for sequence $seqName";
	       return 0;
	    }
	    unless ($aa eq " " || $aa eq "?") {	# count
	       my $trueAa = $aa;
	       if ($aa eq ".") {
	          $trueAa = $consAAs[$i];
	       }
	       if (!$uniqNasStatus{$trueAa}) {
	          $uniqNasStatus{$trueAa} = 1;
		  push @uniqNas, $trueAa;					
	       }	
	       $posTotalCount{$i}++;
	       if (!$posNaCount->{$i}->{$trueAa}) {
	          $posNaCount->{$i}->{$trueAa} = 0;
	       }
	       $posNaCount->{$i}->{$trueAa}++;
	    } # end unless $aa (nuc) is a blank or a ?
			
	    unless ($aa eq "." || $aa eq " " || $aa eq "?") {
	       $resultFlag = 1;
	       if (!$grpAaHash->{$i}->{$aa}->{$grp}) {
	          $grpAaHash->{$i}->{$aa}->{$grp} = 0;
	       }
	       $grpAaHash->{$i}->{$aa}->{$grp}++;
	       if (!$varSiteStatus{$i}) {
	          $varSiteStatus{$i} = 1;
	          push @varSites, $i;
	       }
	    } #end unless aa (nuc) is a gap, blank or ?
	 } #end of else
      } #end of for each element (nuc) 
	
      foreach my $aa (keys %{$grpAaHash->{$i}}) {
         foreach my $grp (keys %{$grpAaHash->{$i}->{$aa}}) {
	    if ($grpAaHash->{$i}->{$aa}->{$grp} == 1) {
	       unless ($privateMutHash->{$i}->{$aa}->{$grp}) {
	          $privateMutHash->{$i}->{$aa}->{$grp} = 1;	# label the mutation of the private site
	       }	
				
	       if (!$privateSiteHash{$i}) {
	          $privateSiteHash{$i} = 1;
		  push @privateSites, $i;	# array of private sites
	       }
	    } elsif ($grpAaHash->{$i}->{$aa}->{$grp} > 1) {
	       unless ($infoMutHash->{$i}->{$aa}->{$grp}) {
	          $infoMutHash->{$i}->{$aa}->{$grp} = 1;	# label the mutation of the informative site
	       }					
	
	       if (!$infoSiteHash{$i}) {
	          $infoSiteHash{$i} = 1;
	       push @informativeSites, $i;	# array of informative sites
	       }
	    } else {
	       die "Something wrong here!\n";
	    }
	 } #end of for each group hash
      } #end of for each aa (nuc) in a group

      if ($infoSiteHash{$i}) {	# informative site
         if (keys %{$infoMutHash->{$i}} == 1) { # only one mutation other than consensus
            my ($mut) = keys %{$infoMutHash->{$i}};
	    if ($consAAs[$i] eq "-") {
	       $gapOnlyInsiteCount++;
	    } elsif ($mut eq "-") {
	       $gapOnlyInsiteCount++;
	    }
	 }
      }	
	
      if ($privateSiteHash{$i}) {	# private site
         if (keys %{$privateMutHash->{$i}} == 1) { # only one mutation other than consensus
	    my ($mut) = keys %{$privateMutHash->{$i}};
	    if ($consAAs[$i] eq "-") {
	       $gapOnlyPriSiteCount++;
	    } elsif ($mut eq "-") {
	       $gapOnlyPriSiteCount++;
	    }
	 }
      }	
   } #end of for each aa (nuc) in alignment

   ### Output file section

   if (!$resultFlag) {
      print "There are no informative nor private sites in the sequence alignment.\n";
   } else {
      # write to alignment display file
      WriteAlnDisplay ($alnDisplayFile, \@seqNames, $seqArr, $element, $seqLen, $maxLen, $alnCondenseDisplayFile);
      $maxLen += 6;
      if (@varSites) {		
         WriteAlnVarSites ($varSitesAlnFile, $maxLen, \@varSites, \@seqNames, $seqArr, $element,  $seqLen, \%seqGrp);
         WriteTabVarSites ($varSitesTabFile, \@varSites, \@seqNames, $seqArr, $element, \%seqGrp);
      }
	
      if (@informativeSites) {
         my $param =
            GetParams ($maxLen, $element, \@seqNames, \@informativeSites, $seqArr,
                       \%infoMutStatus, \@infoMutant, $infoMutHash, \%seqGrp, $datatype);
	 WriteAlnInsites ($alnFile, $maxLen, \@informativeSites, $param, $sortRadio, $seqLen);
	 my $TAB;
	 open ($TAB, ">$tabFile") or die "Couldn't open $tabFile: $!\n";
	 print $TAB "\tTotal_informative\tInformative(noGaps)\tAmbiguities";
         WriteTabInsites ($TAB, \@informativeSites, $gapOnlyInsiteCount, $param, $sortRadio,
                          $seqArr, $element, $datatype, \@nas, \@infoMutant, \%infoMutStatus, $infoMutHash);
	 close $TAB;		
      }
	
      if (@privateSites) {
         my $param =
            GetParams ($maxLen, $element, \@seqNames, \@privateSites, $seqArr,
                       \%privateMutStatus, \@privateMutant, $privateMutHash, \%seqGrp, $datatype);
	 WriteAlnInsites ($privateAlnFile, $maxLen, \@privateSites, $param, $sortRadio, $seqLen);
	 my $TAB;
	 open ($TAB, ">$privateTabFile") or die "Couldn't open $privateTabFile: $!\n";
	 print $TAB "\tTotal_private\tPrivate(noGaps)\tAmbiguities";
	 WriteTabInsites ($TAB, \@privateSites, $gapOnlyPriSiteCount, $param, $sortRadio, $seqArr,
                          $element, $datatype, \@nas, \@privateMutant, \%privateMutStatus, $privateMutHash);
	 close $TAB;		
      }
      
      # write to alignment summary file	
      open STAT, ">$alnStatFile" or die "couldn't open $alnStatFile: $!\n";
      print STAT "Position";
      foreach my $na (sort @uniqNas) {
         print STAT "\t$na";
      }
      print STAT "\tTotal";
      foreach my $na (sort @uniqNas) {
         print STAT "\t$na";
      }
      print STAT "\t1st Freq.\tFreq.\t2nd Freq.\tFreq.\t3rd Freq.\tFreq.\t4th Freq.\tFreq.\n";
      for (my $i = 0; $i < $seqLen; $i++) {
         my $pos = $i + 1;
	 print STAT $pos;
	 foreach my $na (sort @uniqNas) {
            if ($posNaCount->{$i}->{$na}) {
	       print STAT "\t", $posNaCount->{$i}->{$na};
	    } else {
	       print STAT "\t0";
	    }
	 }
	 if ($posTotalCount{$i}) {
	    print STAT "\t", $posTotalCount{$i};
	 } else {
	    print STAT "\t0";
	 }
	 foreach my $na (sort @uniqNas) {
	    if ($posNaCount->{$i}->{$na}) {
	       my $freq = int ($posNaCount->{$i}->{$na} / $posTotalCount{$i} * 10000 + 0.5) / 10000;
	       print STAT "\t", $freq;
	    } else {
	       print STAT "\t0.0000";
	    }
	 }
	 my $idx = 0;
	 foreach my $na (sort {$posNaCount->{$i}->{$b} <=> $posNaCount->{$i}->{$a}} keys %{$posNaCount->{$i}}) {
	    $idx++;
	    last if ($idx > 4);
	    print STAT "\t$na";
	    my $freq = int ($posNaCount->{$i}->{$na} / $posTotalCount{$i} * 10000 + 0.5) / 10000;
	    print STAT "\t", $freq;
	 }
	 print STAT "\n";
      } # end of for-each-position statistics file
      print STAT "\n";
      close STAT;
   }  ### end of IF there is information to print

   ### Paul's code, sans CGI stuff

   my ( $no_informative_sites ) = (!@informativeSites);
   my ( $job_id ) = ( $id );

   if( $VERBOSE ) {
      print "JOB ID IS $job_id\n";
      if( $no_informative_sites ) {
         print "NO INFORMATIVE SITES\n";
      }
   }

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

## Ripped off verbatim from Insites.pm
## To Do:  remove unused portions, e.g. protein stuff.
## maybe normalize spacing.
## TAH 10/15
sub WriteAlnDisplay {
	my ($alignDisplayFile, $seqNames, $seqArr, $element, $seqLen, $maxLen, $uniqAlignDisplayFile) = @_;
	my %seqCount = my %seqStatus = ();
	$maxLen += 10;
	open ALNDISPLAY, ">$alignDisplayFile" or die "couldn't open $alignDisplayFile: $!\n";
	for (my $i = 0; $i < $element; $i++) {
		my $seq = join ("", @{$seqArr->[$i]});
		$seq =~ s/ /\-/g;		
		printf ALNDISPLAY "%-".$maxLen."s", $seqNames->[$i];
		print ALNDISPLAY "$seq\n";
	}
	close ALNDISPLAY;
	
	open UNIQDISPLAY, ">$uniqAlignDisplayFile" or die "couldn't open $uniqAlignDisplayFile: $!\n";	
	for (my $i = 0; $i < $element; $i++) {
		my $seq = join ("", @{$seqArr->[$i]});
		$seq =~ s/ /\-/g;
		if ($i == 0) {
			printf UNIQDISPLAY "%-".$maxLen."s", $seqNames->[$i];
			print UNIQDISPLAY "$seq\n";
		}else {
			if (!$seqStatus{$seq}) {
				$seqStatus{$seq} = $seqNames->[$i];
			}
			$seqCount{$seq}++;
		}
	}
	my $idx = 0;
	foreach my $seq (sort {$seqCount{$b} <=> $seqCount{$a}} keys %seqCount) {
		$idx++;
		my $name = $seqStatus{$seq}."_".$seqCount{$seq};
		printf UNIQDISPLAY "%-".$maxLen."s", $name;
		print UNIQDISPLAY "$seq\n";
	}
	close UNIQDISPLAY;
}

## Ripped off verbatim from Insites.pm
## To Do:  remove unused portions, e.g. protein stuff.
## maybe normalize spacing.
## TAH 10/15
sub WriteAlnVarSites {
	my ($varSitesAlnFile, $maxLen, $varSitesRef, $seqNamesRef, $seqArr, $element,  $seqLen, $seqGrp) = @_;
	
	my $digits = 1;
	while (1) {
		if ($seqLen =~ /^\d{$digits}$/) {
			last;
		}
		$digits++
	}
	
	open (ALN, ">$varSitesAlnFile") or die "Couldn't open $varSitesAlnFile: $!\n";
	
	printf ALN "%".$maxLen."s", "";
	
	for (my $i = $digits; $i > 0; $i--) {
		foreach my $position (@$varSitesRef) {
			my $aaIndex = $position + 1;
			my $pattern = "(\\d)";
			for (my $j = 1; $j < $i; $j++) {
				$pattern .= "\\d";		
			}
			if ($aaIndex =~ /$pattern$/) {				
				print ALN $1;
			}else {
				print ALN " ";
			}
		}
		print ALN "\n";
		if ($i == 1) {
			print ALN "\n";
		}else {
			printf ALN "%".$maxLen."s", "";
		}	
	}
	my $grp = "";
	for (my $i = 0; $i < $element; $i++) {
		if ($i == 0) {
			printf ALN "%-".$maxLen."s", $seqNamesRef->[$i];
			foreach my $varSite (@{$varSitesRef}) {
				print ALN $seqArr->[$i]->[$varSite];
			}
			print ALN "\n";
		}else {
			if (!$grp || $grp && ($seqGrp->{$seqNamesRef->[$i]} ne $grp)) {
				$grp = $seqGrp->{$seqNamesRef->[$i]};
				print ALN "\n";
			}
			printf ALN "%-".$maxLen."s", $seqNamesRef->[$i];
			foreach my $varSite (@{$varSitesRef}) {
				print ALN $seqArr->[$i]->[$varSite];
			}
			print ALN "\n";
		}		
	}
	close ALN;
}

## Ripped off verbatim from Insites.pm
## To Do:  remove unused portions, e.g. protein stuff.
## maybe normalize spacing.
## TAH 10/15
sub WriteTabVarSites {
	my ($varSitesTabFile, $varSitesRef, $seqNamesRef, $seqArr, $element, $seqGrp) = @_;
	open TAB, ">$varSitesTabFile" or die "couldn't open $varSitesTabFile: $!\n";
	foreach my $site (@$varSitesRef) {
		my $pos = $site + 1;
		print TAB "\t$pos";
	}
	print TAB "\n";
	my $grp = "";
	for (my $i = 0; $i < $element; $i++) {
		if ($i == 0) {
			print TAB $seqNamesRef->[$i];
			foreach my $varSite (@{$varSitesRef}) {
				print TAB "\t$seqArr->[$i]->[$varSite]";
			}
			print TAB "\n";
		}else {
			if (!$grp || $grp && ($seqGrp->{$seqNamesRef->[$i]} ne $grp)) {
				$grp = $seqGrp->{$seqNamesRef->[$i]};
				print TAB "\n";
			}
			print TAB $seqNamesRef->[$i];
			foreach my $varSite (@{$varSitesRef}) {
				print TAB "\t$seqArr->[$i]->[$varSite]";
			}
			print TAB "\n";
		}		
		
	}
	close TAB;
}

## Ripped off verbatim from Insites.pm
## To Do:  remove unused portions, e.g. protein stuff.
## maybe normalize spacing.
## TAH 10/15
sub GetParams {
	my ($maxLen, $element, $seqnamesRef, $sitesRef, $seqArr, $mutStatusRef, $mutantRef, $aaHash, $seqGrp, $datatype) = @_;
	my (@alnLines, @lines, $param);
	my $grp = "";
	for (my $i = 0; $i < $element; $i++) {
		my $line = "";
		my $alnLine = my $seqName = $seqnamesRef->[$i];
		my $len = length $seqnamesRef->[$i];
				
		for (my $j = $len; $j < $maxLen; $j++) {
			$alnLine .= " ";
		}
		my $infoCount = my $gapInfoCount = my $ambiCount = 0;
		foreach my $position (@$sitesRef) {
			my $aa = $seqArr->[$i]->[$position];
			if ($i == 0) {
				if (!$mutStatusRef->{$aa}) {
					$mutStatusRef->{$aa} = 1;
					push @$mutantRef, $aa;
				}
			}else {
				if (!$grp || ($grp && $seqGrp->{$seqName} ne $grp)) {
					$grp = $seqGrp->{$seqName};
					push @alnLines, "";
					push @lines, "";
				}
				unless ($aa eq "." || $aa eq " " || $aa eq "?") {
					if (!$mutStatusRef->{$aa}) {
						$mutStatusRef->{$aa} = 1;
						push @$mutantRef, $aa;
					}
					if (!$aaHash->{$position}->{$aa}->{$grp}) {	# private site, change aa to consensus for display
						#$aa = ".";
					}else {					
						$infoCount++;
						if ($aa eq "-") {
							$gapInfoCount++;
						}else {
							if ($datatype eq 'nt' && $aa !~ /[ACGTacgt\-]/) {
								$ambiCount++;
							}elsif ($datatype eq 'aa' && $aa eq 'X') {
								$ambiCount++;
							}
						}
					}
				}
			}
			$alnLine .= $aa;
			$line .= "\t".$aa;
		}
		
		push @alnLines, $alnLine;
		my $noGapInfoCount = $infoCount - $gapInfoCount;
		my $tabLine;
		if ($line) {
			$tabLine = $seqnamesRef->[$i]."\t".$infoCount."\t".$noGapInfoCount."\t".$ambiCount.$line;
		}
		push @lines, $tabLine;
	}
	$param->{lines} = \@lines;
	$param->{alnLines} = \@alnLines;
	return $param;
}

## Ripped off verbatim from Insites.pm
## To Do:  remove unused portions, e.g. protein stuff.
## maybe normalize spacing.
## TAH 10/15
sub WriteAlnInsites {
	my ($alnFile, $maxLen, $informativeSitesRef, $param, $sortRadio, $nchar) = @_;
	my $firstAlnLine = shift @{$param->{alnLines}};
	my @alnOutputs = @{$param->{alnLines}};
	
	my $digits = 1;
	while (1) {
		if ($nchar =~ /^\d{$digits}$/) {
			last;
		}
		$digits++
	}
	
	open (ALN, ">$alnFile") or die "Couldn't open $alnFile: $!\n";
	
	printf ALN "%".$maxLen."s", "";
	
	for (my $i = $digits; $i > 0; $i--) {
		foreach my $position (@$informativeSitesRef) {
			my $aaIndex = $position + 1;
			my $pattern = "(\\d)";
			for (my $j = 1; $j < $i; $j++) {
				$pattern .= "\\d";		
			}
			if ($aaIndex =~ /$pattern$/) {				
				print ALN $1;
			}else {
				print ALN " ";
			}
		}
		print ALN "\n";
		if ($i == 1) {
			print ALN "\n";
		}else {
			printf ALN "%".$maxLen."s", "";
		}	
	}
	
	print ALN $firstAlnLine,"\n";
	
	if ($sortRadio eq "y") {
		my $alnSortFields;
		
		for (my $i = $maxLen+1; $i <= @$informativeSitesRef+$maxLen; $i++) {
			push @$alnSortFields, $i;
		}
		my @alnSorted = fieldsort '', $alnSortFields, @{$param->{alnLines}};
		@alnOutputs = reverse @alnSorted;
	}
	
	foreach my $line (@alnOutputs) {
		$line =~ s/\t//g;
		print ALN $line,"\n";
	}
	close ALN;
}

## Ripped off verbatim from Insites.pm
## To Do:  remove unused portions, e.g. protein stuff.
## maybe normalize spacing.
## TAH 10/15
sub WriteTabInsites {
	my ($TAB, $informativeSitesRef, $gapOnlySiteCount, $param, $sortRadio, $seqArr, $element, $datatype, $nasRef, $infoMutantRef, $infoMutStatusRef, $posMutHash) = @_;
	my $informativeCount = scalar @$informativeSitesRef;
	my $notGapOnlyInSitesCount = $informativeCount - $gapOnlySiteCount;
	my $alignAmbiCount = 0;
	foreach my $site (@$informativeSitesRef) {
		my $position = $site + 1;
		print $TAB "\t", $position;
		foreach my $na (keys %{$posMutHash->{$site}}) {
			if ($datatype eq 'nt' && $na !~ /[ACGTacgt\-]/) {
				$alignAmbiCount++;
				last;
			}elsif ($datatype eq 'aa' && $na eq 'X') {
				$alignAmbiCount++;
				last;
			}
		}
	}
	print $TAB "\n";
		
	my $firstLine = shift @{$param->{lines}};
	print $TAB $firstLine,"\n";
	my @tabOutputs = @{$param->{lines}}; 
	
	if ($sortRadio eq "y") {
		my $tabSortFields;
		
		for (my $i = 5; $i <= $informativeCount+4; $i++) {
			push @$tabSortFields, $i;
		}	
		my @tabSorted = fieldsort '\t', $tabSortFields, @{$param->{lines}};
		@tabOutputs = reverse @tabSorted;
	}
	
	foreach my $line (@tabOutputs) {
		print $TAB $line,"\n";
	}
	
	print $TAB "\nAlignment\t$informativeCount\t$notGapOnlyInSitesCount\t$alignAmbiCount\n";
	
	# calculate the total na/aa at each informative site
	my (%totalNAcount, $naCount, $infoSiteMutation, $totalMutation, %mutationStatus, @mutations);
	foreach my $position (@$informativeSitesRef) {
		my $count = 0;
		my $consAa = $seqArr->[0]->[$position];
		for (my $i = 1; $i < $element; $i++) {	# exclude CON
			my $aa = $seqArr->[$i]->[$position];
	
			unless ($aa eq " ") {
				unless ($aa eq "?") {
					if ($aa eq ".") {
						$aa = $consAa;
					}
					if (!$naCount->{$position}->{$aa}) {
						$naCount->{$position}->{$aa} = 0;
					}
					$naCount->{$position}->{$aa}++;
				}			
				$count++;
			}
		}
		$totalNAcount{$position} = $count;
		
		# calculate mutation information (will exclude gap to na and na to gap)
		unless ($consAa eq "-") {
			foreach my $na (sort keys %{$posMutHash->{$position}}) {
				unless ($na eq $consAa || $na eq "-") {
					my $mutation = $consAa."-".$na;
					$infoSiteMutation->{$position}->{$mutation} = 1;
					if (!$totalMutation->{$mutation}) {
						$totalMutation->{$mutation} = 0;
					}
					$totalMutation->{$mutation}++;
					if (!$mutationStatus{$mutation}) {
						$mutationStatus{$mutation} = 1;
						push @mutations, $mutation;
					}
				}
			}
		}
	}
	
	foreach my $na (sort @$infoMutantRef) {
		print $TAB "\t\t\t$na";
		foreach my $position (@$informativeSitesRef) {
			my $count = 0;
			if ($naCount->{$position}->{$na}) {
				$count = $naCount->{$position}->{$na};
			}				
			print $TAB "\t",$count;
		}
		print $TAB "\n";
	}
	
	foreach my $an (@$nasRef) {
		unless ($infoMutStatusRef->{$an}) {
			print $TAB "\t\t\t$an";
			foreach my $position (@$informativeSitesRef) {
				print $TAB "\t0";
			}
			print $TAB "\n";
		}
	}
	
	print $TAB "\t\t\tTotal";
	foreach my $position (@$informativeSitesRef) {
		print $TAB "\t",$totalNAcount{$position};
	}
	print $TAB "\n\n";
	
	# write mutation information
	print $TAB "\t\t\tTypes of Mutation";
	foreach (@$informativeSitesRef) {
		print $TAB "\t";
	}
	print $TAB "\tTotal\n";

	# for all possible mutaion types including ambiguities
	my (@mutation_types, @all_nas);
	if ($datatype eq "nt") {		
		my @ambi_nas = qw (R Y K M S W B D H V N);
		push @all_nas, @$nasRef, @ambi_nas;
	}else {
		push @all_nas, @$nasRef, "X";
	}
	for (my $i = 0; $i < @all_nas; $i++) {
		for (my $j = 0; $j < @all_nas; $j++) {
			unless ($i == $j) {
				my $mut_type = $all_nas[$i]."-".$all_nas[$j];
				push @mutation_types, $mut_type;
			}
		}
	}
	foreach my $mutation (@mutation_types) {
		print $TAB "\t\t\t$mutation";
		foreach my $position (@$informativeSitesRef) {
			my $count = 0;
			if ($infoSiteMutation->{$position}->{$mutation}) {
				$count = $infoSiteMutation->{$position}->{$mutation};
			}
			print $TAB "\t",$count;
		}
		if ($totalMutation->{$mutation}) {
			print $TAB "\t", $totalMutation->{$mutation}, "\n";
		}else {
			print $TAB "\t0\n";
		}		
	}	
}


### main routine
runInSitesOffline( @ARGV );

1;


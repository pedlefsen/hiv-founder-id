#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) runPhyML
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##  Description:
##
##      Script for running PhyML locally to get the pairwise
##      diversities and the newick-formatted estimated tree.
##
##      PhyML is available on github at
##      https://github.com/stephaneguindon/phyml
##      and executables can be downloaded from
##      http://www.atgc-montpellier.fr/phyml/versions.php
##
##      Modified from diver.pl (with code from Diver.pm), all by
##      Wenjie Deng, version as of September 17, 2015.
##
##      Note that a patch is needed to get phyml to print the distance matrix. 
##      See https://groups.google.com/forum/#!topic/phyml-forum/yeIwFkwsI8c
##      I ended up having problems with the phyml on github so I got it from
##      https://github.com/stephaneguindon/phyml-downloads/releases/download/stable/phyml-20120412.tar.gz
##      
###******************************************************************************

use Getopt::Std; # for getopts
use File::Path qw( make_path );
use Path::Tiny;
use Statistics::Descriptive;
use Readonly;

# use Bio::TreeIO;
# use Bio::Tree::TreeFunctionsI;
# use Net::SMTP;

use strict;
use vars qw( $opt_D $opt_V );
use vars qw( $VERBOSE $DEBUG );

## NOTE THAT THIS NEEDS TO BE PATCHED; SEE ABOVE.
use constant PHYML_EXECUTABLE => "./phyml";

Readonly my %PHYML_OPTIONS => (
  datatype       => "nt",	# sequence data type ('nt' or 'aa')
  bootstrap      => 0,	        # 0 means neither aLRT nor bootstrap
  subModel       => "GTR",	# substitution model
  freqRadio      => "m",	# equilibrium frequencies
  ttRatio        => "e",	# transition/transversion ratio
  proportion     => "e",	# proportion of invariable sites
  subRateCat     => 4,	# number of substitution rate categories
  gamma          => "e",	# Gamma distribution parameter
  treeImprovType => "NNI",	# type of tree improvement
  opt            => "tlr"	# optimise topology, branch lengths or rate parameters
        );
  
sub runPhyML {
  @ARGV = @_;

  sub runPhyML_usage {
    print "\trunPhyML [-DV] <input_fasta_file> [<output_dir>]\n";
    exit;
  }

  # This means -D and -V are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # But first reset the opt vars.
  ( $opt_D, $opt_V ) = ();
  if( not getopts('DV') ) {
    runPhyML_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;
  
  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }
  
  my $input_fasta_file = shift @ARGV || runPhyML_usage();
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
  make_path( $output_path_dir );

  my $phylipFile = $output_path_dir . '/' . $input_fasta_file_short_nosuffix .'.phylip';	# reformatted for phyml
  my $phymlOutFile = $phylipFile . "_phyml.out"; # phyml's redirected standard output with distance matrix (phyml patched to show it)
  my $errFile = $phylipFile . "_phyml.err"; # pairwise (not tree-based) diversities in a table format.
  my $pwDiversityFile = $phylipFile . "_phyml_pwdiversity.txt"; # pairwise (not tree-based) diversities in a table format.

  ## First, using code from Wenjie's Fasta2Phylip.pl script (from http://indra.mullins.microbiol.washington.edu/perlscript/docs/Sequence.html)

  # Create a temporary file for the conversion via "unix", whatever Wenjie meant by that.
  my $unixFile = $phylipFile .".unix";
  ConvertToUnix( $input_fasta_file, $unixFile );
  my @seqNames = ChangetoPhylip( $unixFile, $phylipFile );
  my $seqCount = scalar( @seqNames );
  unlink( $unixFile );# Removing the temporary file.

  my @cmd = (
              PHYML_EXECUTABLE,
             '-i', $phylipFile,
             '-d', $PHYML_OPTIONS{ "datatype" },
             '-q', # indicates that the file is in a sequential (not interleaved) format.
             '-b', $PHYML_OPTIONS{ "bootstrap" },
             '-m', $PHYML_OPTIONS{ "subModel" },
             '-v', $PHYML_OPTIONS{ "proportion" },
             '-c', $PHYML_OPTIONS{ "subRateCat" },
             '-o', $PHYML_OPTIONS{ "opt" },
             '-a', $PHYML_OPTIONS{ "gamma" },
             '-f', $PHYML_OPTIONS{ "freqRadio" },
             '-t', $PHYML_OPTIONS{ "ttRatio" },
             '-s', $PHYML_OPTIONS{ "treeImprovType" }
  );
  
  # Run PhyML
  if( $VERBOSE ) {
    print( "Running PhyML using command:\n" );
    print( join( ' ', @cmd ), "\n" );
  }
  # run it; redirect STDOUT and STDERR
  my $cmd_str  = join( ' ', @cmd );
  my $phyml_result = system( "$cmd_str 1>$phymlOutFile 2>$errFile" );
  if( $VERBOSE ) {
    print( "Done running PhyML." );
  }
  if( -s $errFile ) {
    open ERR, $errFile;
    my @lines = <ERR>;
    my $errMsg = join('', @lines);
    close ERR;
    die( "Error running PhyML: $errMsg" );
  }

  # Get pairwise distances.
  my $pwDistHashRef;
  my( @alternative_seqNames );
  my( @seqNamesWithDists );
  my $flag = my $count = 0;
  if( $VERBOSE ) {
    print( "Extracting pairwise distances from PhyML output in file $phymlOutFile\n" );
  }
  if( $VERBOSE ) {
    print( "Opening file $phymlOutFile for reading.." );
  }
  open( IN, $phymlOutFile ) or die( "Couldn't open $phymlOutFile: $!\n" );
  if( $VERBOSE ) {
    print( ".done\n" );
  }
  if( $VERBOSE ) {
    print( "Extracting pairwise distances from PhyML output.\n" );
  }
  while( my $line = <IN> ) {
    chomp $line;
    next if ($line =~ /^\s*$/);
    if( $line =~ /Building BioNJ tree/ ) {
      $flag = 1;
      next;
    }
    if( $flag && $line =~ /^(\d+)$/ ) {
      $count = $1;
      next;
    }
    if( $line =~ /^\./ ) {
      $flag = 0;
    }
    if( $flag ) {
      if( $line =~ /^(.+)$/ ) {
        my $seqnName = $1;				
        push @seqNamesWithDists, $seqnName;
        $seqnName =~ /^(\S+)\s+/;
        my $seqName = $1;
        unless( $seqName eq $seqNames[ scalar( @alternative_seqNames ) ] ) {
          warn( "seqName changed: was " . $seqNames[ scalar( @alternative_seqNames ) ] . "; is now " . $seqName );
          $seqNames[ scalar( @alternative_seqNames ) ] = $seqName;
        }
        push @alternative_seqNames, $seqName;
      }
    }
  }
  if( $VERBOSE ) {
    print( "Done.\nClosing file $phymlOutFile.." );
  }
  close IN;
  if( $VERBOSE ) {
    print( ".done\n" );
  }
  foreach my $seqNameWithDists ( @seqNamesWithDists ) {
    $seqNameWithDists =~ /^(\S+)\s+(\S+)(.*)$/;
    my $name = $1;
    my $dists = $2.$3;	
    my @distances = split /\s+/, $dists;
    for( my $i = 0; $i < scalar( @distances ); $i++ ) {
      unless( $name eq $seqNames[ $i ] ) {					
        if( $distances[ $i ] eq '-' ) { # This means the sequences are nonoverlapping (eg right half and left half)
          next; # We just don't include these in the calculation.  PhyML pairwise diffs exclude gaps.
        }
        $pwDistHashRef->{ $name }->{ $seqNames[ $i ] } =
          $distances[ $i ];
      }
    }		
  }
  
  my @diversity;
  if( scalar @seqNames == 1 ) {	# for the case of only one sequence in the group
    push @diversity, 0;
  } else {	# at least two sequences in the group
    my ( $firstName, $secondName, $i, $j );
    for( $i = 0; $i < scalar( @seqNames ); $i++ ) {	# print distance between same group
      $firstName = $seqNames[ $i ];
      if( !defined( $pwDistHashRef->{ $firstName } ) ) {
        print "Warning: no distance values from $firstName\n";
      } else {
        for( $j = ( $i + 1 ); $j < scalar( @seqNames ); $j++ ) {
          $secondName = $seqNames[ $j ];
          if( $firstName eq $secondName ) {
            next;
          }
          if( !defined( $pwDistHashRef->{ $firstName }->{ $secondName } ) ) {
            print "Warning: no distance value between $firstName and $secondName\n";
          } else {
            push @diversity, $pwDistHashRef->{ $firstName }->{ $secondName };
            print $firstName, " ", $secondName, " ", $pwDistHashRef->{ $firstName }->{ $secondName }, "\n";
          }
        } # End foreach $j
      } # End if $firstName is in the hash
    } # End foreach $i
  } # Else if at least two sequences in the group

  print( "There are ", scalar( @diversity ), " pairwise distances.\n" );
		
  if( $VERBOSE ) {
    print( "Opening file $pwDiversityFile for writing.." );
  }
  open OUT, ">$pwDiversityFile" || die "Couldn't open diversity output file $pwDiversityFile: $!\n";
  if( $VERBOSE ) {
    print( ".done\n" );
  }
  if( $VERBOSE ) {
    print( "Writing diversity statistics to file $pwDiversityFile.." );
  }
  print OUT "Number of sequences\tMean\tStandard error\tMin\tQ1\tMedian\tQ3\tMax\n";
  if( scalar( @diversity ) > 0 ) {
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data( @diversity );
    my $mean = $stat->mean();
    my $sqrt_n = ( scalar( @diversity ) ** 0.5 );
    my $sem = ( $stat->standard_deviation() / $sqrt_n );
    my $min = $stat->min();
    my $q1 = $stat->quantile( 1 );
    my $median = $stat->quantile( 2 );
    my $q3 = $stat->quantile( 3 );
    my $max = $stat->max();
    my @diversityResult = ( $mean, $sem, $min, $q1, $median, $q3, $max );
  
    print OUT scalar( @seqNames ), "\t"; 
    print OUT join( "\t", @diversityResult ), "\n";
  } else {
    print OUT "0\t0\t0\t0\t0\t0\t0\t0\n";
  }

  if( $VERBOSE ) {
    print( "Done.\nClosing file $phymlOutFile.." );
  }
  close OUT;
  if( $VERBOSE ) {
    print( ".done\n" );
  }

  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  return 0;
} # runPhyML(..)

######################################################################################
###### From Wenjie Deng's Fasta2Phylip.pl script (from http://indra.mullins.microbiol.washington.edu/perlscript/docs/Sequence.html)
######
sub ConvertToUnix {
	my ($infile, $unixFile) = @_;
	open (IN, $infile) or die "Couldn't open $infile: $!\n";
	open (OUT, ">$unixFile") or die "Couldn't open $unixFile: $!\n";
	my @buffer = <IN>;
	close IN;
	my $line = "";
	foreach my $element (@buffer) {
		$line .= $element;
	}
	if ($line =~ /\r\n/) {
		$line =~ s/\r//g;
	}elsif ($line =~ /\r/) {
		$line =~ s/\r/\n/g;
	}
	print OUT $line;	
	close OUT;	
}
sub ChangetoPhylip {
	my ($unixFile, $phylipFile) = @_;
	my $seqCount = 0;
	my $seq = my $seqName = "";
	open IN, $unixFile or die "Couldn't open $unixFile\n";
	while (my $line = <IN>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		if ($line =~ /^>/) {
			$seqCount++;
		}elsif ($seqCount == 1) {
			$seq .= $line;
		}
	}
	close IN;
	my $seqLen = length $seq;
	
	open(IN, $unixFile) || die "Can't open $unixFile\n";
	open(OUT, ">$phylipFile") || die "Cant open $phylipFile\n";
	print OUT $seqCount," ",$seqLen,"\n";
	$seqCount = 0;
	$seq = "";
        my @seqNames;
	while(my $line = <IN>) {
		chomp $line;	
		next if($line =~ /^\s*$/);
	
		if($line =~ /^>(\S+)/) {
			if ($seqCount) {
				my $len = length $seq;
				if ($len == $seqLen) {
					print OUT "$seqName\t$seq\n";
					$seq = $seqName = "";
				}else {
					unlink $unixFile;
					unlink $phylipFile;
					die "Error: the sequence length of $seqName is not same as others.\n";
				}
			}	
			$seqName = $seqCount; # PAUL CHANGED; USING ONLY THE index is safest
                        push @seqNames, $seqName;
			$seqCount++;
		}else {
			$seq .= $line;		
		}		
	}
	close IN;
	# check the length of last sequence
	my $len = length $seq;
	if ($len == $seqLen) {
		print OUT "$seqName\t$seq\n";
	}else {
		unlink $unixFile;
		unlink $phylipFile;
		die "Error: the sequence length of $seqName is not same as others.\n";
	}	
	close IN;
	close OUT;
        ## PAUL CHANGED; ADDED RETURN STATEMENT.
        die unless( scalar( @seqNames ) == $seqCount );
        return( @seqNames );
}
######################################################################################

runPhyML( @ARGV );

1;


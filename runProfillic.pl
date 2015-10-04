#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) runPhyML
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##  Description:
##
##      Script for running profillic locally to get the profile estimated
##      and alignment profiles generated.
##
##      profillic and profuse are available on github at
##      https://github.com/galosh/profillic
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

# For now I'm not actually running the trainer, just using the given alignment to make a profile (hmmer-style).
#use constant PROFILLIC_EXECUTABLE => "profillic/dist/profillic_DNA_CQA";
use constant ALIGNEDFASTATOPROFILE_EXECUTABLE => "profuse/dist/alignedFastaToProfile_DNA";
use constant PROFILETOALIGNMENTPROFILE_EXECUTABLE => "profuse/dist/profileToAlignmentProfile_AA";

## ERE I AM HAVE TO SET IT UP NOT WORRYING ABOUT PROFILLIC (training) JUST YET
Readonly my %PROFILLIC_OPTIONS => (
  random_seed => 98103,	                # use my zip code as a default.  Not much worse than any other.
  even_starting_profile_multiple => -1	# start at the starting profile.
        );
  
sub runProfillic {
  @ARGV = @_;

  sub runProfillic_usage {
    print "\trunProfillic [-DV] <input_fasta_file> [<output_dir>]\n";
    exit;
  }

  # This means -D and -V are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # But first reset the opt vars.
  ( $opt_D, $opt_V ) = ();
  if( not getopts('DV') ) {
    runProfillic_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;
  
  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }
  
  my $input_fasta_file = shift @ARGV || runProfillic_usage();
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

  ## We have to first create an ungapped version.  Here we are lazily not protecting the sequence names from also losing dashes.
  my $input_fasta_file_contents = path( $input_fasta_file )->slurp();
  my $ungapped_fasta_file_contents = ( $input_fasta_file_contents =~ s/-//g );
  # Now write it out to a temporary location in the output dir.
  my $ungapped_fasta_file = "${output_path_dir_for_input_fasta_file}/${input_fasta_file_short_nosuffix}_ungapped${input_fasta_file_suffix}";
  if( $VERBOSE ) {
    print( "Writing out ungapped input file \"$ungapped_fasta_file\".." );
  }
  if( $VERBOSE ) { print "Opening file \"$ungapped_fasta_file for writing..\n"; }
  unless( open ungapped_fasta_fileFH, ">$ungapped_fasta_file" ) {
    warn "Unable to open output file \"$ungapped_fasta_file: $!\n";
    return 1;
    }
  print ungapped_fasta_fileFH $ungapped_fasta_file_contents;
  close( ungapped_fasta_fileFH );

  my @profillic_cmd = (
             PROFILLIC_EXECUTABLE,
             $ungapped_fasta_file,
             '--starting_profile', $starting_profile, # ERE I AM! this one requries running alignedFastaToProfile_DNA, which requires having a consensus.  I could get it in perl I know.
             '--even_starting_profile_multiple', $PROFILLIC_OPTIONS{ "even_starting_profile_multiple" }
             '-b', $PROFILLIC_OPTIONS{ "bootstrap" },
             '-m', $PROFILLIC_OPTIONS{ "subModel" },
             '-v', $PROFILLIC_OPTIONS{ "proportion" },
             '-c', $PROFILLIC_OPTIONS{ "subRateCat" },
             '-o', $PROFILLIC_OPTIONS{ "opt" },
             '-a', $PROFILLIC_OPTIONS{ "gamma" },
             '-f', $PROFILLIC_OPTIONS{ "freqRadio" },
             '-t', $PROFILLIC_OPTIONS{ "ttRatio" },
             '-s', $PROFILLIC_OPTIONS{ "treeImprovType" }
  );
  
  # Run Profillic
  if( $VERBOSE ) {
    print( "Running Profillic using command:\n" );
    print( join( ' ', @cmd ), "\n" );
  }
  # run it; redirect STDOUT and STDERR
  system( join( ' ', @cmd ) . "1>$profillicOutFile 2>$errFile" );
  if( $VERBOSE ) {
    print( "Done running Profillic." );
  }
  if( -s $errFile ) {
    open ERR, $errFile;
    my @lines = <ERR>;
    my $errMsg = join('', @lines);
    close ERR;
    die( "Error running Profillic: $errMsg" );
  }


  
  my @profillic_cmd = (
             PROFILLIC_EXECUTABLE,
             $ungapped_fasta_file,
             '--starting_profile', $starting_profile, # ERE I AM! this one requries running alignedFastaToProfile_DNA, which requires having a consensus.  I could get it in perl I know.
             '--even_starting_profile_multiple', $PROFILLIC_OPTIONS{ "even_starting_profile_multiple" }
             '-b', $PROFILLIC_OPTIONS{ "bootstrap" },
             '-m', $PROFILLIC_OPTIONS{ "subModel" },
             '-v', $PROFILLIC_OPTIONS{ "proportion" },
             '-c', $PROFILLIC_OPTIONS{ "subRateCat" },
             '-o', $PROFILLIC_OPTIONS{ "opt" },
             '-a', $PROFILLIC_OPTIONS{ "gamma" },
             '-f', $PROFILLIC_OPTIONS{ "freqRadio" },
             '-t', $PROFILLIC_OPTIONS{ "ttRatio" },
             '-s', $PROFILLIC_OPTIONS{ "treeImprovType" }
  );
  
  # Run Profillic
  if( $VERBOSE ) {
    print( "Running Profillic using command:\n" );
    print( join( ' ', @cmd ), "\n" );
  }
  # run it; redirect STDOUT and STDERR
  system( join( ' ', @cmd ) . "1>$profillicOutFile 2>$errFile" );
  if( $VERBOSE ) {
    print( "Done running Profillic." );
  }
  if( -s $errFile ) {
    open ERR, $errFile;
    my @lines = <ERR>;
    my $errMsg = join('', @lines);
    close ERR;
    die( "Error running Profillic: $errMsg" );
  }

  # Get pairwise distances.
  my $pwDistHashRef;
  my( @seqNames, @seqNamesWithDists );
  my $flag = my $count = 0;
  if( $VERBOSE ) {
    print( "Extracting pairwise distances from Profillic output in file $profillicOutFile\n" );
  }
  if( $VERBOSE ) {
    print( "Opening file $profillicOutFile for reading.." );
  }
  open( IN, $profillicOutFile ) or die( "Couldn't open $profillicOutFile: $!\n" );
  if( $VERBOSE ) {
    print( ".done\n" );
  }
  if( $VERBOSE ) {
    print( "Extracting pairwise distances from Profillic output.\n" );
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
      if( $line =~ /^(.*)$/ ) {
        my $seqnName = $1;				
        push @seqNamesWithDists, $seqnName;
        $seqnName =~ /^(\S+)\s+/;
        my $seqName = $1;
        push @seqNames, $seqName;
      }
    }
  }
  if( $VERBOSE ) {
    print( "Done.\nClosing file $profillicOutFile.." );
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
          next; # We just don't include these in the calculation.  Profillic pairwise diffs exclude gaps.
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
        warn "Warning: no distance values from $firstName\n"
      } else {
        for( $j = ( $i + 1 ); $j < scalar( @seqNames ); $j++ ) {
          $secondName = $seqNames[ $j ];
          if( !defined( $pwDistHashRef->{ $firstName }->{ $secondName } ) ) {
            warn "Warning: no distance value between $firstName and $secondName\n"
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

  print OUT "Number of sequences\tMean\tStandard error\tMin\tQ1\tMedian\tQ3\tMax\n";
  print OUT scalar( @seqNames ), "\t"; 
  print OUT join( "\t", @diversityResult ), "\n";

  if( $VERBOSE ) {
    print( "Done.\nClosing file $profillicOutFile.." );
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
} # runProfillic(..)

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
			$seqName = $seqCount . $1; # PAUL CHANGED; ADDED THE index
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
}
######################################################################################

runProfillic( @ARGV );

1;


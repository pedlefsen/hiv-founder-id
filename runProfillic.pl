#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) runProfillic
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
use vars qw( $opt_D $opt_V $opt_t );
use vars qw( $VERBOSE $DEBUG );

# For now I'm not actually running the trainer, just using the given alignment to make a profile (hmmer-style).
use constant PROFILLIC_EXECUTABLE => "profillic/dist/profillic_DNA_CQA";
use constant ALIGNEDFASTATOPROFILE_EXECUTABLE => "profuse/dist/alignedFastaToProfile_DNA";
use constant PROFILETOALIGNMENTPROFILE_EXECUTABLE => "profuse/dist/profileToAlignmentProfile_DNA";

## ERE I AM HAVE TO SET IT UP NOT WORRYING ABOUT PROFILLIC (training) JUST YET
Readonly my %PROFILLIC_OPTIONS => (
  random_seed => 98103,	                # use my zip code as a default.  Not much worse than any other.
  even_starting_profile_multiple => -1	# start at the starting profile.
        );
  
sub runProfillic {
  @ARGV = @_;

  sub runProfillic_usage {
    print "\trunProfillic [-DVt] <input_fasta_file> <output_alignment_profile_files_list_file> [<output_dir>]\n";
    exit;
  }

  # This means -D and -V are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # opt_t means train the profile (using profillic_DNA_CQA)
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_t ) = ();
  if( not getopts('DVt') ) {
    runProfillic_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;

  my $train_the_profile = $opt_t;
  
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

  my $output_alignment_profile_files_list_file = shift @ARGV || runProfillic_usage();
  
  my $output_path_dir = shift @ARGV ||
    $input_fasta_file_path . "/" . $input_fasta_file_short_nosuffix . "_hiv-founder-id_resultsDir";
  # Remove the trailing "/" if any
  if( defined( $output_path_dir ) ) {
    ( $output_path_dir ) = ( $output_path_dir =~ /^(.*[^\/])\/*$/ );
  }
  make_path( $output_path_dir );

  ## STEP 1: Make a consensus sequence.
  if( $VERBOSE ) {
    print "Calling R to create consensus sequence fasta file..";
  }
  my $R_output = `export computeConsensusSequenceFromAlignedFasta_inputFilename="$input_fasta_file"; export computeConsensusSequenceFromAlignedFasta_outputDir="$output_path_dir"; R -f computeConsensusSequenceFromAlignedFasta.R --vanilla --slave`;
  # The output has the file name of the consensus file.
  if( $DEBUG ) {
    print( "GOT: $R_output\n" );
  }
  if( $VERBOSE ) {
    print ".done.\n";
  }
  # Parse it to get the consensus sequence filename.
  my ( $consensus_fasta_file ) = ( $R_output =~ /\"([^\"]+)\"/ );
  if( $DEBUG ) {
    print "consensus file: $consensus_fasta_file\n";
  }

  ## STEP 2: Create a starting profile.
  my $starting_profile_file = "${output_path_dir}/${input_fasta_file_short_nosuffix}_AlignedFastaToProfile.prof";
  my $alignedfastatoprofileOutFile = "${output_path_dir}/${input_fasta_file_short_nosuffix}_AlignedFastaToProfile.out";
  my $alignedfastatoprofileErrFile = "${output_path_dir}/${input_fasta_file_short_nosuffix}_AlignedFastaToProfile.err";
  my @alignedfastatoprofile_cmd = (
             ALIGNEDFASTATOPROFILE_EXECUTABLE,
             $input_fasta_file,
             $consensus_fasta_file,
             $starting_profile_file
  );
  
  # Run AlignedFastaToProfile
  if( $VERBOSE ) {
    print( "Running AlignedFastaToProfile using command:\n" );
    print( join( ' ', @alignedfastatoprofile_cmd ), "\n" );
  }
  # run it; redirect STDOUT and STDERR
  system( join( ' ', @alignedfastatoprofile_cmd ) . " 1>$alignedfastatoprofileOutFile 2>$alignedfastatoprofileErrFile" );
  if( $VERBOSE ) {
    print( "Done running AlignedFastaToProfile." );
  }
  if( -s $alignedfastatoprofileErrFile ) {
    open ERR, $alignedfastatoprofileErrFile;
    my @lines = <ERR>;
    my $errMsg = join('', @lines);
    close ERR;
    die( "Error running AlignedFastaToProfile: $errMsg" );
  }

  ## STEP 4: Create an ungapped version
  ## TODO: DEHACKIFY! Here we are lazily not protecting the sequence names from also losing dashes.
  my $input_fasta_file_contents = path( $input_fasta_file )->slurp();
  $input_fasta_file_contents =~ s/-//g;
  
  # Now write it out to a temporary location in the output dir.
  my $ungapped_fasta_file = "${output_path_dir}/${input_fasta_file_short_nosuffix}_ungapped${input_fasta_file_suffix}";
  if( $VERBOSE ) {
    print( "Writing out ungapped input file \"$ungapped_fasta_file\".." );
  }
  if( $VERBOSE ) { print "Opening file \"$ungapped_fasta_file for writing..\n"; }
  unless( open ungapped_fasta_fileFH, ">$ungapped_fasta_file" ) {
    warn "Unable to open output file \"$ungapped_fasta_file: $!\n";
    return 1;
    }
  print ungapped_fasta_fileFH $input_fasta_file_contents;
  close( ungapped_fasta_fileFH );
    
  ## STEP 5: Ensure the profile is at a/the maximum-likelihood point.
  my $trained_profile_file = $starting_profile_file;
  if( $train_the_profile ) {
    $trained_profile_file = "${output_path_dir}/${input_fasta_file_short_nosuffix}_Profillic.prof";
    my $profillicOutFile = "${output_path_dir}/${input_fasta_file_short_nosuffix}_Profillic.out";
    my $profillicErrFile = "${output_path_dir}/${input_fasta_file_short_nosuffix}_Profillic.err";
    my @profillic_cmd = (
               PROFILLIC_EXECUTABLE,
               $trained_profile_file,
               $ungapped_fasta_file,
               '--starting_profile', $starting_profile_file,
               '--even_starting_profile_multiple', $PROFILLIC_OPTIONS{ "even_starting_profile_multiple" },
               '--random_seed', $PROFILLIC_OPTIONS{ "random_seed" }
      );
    
    # Run Profillic
    if( $VERBOSE ) {
      print( "Running Profillic using command:\n" );
      print( join( ' ', @profillic_cmd ), "\n" );
    }
    # run it; redirect STDOUT and STDERR
    system( join( ' ', @profillic_cmd ) . " 1>$profillicOutFile 2>$profillicErrFile" );
    if( $VERBOSE ) {
      print( "Done running Profillic." );
    }
    if( -s $profillicErrFile ) {
      open ERR, $profillicErrFile;
      my @lines = <ERR>;
      my $errMsg = join('', @lines);
      close ERR;
      die( "Error running Profillic: $errMsg" );
      }
  } # End if $train_the_profile

  ## STEP 6: Create Alignment Profiles
  my $alignment_profiles_output_file_prefix = "${output_path_dir}/${input_fasta_file_short_nosuffix}_profileToAlignmentProfile";
  my $profiletoalignmentprofileOutFile = $output_alignment_profile_files_list_file; #"${output_path_dir}/${input_fasta_file_short_nosuffix}_profileToAlignmentProfile.out";
  my $profiletoalignmentprofileErrFile = "${output_path_dir}/${input_fasta_file_short_nosuffix}_profileToAlignmentProfile.err";
  my @profiletoalignmentprofile_cmd = (
    PROFILETOALIGNMENTPROFILE_EXECUTABLE,
    $trained_profile_file,
    $ungapped_fasta_file,
    $alignment_profiles_output_file_prefix,
    "--individual" # means create individual alignment profiles, one per entry in the given fasta file.
  );
  
  # Run profileToAlignmentProfile
  if( $VERBOSE ) {
    print( "Running profileToAlignmentProfile using command:\n" );
    print( join( ' ', @profiletoalignmentprofile_cmd ), "\n" );
  }
  # run it; redirect STDOUT and STDERR
  system( join( ' ', @profiletoalignmentprofile_cmd ) . " 1>$profiletoalignmentprofileOutFile 2>$profiletoalignmentprofileErrFile" );
  if( $VERBOSE ) {
    print( "Done running profileToAlignmentProfile." );
  }
  # This STDERR output is just verbose feedback; write it out to STDOUT.
  if( $VERBOSE && ( -s $profiletoalignmentprofileErrFile ) ) {
    open ERR, $profiletoalignmentprofileErrFile;
    my @lines = <ERR>;
    my $errMsg = join('', @lines);
    close ERR;
    print( "profileToAlignmentProfile output:\n$errMsg" );
  }

  if( $VERBOSE ) {
    select STDOUT;
    $| = $old_autoflush;
  }
  
  if( $VERBOSE ) {
    print( "Alignment profiles are listed in $output_alignment_profile_files_list_file\n" );
  }
  return 0;
} # runProfillic(..)

runProfillic( @ARGV );

1;


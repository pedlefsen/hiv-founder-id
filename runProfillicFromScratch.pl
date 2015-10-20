#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) runProfillicFromScratch
##  Author:
##      Paul Thatcher Edlefsen   pedlefse@fredhutch.org
##  Description:
##
##      Script for running profillic locally to get the profile estimated
##      and alignment profiles generated, starting from unaligned sequences.
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
use vars qw( $opt_D $opt_V $opt_e );
use vars qw( $VERBOSE $DEBUG );

# For now I'm not actually running the trainer, just using the given alignment to make a profile (hmmer-style).
use constant PROFILLIC_EXECUTABLE => "profillic/dist/profillic_DNA_CQA";
use constant ALIGNEDFASTATOPROFILE_EXECUTABLE => "profuse/dist/alignedFastaToProfile_DNA";
use constant PROFILETOALIGNMENTPROFILE_EXECUTABLE => "profuse/dist/profileToAlignmentProfile_DNA";

Readonly my %PROFILLIC_OPTIONS => (
  random_seed => 98103,	                # use my zip code as a default.  Not much worse than any other.
  even_starting_profile_multiple => -1	# start at the starting profile.
        );
  
sub runProfillicFromScratch {
  @ARGV = @_;

  sub runProfillicFromScratch_usage {
    print "\trunProfillicFromScratch [-DVe] <input_fasta_file> <output_alignment_profile_files_list_file> [<output_dir>]\n";
    exit;
  }

  # This means -D and -V are ok, but nothin' else.
  # opt_D means print debugging output.
  # opt_V means be verbose.
  # opt_e means compute the self-entropy of the profile
  # But first reset the opt vars.
  ( $opt_D, $opt_V, $opt_e ) = ();
  if( not getopts('DVe') ) {
    runProfillicFromScratch_usage();
  }
  
  $DEBUG ||= $opt_D;
  $VERBOSE ||= $opt_V;

  my $compute_self_entropy = $opt_e;
  
  my $old_autoflush;
  if( $VERBOSE ) {
    select STDOUT;
    $old_autoflush = $|;
    $| = 1; # Autoflush.
  }
  
  my $input_fasta_file = shift @ARGV || runProfillicFromScratch_usage();
  my ( $input_fasta_file_path, $input_fasta_file_short ) =
    ( $input_fasta_file =~ /^(.*?)\/([^\/]+)$/ );
  unless( $input_fasta_file_short ) {
    $input_fasta_file_short = $input_fasta_file;
    $input_fasta_file_path = ".";
  }
  my ( $input_fasta_file_short_nosuffix, $input_fasta_file_suffix ) =
    ( $input_fasta_file_short =~ /^([^\.]+)(\..+)?$/ );

  my $output_alignment_profile_files_list_file = shift @ARGV || runProfillicFromScratch_usage();
  
  my $output_path_dir = shift @ARGV ||
    $input_fasta_file_path . "/" . $input_fasta_file_short_nosuffix . "_hiv-founder-id_resultsDir";
  # Remove the trailing "/" if any
  if( defined( $output_path_dir ) ) {
    ( $output_path_dir ) = ( $output_path_dir =~ /^(.*[^\/])\/*$/ );
  }
  make_path( $output_path_dir );

  ## STEP 1: Create an ungapped version of the input file
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
    
  ## STEP 2: Use profillic to train the profile to a/the maximum-likelihood point.
  my $trained_profile_file = "${output_path_dir}/${input_fasta_file_short_nosuffix}_Profillic.prof";
  my $profillicOutFile = "${output_path_dir}/${input_fasta_file_short_nosuffix}_Profillic.out";
  my $profillicErrFile = "${output_path_dir}/${input_fasta_file_short_nosuffix}_Profillic.err";
  my @profillic_cmd = (
               PROFILLIC_EXECUTABLE,
               $trained_profile_file,
               $ungapped_fasta_file,
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
  
  # Step 3: (maybe) calculate profile self-entropy
  if( $compute_self_entropy ) {
    my $profilecrossentropyOutFile = "${output_path_dir}/${input_fasta_file_short_nosuffix}_ProfileCrossEntropy.out";
    my $profilecrossentropyErrFile = "${output_path_dir}/${input_fasta_file_short_nosuffix}_ProfileCrossEntropy.err";
    if( $VERBOSE ) {
      print( "Running ProfileCrossEntropy using command:\n" );
      print( join( ' ', @profilecrossentropy_cmd ), "\n" );
    }
    # run it; redirect STDOUT and STDERR
    system( join( ' ', @profilecrossentropy_cmd ) . " 1>$profilecrossentropyOutFile 2>$profilecrossentropyErrFile" );
    if( $VERBOSE ) {
      print( "Done running ProfileCrossEntropy." );
    }
    if( -s $profilecrossentropyErrFile ) {
      open ERR, $profilecrossentropyErrFile;
      my @lines = <ERR>;
      my $errMsg = join('', @lines);
      close ERR;
      die( "Error running ProfileCrossEntropy: $errMsg" );
    }
    if( -s $profilecrossentropyOutFile ) {
      open PROFILE_CROSS_ENTROPY, $profilecrossentropyErrFile;
      my @lines = <PROFILE_CROSS_ENTROPY>;
      my $crossEntropy = join( '', @lines );
      close PROFILE_CROSS_ENTROPY;

      print $crossEntropy, "\n";
    }
  } # End if $compute_self_entropy

  ## STEP 4: Create Alignment Profiles
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
} # runProfillicFromScratch(..)

runProfillicFromScratch( @ARGV );

1;


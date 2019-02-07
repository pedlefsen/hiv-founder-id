#!/usr/bin/perl -w

#{{{ Load Libraries
# set foldmethod=marker in vim. markers are the defaults: {{{ and }}}
use Getopt::Std; # for getopts
use File::Path qw( make_path );
use File::Temp qw / :POSIX /;
use Path::Tiny;
use Try::Tiny;
use Cwd;
use File::Spec;

# For screenscraping
use WWW::Mechanize;
use LWP::Simple;

# Probably just because PL is noob.
use IO::Socket::SSL qw();
use Data::Dumper;

use strict;
use vars qw( $opt_D $opt_V $opt_e $opt_p $opt_P $opt_R $opt_x );
use vars qw( $VERBOSE $DEBUG );
# }}}

my $emailAddress = "pedlefsen\@gmail.com";

# {{{ Documentation for %opt:
# i => input file
# o => output file name (for this script)
# e => ElimDupes deduplicated file name (name of output file produced by elimdupes). Populated by parseARGV.
# d => ElimDupes duplicates table list produced by ElimDupes. Populated by parseARGV.
# f => Filename for the final output if recombination was found. Populated by parseARGV.
# u => Number of duplicates removed. Populated by ElimDupes.
# V => Verbose or not?
# j => job_id. Populated by runRAPROnline.
# r => filename with the stats of the recombinants. Populated by parseARGV and used be retrieveResults.
# n => number of recombinants removed
# a => Active file -> the one that will be operated on next. }}}

sub parseARGV { #{{{
  @ARGV = @_;
  ( $opt_V ) = ();
  my %opt = ();
  if( not getopts("Vi:o:", \%opt) ) {
    print "\nERROR";
    print "\nrunRAPROnline.pl [-V] -i <input_fasta_file> -o <output_fasta_file>";
    print "\nFile names must end in lowercase .fasta\n";
    exit;
  }

  my $cwd = getcwd();
  my $rel_input = File::Spec->abs2rel($opt{'i'}, $cwd);
  my $rel_output = File::Spec->abs2rel($opt{'o'}, $cwd);

  if ( $opt{'V'} ){
    print "\n\nFor some reason absolute paths leads to errors. Remapping:\n$opt{'i'}\n$rel_input\n\n$opt{'o'}\n$rel_output\n\n";
  }
  
  $opt{'i'} = $rel_input;
  $opt{'o'} = $rel_output;

  unless (-e $opt{'i'}){
    print "\nERROR: Specified input file does not exist.\n";
    exit;
  }

  if (index($opt{'o'}, ".fasta") == -1){
    print "\nERROR: Specified output file does not end with .fasta (lowercase).\n";
    exit;
  }

  if (index($opt{'i'}, ".fasta") == -1){
    print "\nERROR: Specified input file does not end with .fasta (lowercase).\n";
    exit;
  }

  $opt{'e'} = $opt{'o'} =~ s/.fasta/_elimdupes.fasta/gr;
  $opt{'d'} = $opt{'o'} =~ s/.fasta/_elimdupes.csv/gr;
  $opt{'r'} = $opt{'o'} =~ s/.fasta/_rec_list.tab/gr;
  $opt{'f'} = $opt{'o'} =~ s/.fasta/_RAPR.fasta/gr;
  $opt{'a'} = $opt{'i'};

  $opt{'n'} = -1;
  return %opt;
} #}}} parseARGV

sub runElimDupes { #{{{
  my %opt = @_;

  if ($opt{'V'}){
    foreach my $key (sort(keys %opt)){
      print "In \%opt, $key has value $opt{$key}\n";
    }
  }

  my $ElimDupes_command = "Rscript ElimDupes.R --verbose --input_file=$opt{'a'} --output_file=$opt{'e'}";

  my $ElimDupes_result = `$ElimDupes_command`;

  if ($opt{'V'}){
    print "\n$ElimDupes_command\n";
    print "\n$ElimDupes_result\n";
  }

  my @number_of_dups_removed_string = $ElimDupes_result =~ /^\[1\] "Number of duplicates removed: [0-9]*"/mg;
  my $size = @number_of_dups_removed_string;
  my $number_of_dups_removed = $number_of_dups_removed_string[0];
  if ($size != 1){
    print "\nOutput from ElimDupes.R not parsed correctly\n";
    exit;
  } else {
    $number_of_dups_removed =~ s/^\[1\] "Number of duplicates removed: ([0-9]*)"/$1/g;
    # always put it in the ElimDupes format even if there are no dupes, so skip the test.
    #if ($number_of_dups_removed > 0){
      $opt{'a'} = $opt{'e'};
    #}
  }
  return %opt;
} #}}}

sub runRAPROnline { #{{{

  my %opt = @_;

  my $job_id = "NOTRUN";

  my $output_fasta_file = $opt{'o'};
        my ( $output_fasta_file_path, $output_fasta_file_short ) =
          ( $output_fasta_file =~ /^(.*?)\/([^\/]+)$/ );

  my $input_fasta_file = $opt{'a'};
        my ( $input_fasta_file_path, $input_fasta_file_short ) =
          ( $input_fasta_file =~ /^(.*?)\/([^\/]+)$/ );
        unless( $input_fasta_file_short ) {
            $input_fasta_file_short = $input_fasta_file;
            $input_fasta_file_path = ".";
          }
        my ( $input_fasta_file_short_nosuffix, $input_fasta_file_suffix ) =
          ( $input_fasta_file_short =~ /^([^\.]+)(\..+)?$/ );

  $VERBOSE = $opt{'V'};

  my $R_output;

  # RAP has a problem with certain characters in fasta headers. 
  ## Step 1: Set up the table mapping sequence names to numbers.
      if( $VERBOSE ) {
        print "Reading sequence names from file \"", $input_fasta_file, "\"..";
      }
      my $input_fasta_file_contents =
           path( $input_fasta_file )->slurp();
      if( $VERBOSE ) {
        print ".done\n";
      }
      if( $DEBUG ) {
        #print $input_fasta_file_contents;
      }
      my ( @seq_names ) =  ( $input_fasta_file_contents =~ /\n?>[ \t]*(.+) *\n/g );

  ## STEP 2: Rename the seqs to just their numbers.
  my ( $input_fasta_file_numbered ) =
    ( $input_fasta_file =~ /^(.+)$input_fasta_file_suffix$/ );
  $input_fasta_file_numbered .= "_numbered$input_fasta_file_suffix";
  if( $VERBOSE ) {
    print "Calling R to create a version of the fasta file in which seqs are numbered instead of named ($input_fasta_file_numbered)..";
  }

  $R_output = `export computeConsensusSequenceFromAlignedFasta_inputFilename="$input_fasta_file"; export computeConsensusSequenceFromAlignedFasta_outputFilename="$input_fasta_file_numbered"; export computeConsensusSequenceFromAlignedFasta_includeFullAlignment="TRUE"; export  computeConsensusSequenceFromAlignedFasta_includeConsensus="FALSE"; export computeConsensusSequenceFromAlignedFasta_useSeqeunceNumbersAsNames="TRUE"; export computeConsensusSequenceFromAlignedFasta_addSuffixToSeqNames="_1"; R -f computeConsensusSequenceFromAlignedFasta.R --vanilla --slave`;
  # The output has the file name of the "consensus file" which is not in fact the consensus but the fasta with renamed seqs.
  if( $DEBUG ) {
    print( "GOT: $R_output\n" );
  }
  if( $VERBOSE ) {
    print ".done.\n";
  }
  ## Parse it to get the filename.
  #my ( $fasta_file_readyForRAP ) = ( $R_output =~ /\"([^\"]+)\"/ );
  my ( $fasta_file_readyForRAP ) = $input_fasta_file_numbered;
  if( $VERBOSE ) {
    print "Fasta file almost ready for RAP: $fasta_file_readyForRAP\n";
  }
  
  # Also fix the filename.
  my $fasta_file_readyForLANL = File::Temp::tempnam( ".", "tmp" ) . ".fasta";
  # remove leading dots and slashes from filename:
  ( $fasta_file_readyForLANL ) = ( $fasta_file_readyForLANL =~ /^\.\/(.+)$/ );
  `cp $fasta_file_readyForRAP $fasta_file_readyForLANL`;

  if( $VERBOSE ) {
    print "Fasta file finally ready for RAPR on LANL: $fasta_file_readyForLANL\n";
  }
  my $mech = WWW::Mechanize->new( ssl_opts => {
      SSL_verify_mode => IO::Socket::SSL::SSL_VERIFY_NONE,
      verify_hostname => 0, # this key is likely going to be removed in future LWP >6.04
    } );
  
  $mech->get( "https://www.hiv.lanl.gov/content/sequence/RAP2017/rap.html" );
  
  # {{{ Verbose output
  if( $opt{'V'} ) {
      print "=" x50;
      print "\n\nforms\n\n";
      print "=" x50;
      print "\n\n";
      
      foreach my $form ($mech->forms) {
          print "$form";
          print "\n";
      }
      print "\n";
  } #}}}
  
  my @forms = $mech->forms;
  
  # {{{ Verbose output
  if( $opt{'V'} ) {
      foreach my $form (@forms) {
        my @inputfields = $form->param;
        print "=" x50;
        print "\n\n$form\n\n";
        print "=" x50;
        print "\n\nForm Input List:\n\n";
        print Dumper \@inputfields;
        print "\n\nForm Input Details:\n\n";
        print Dumper $form->inputs;
      }
  } #}}}
  
  my %fields = (
    seq_input => '',
    alignmentFile => $fasta_file_readyForLANL,
    fdr => '0.1',
    ntf => '0',#'0', # 1 tf
    tp => '0', # no time points available
    elimdupes => '2',#'2',
    alwaysemail => '1'
  );
  # Additional Fields that do not have to be specified: {{{
  #  GroupDetermining => 'col',
  #  ColumnNum => '1',
  #  ColumnDelimitor => '.',
  #  NumOfChars => '1',
  #  csv_input => '',
  #  GroupsFile => '', }}}
  if( $opt{'V'} ){
    print "\nFields to be submitted:\n";
    print map { "$_ => $fields{$_}\n" } keys %fields;
  }
  
  my $response = $mech->submit_form(  form_number => 2,
                                  fields    => \%fields
                       );
  
  # {{{ Verbose output
  if( $opt{'V'} ) {
    print "=" x50;
    print "\n\nResult of submitting:\n\n";
    print "=" x50;
    print "\n\n";
  } #}}}
  
  if ($response->is_success) {
    my $result = $response->decoded_content;
  
    if ( $result =~ /If you want to proceed, please enter the email address below/ ) {
      if ($opt{'V'}){
        print "\nSubmitting Email Adress.\n";
      }
  
      # extract job id
      my @forms = $mech->forms;
      my $c_form = $forms[1];
      my @names = $c_form->param;
      $job_id = $c_form->param( "JOB_ID" );
      $job_id =~ s/\/tmp\/RAP\///ig;
  
      if ( $opt{'V'} ){
            print Dumper @names;
            print "$job_id\n";
      }
  
      $mech->submit_form( form_number => 2,
                          fields      => {
                                          titleFromUser => "$input_fasta_file",
                                          EMAIL => $emailAddress,
                                          EMAIL2 => $emailAddress
                                         }
      );
      if ($response->is_success) {
        my $result = $response->decoded_content;
        if ( $opt{'V'} ){
          print $result;
        }
      } else {
        print "\nERROR: submission of email address failed: $response->status_line \n";
        exit;
      }
  
    } else {
      print "\nError: Expecting a request for an email address. Instead, the result of no recombinants were already obtained.\n";
      exit;
    }
  } else {
    print "\nError: Submission of input file failed: $response->status_line \n";
    exit;
  }
  $opt{'j'} = $job_id;
  
  # Store these for renaming seqs.
  $opt{'seq_names'} = \@seq_names;
  return %opt;
} #}}}

sub retrieveResult { #{{{
  
  my %opt = @_;

  my $job_id = $opt{'j'};
  my $url_to_grab = "https://www.hiv.lanl.gov/tmp/RAP/$job_id/out.html";

  if ( $opt{'V'} ){
    print "\n$job_id\n";
    print "\n$url_to_grab\n";
  }

  my $mech = WWW::Mechanize->new( ssl_opts => {
      SSL_verify_mode => IO::Socket::SSL::SSL_VERIFY_NONE,
      verify_hostname => 0, # this key is likely going to be removed in future LWP >6.04
    } );

  my $times_checked = 0;

  while ( $times_checked < 10){
    $mech->get( $url_to_grab );
    my $response = $mech->content;

    if ( $opt{'V'} ){ #{{{ Verbose output
      print "=" x50;
      print "\n\nRESPONSE after retrieving the result.\n\n";
      print "=" x50;
      print "\n\n";
      print $response;
      print "\n\n";
      print "=" x50;
      print "\n\n";
    } # }}}
    
    if ( $response =~ /List of Recombinants/ ){
      my $url_with_recs = "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/tmp/RAP/$job_id/sample.recs.txt";
      my $RAP_output_file = $opt{'r'};
      $mech->get( $url_with_recs, ":content_file" => $RAP_output_file );
      
      # Fix the numbers back to sequence names.
      ## TODO: DEHACKIFY? KEEP A "BEFORE" VERSION.
      `cp $RAP_output_file ${RAP_output_file}_beforeFixingNumbersToSequenceNames.tab`;
      
      my $RAP_output_file_contents = path( $RAP_output_file )->slurp();
      my @RAP_output_file_lines = split( "\n", $RAP_output_file_contents );
      my @seq_names = @{ $opt{'seq_names'} };
      my ( @recombinant_line, $rl0, $rl4, $rl6 );
      for( my $line_i = 1; $line_i < scalar( @RAP_output_file_lines ); $line_i++ ) {
        ## TODO: REMOVE!
        print "LINE $line_i: $RAP_output_file_lines[ $line_i ]\n";
          @recombinant_line = split( /\t/, $RAP_output_file_lines[ $line_i ] );
          # Relevant columns are 0, 4, 6.
          ( $rl0 ) = ( $recombinant_line[ 0 ] =~ /^(\d+)_/ ); # strip off the suffix if there is one.
          ( $rl4 ) = ( $recombinant_line[ 4 ] =~ /^(\d+)_/ ); # strip off the suffix if there is one.
          ( $rl6 ) = ( $recombinant_line[ 6 ] =~ /^(\d+)_/ ); # strip off the suffix if there is one.
          $recombinant_line[ 0 ] = $seq_names[ $rl0 - 1 ];
          $recombinant_line[ 4 ] = $seq_names[ $rl4 - 1 ];
          $recombinant_line[ 6 ] = $seq_names[ $rl6 - 1 ];
          $RAP_output_file_lines[ $line_i ] = join( "\t", @recombinant_line );
        }
      $RAP_output_file_contents = join( "\n", @RAP_output_file_lines ) . "\n";

        if( $VERBOSE ) {
          print( "Writing out RAP output file \"$RAP_output_file\".." );
        }
        if( $VERBOSE ) { print "Opening file \"$RAP_output_file\" for writing..\n"; }
        unless( open RAP_output_fileFH, ">$RAP_output_file" ) {
            warn "Unable to open output file \"$RAP_output_file\": $!\n";
            return 1;
          }
        print RAP_output_fileFH $RAP_output_file_contents;
      close( RAP_output_fileFH );
  
      print "Recombinants identified ($RAP_output_file)\n";

      # NOW call ElimDupes.R and make it remove the sequences in the recombined list.
      my $second_ElimDupes_command = "Rscript ElimDupes.R --verbose --input_file=$opt{'a'} --output_file=$opt{'f'} --reduplicate --duplicates_file=$opt{'d'} --but_first_remove=$RAP_output_file";
      if ($opt{'V'}){
        print "\n\n\n\n";
        print $second_ElimDupes_command;
        print "\n\n\n\n";
      }
      my $second_ElimDupes_result = `$second_ElimDupes_command`;
      
      if ( $opt{'V'} ){
        print "\n\n\nSecond ElimDupes Result\n\n\n";
        print $second_ElimDupes_result;
        print "\n\n";
      }

      $opt{'a'} = $opt{'f'};
      $opt{'n'} = `wc -l < $opt{'r'}` - 1;
      print "\nTotal number of recombinants found: $opt{'n'}\n";
      print "\nReduplicated file: $opt{'f'}\n";

      $times_checked = 100;
    } elsif ( $response =~ /No recombinants were found in the submitted alignment/ ) {
      print "\nNo recombinants were found\n";
      $opt{'a'} = $opt{'i'};
      $opt{'n'} = 0;
      $times_checked = 100;
    } elsif ( $response =~ /We encountered a problem with your input. Please consider these tips/ || $response =~ /Sorry!/  ) {
      print "\nError occurred. LANL rejected form\n";
      $opt{'n'} = -1;
      $times_checked = 100;
    } else {
      $times_checked = $times_checked + 1;
      if ($opt{'V'}){
        print "\nStill waiting for the result - sleeping for the $times_checked\'th time.\n";
        print "Checking the url $url_to_grab\n";
      }
      sleep 5 + 10*$times_checked;
    }
  }
  return(%opt);
} #}}}

# main loop {{{
# DEBUG {{{
#my %opt = ();
#print "\n\n\nWARNING\n\n\nOVERRIDING opt with DEBUG settings\n\n\n";
#%opt = (
#  e => '/tmp/rapr_samp_out_elimdupes.fasta',
#  f => '/tmp/rapr_samp_out_RAPR.fasta',
#  d => '/tmp/rapr_samp_out_elimdupes.csv',
#  j => 'aeATy0Xc3v',
#  V => 0,
#  a => '/tmp/rapr_samp_out_elimdupes.fasta',
#  u => 5,
#  o => '/tmp/rapr_samp_out.fasta',
#  r => '/tmp/rapr_samp_out_rec_list.tab',
#  i => 'tests/data/rapr_recomb.fasta'
#);
#print Dumper %opt;
# }}}

my %opt = parseARGV( @ARGV );
if ( $opt{'V'} ){
  print "\nValue of \%opt after parseARGV:\n";
  print map { "$_ => $opt{$_}\n" } keys %opt;
}

%opt = runElimDupes( %opt );

if ( $opt{'V'} ){
  print "\nValue of \%opt after runElimDupes:\n";
  print map { "$_ => $opt{$_}\n" } keys %opt;
}

%opt = runRAPROnline( %opt );
if ( $opt{'V'} ){
  print "\nValue of \%opt after runRAPROnline:\n";
  print map { "$_ => $opt{$_}\n" } keys %opt;
}

%opt = retrieveResult( %opt );
if ( $opt{'V'} ){
  print "\nValue of \%opt after retrieveResult:\n";
  print map { "$_ => $opt{$_}\n" } keys %opt;
}

1; #}}}


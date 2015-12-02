source( "removeDuplicateSequencesFromAlignedFasta_safetosource.R" )

## Here is where the action is.
input.fasta.file <- Sys.getenv( "removeDuplicateSequencesFromAlignedFasta_inputFilename" );
output.fasta.file <- Sys.getenv( "removeDuplicateSequencesFromAlignedFasta_outputFastaFilename" );
if( nchar( output.fasta.file ) == 0 ) {
    output.fasta.file <- NULL;
}
output.dir <- Sys.getenv( "removeDuplicateSequencesFromAlignedFasta_outputDir" );
if( nchar( output.dir ) == 0 ) {
    output.dir <- NULL;
}
output.table.file <- Sys.getenv( "removeDuplicateSequencesFromAlignedFasta_outputTableFilename" );
if( nchar( output.table.file ) == 0 ) {
    output.table.file <- NULL;
}
include.gaps.in.Hamming <- Sys.getenv( "runPoissonFitter_includeGapsInHamming" );
if( ( include.gaps.in.Hamming == "" ) || ( toupper( include.gaps.in.Hamming ) == "F" ) || ( toupper( include.gaps.in.Hamming ) == "FALSE" ) || ( include.gaps.in.Hamming == "0" ) ) {
    include.gaps.in.Hamming = FALSE;
} else {
    include.gaps.in.Hamming = TRUE;
}
increase.threshold.to.ensure.max <- Sys.getenv( "removeDuplicateSequencesFromAlignedFasta_increaseThresholdToEnsureMax" );
if( nchar( increase.threshold.to.ensure.max ) == 0 ) {
    increase.threshold.to.ensure.max <- NULL;
}

## TODO: REMOVE
# warning( paste( "aligned fasta input file:", input.fasta.file ) );
# if( !is.null( output.dir ) ) {
#     warning( paste( "consensus fasta output dir:", output.dir ) );
# }
# if( !is.null( output.fasta.file ) ) {
#     warning( paste( "consensus fasta output file:", output.fasta.file ) );
# }
if( file.exists( input.fasta.file ) ) {
    print( removeDuplicateSequencesFromAlignedFasta( input.fasta.file, output.dir = output.dir, output.fasta.file = output.fasta.file, output.table.file = output.table.file, include.gaps.in.Hamming = include.gaps.in.Hamming, increase.threshold.to.ensure.max = increase.threshold.to.ensure.max ) );
} else {
    stop( paste( "File does not exist:", input.fasta.file ) );
}

source( "maskSynonymousCodonsInAlignedFasta_safetosource.R" )

## Here is where the action is.
input.codon.fasta.file <- Sys.getenv( "maskSynonymousCodonsInAlignedFasta_inputFilename" );
input.aa.fasta.file <- Sys.getenv( "maskSynonymousCodonsInAlignedFasta_translatedInputFilename" );
if( nchar( input.aa.fasta.file ) == 0 ) {
    input.aa.fasta.file <- NULL;
}
output.fasta.file <- Sys.getenv( "maskSynonymousCodonsInAlignedFasta_outputFilename" );
if( nchar( output.fasta.file ) == 0 ) {
    output.fasta.file <- NULL;
}
output.dir <- Sys.getenv( "maskSynonymousCodonsInAlignedFasta_outputDir" );
if( nchar( output.dir ) == 0 ) {
    output.dir <- NULL;
}

## TODO: REMOVE
# warning( paste( "aligned fasta input file:", input.codon.fasta.file ) );
# if( !is.null( output.dir ) ) {
#     warning( paste( "consensus fasta output dir:", output.dir ) );
# }
# if( !is.null( output.fasta.file ) ) {
#     warning( paste( "consensus fasta output file:", output.fasta.file ) );
# }
if( file.exists( input.codon.fasta.file ) ) {
    print( maskSynonymousCodonsInAlignedFasta( input.codon.fasta.file, input.aa.fasta.file, output.dir = output.dir, output.fasta.file = output.fasta.file ) );
} else {
    stop( paste( "File does not exist:", input.codon.fasta.file ) );
}

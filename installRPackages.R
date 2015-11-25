# Look for, and also install, packages in "~/R/library"
.libPaths( c( "~/R/Library", .libPaths() ) );
source( "https://bioconductor.org/biocLite.R" );

install.packages( pkgs = "binom", dependencies = TRUE, quiet = TRUE );
if( !require( "binom" ) ) {
    stop( "Error loading package \"binom\"" );
}
install.packages( pkgs = "xtable", dependencies = TRUE, quiet = TRUE );
if( !require( "xtable" ) ) {
    stop( "Error loading package \"xtable\"" );
}
install.packages( pkgs = "coin", dependencies = TRUE, quiet = TRUE );
if( !require( "coin" ) ) {
    stop( "Error loading package \"coin\"" );
}
install.packages( pkgs = "ggplot2", dependencies = TRUE, quiet = TRUE );
if( !require( "ggplot2" ) ) {
    stop( "Error loading package \"ggplot2\"" );
}
install.packages( pkgs = "tools", dependencies = TRUE, quiet = TRUE );
if( !require( "tools" ) ) {
    stop( "Error loading package \"tools\"" );
}
install.packages( pkgs = "ade4", dependencies = TRUE, quiet = TRUE );
if( !require( "ade4" ) ) {
    stop( "Error loading package \"ade4\"" );
}
install.packages( pkgs = "ape", dependencies = TRUE, quiet = TRUE );
if( !require( "ape" ) ) {
    stop( "Error loading package \"ape\"" );
}
install.packages( pkgs = "dynamicTreeCut", dependencies = TRUE, quiet = TRUE );
if( !require( "dynamicTreeCut" ) ) {
    stop( "Error loading package \"dynamicTreeCut\"" );
}
install.packages( pkgs = "entropy", dependencies = TRUE, quiet = TRUE );
if( !require( "entropy" ) ) {
    stop( "Error loading package \"entropy\"" );
}
biocLite( pkgs = eval( "seqinr" ) );
if( !require( package = eval( "seqinr" ) ) ) {
    stop( paste( "Error loading package", "seqinr" ) );
}

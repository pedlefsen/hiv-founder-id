# Look for, and also install, packages in "~/R/library"
.libPaths( c( "~/R/Library", .libPaths() ) );
source( "http://bioconductor.org/biocLite.R" );

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
install.packages( pkgs = "ROCR", dependencies = TRUE, quiet = TRUE );
if( !require( "ROCR" ) ) {
    stop( "Error loading package \"ROCR\"" );
}
install.packages( pkgs = "glmnet", dependencies = TRUE, quiet = TRUE );
if( !require( "glmnet" ) ) {
    stop( "Error loading package \"glmnet\"" );
}
install.packages( pkgs = "glmnetUtils", dependencies = TRUE, quiet = TRUE );
if( !require( "glmnetUtils" ) ) {
    stop( "Error loading package \"glmnetUtils\"" );
}

biocLite( pkgs = "seqinr" );
if( !require( package = "seqinr" ) ) {
    stop( paste( "Error loading package", "seqinr" ) );
}
biocLite( pkgs = "Biostrings" );
if( !require( package = "Biostrings" ) ) {
    stop( paste( "Error loading package", "Biostrings" ) );
}

install.packages( pkgs = "optparse", dependencies = TRUE, quiet = TRUE );
if( !require( "optparse" ) ) {
    stop( "Error loading package \"optparse\"" );
}

## TO INSTALL hypermuteR on linux:
# (in R:)
library( "devtools" );
install_github( "philliplab/hypermutR" )
## BUT on my mac I had to rebuild it from source:
# (in shell:)
# git clone https://github.com/philliplab/hypermutR.git
# rm hypermutR/src/*.o hypermutR/src/*.so
# sudo R CMD INSTALL hypermutR
# Test it loads (in R:)
if( !require( package = "hypermutR" ) ) {
    stop( paste( "Error loading package", "hypermutR" ) );
}

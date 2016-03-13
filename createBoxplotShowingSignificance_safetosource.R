
library( "ggplot2" )

## Only works for 2 to 4 groups.  Only tested on 2 and 4 groups.
## test.results.as.list.of.lists should contain top-level entries for levels 1..(n-1) for n groups; each entry [[i]] should contain a sublist with entries for levels (i+1)..n.  Each entry should be a list with a "p.value" element defined.  Each entry must be named using names from the levels of the "group" factor var.
## NOTE BY DEFAULT, assumes the.var has continuous scale.  For discrete scale y, set y.scale.is.discrete = TRUE
createBoxplotShowingSignificance <- function ( group, the.var, the.var.name, test.results.as.list.of.lists, the.var.scale.is.discrete = FALSE, p.value.threshold = 0.05, the.title = NULL )
{
## TODO: LOTS TO DO TO MAKE THIS GENERALLY WORK.  FOR NOW WRITING IT FOR FOUR GROUPS MAX.
# Perhaps it depends on the number of significant values, ultimately; but for now (since we assume four groups) we can just lazily always make room for the significance bars.

  .the.data.frame <- data.frame( outcome.var = the.var, group.var = group );

  .n.groups <- length( levels( group ) );
  stopifnot( .n.groups >= 2 );
  stopifnot( .n.groups <= 4 );
  pairs.heights.map <- matrix( NA, nrow = length( levels( group ) ), ncol = length( levels( group ) ) );
  for( pair.low in 1:( length( levels( group ) ) - 1 ) ) {
      for( pair.high in ( pair.low + 1 ):length( levels( group ) ) ) {
          if( ( pair.high - pair.low ) == 1 ) {
              # Adjacent pairs go on the first level.
              pairs.heights.map[ pair.low, pair.high ] <- 1;
              pairs.heights.map[ pair.high, pair.low ] <- 1;
          } else if( ( pair.high - pair.low ) == 2 ) {
              # Need different levels for the overlapping lines 1,3 and 2,4.
              if( pair.low == 1 ) {
                  stopifnot( pair.high == 3 );
                  # (1,3) gets the second level.
                  pairs.heights.map[ pair.low, pair.high ] <- 2;
                  pairs.heights.map[ pair.high, pair.low ] <- 2;
              } else if( pair.low == 2 ) {
                  stopifnot( pair.high == 4 );
                  # (2,4) gets the third level.
                  pairs.heights.map[ pair.low, pair.high ] <- 3;
                  pairs.heights.map[ pair.high, pair.low ] <- 3;
              } else {
                  stop( "Unexpectedly, pair.low is neither 1 nor 2, but the span is 2" );
              }
          } else {
              stopifnot( pair.low == 1 );
              stopifnot( pair.high == 4 );
              # (1,4) gets the fourth level.
              pairs.heights.map[ pair.low, pair.high ] <- 4;
              pairs.heights.map[ pair.high, pair.low ] <- 4;
          }
      } # End foreach pair.high
  } # End foreach pair.low

  percentage.increment.for.first.comparison.bar <- 2.5;
  percentage.increment.per.comparison.bar <- 10;
  percentage.increment.height.of.comparison.bar <-
     ( percentage.increment.per.comparison.bar * .4 );
  percentage.increment.height.of.comparison.star <-
     ( percentage.increment.per.comparison.bar * .35 );
  .outcome.var.max <- base::max( the.var, na.rm = T );
  .outcome.var.min <- base::min( the.var, na.rm = T );
  
  .y.location.of.adjacent.group.comparison.lines <- ( .outcome.var.max * ( 1 + ( percentage.increment.for.first.comparison.bar / 100 ) ) ); # percentage.increment.per.comparison.bar% higher.
  .num.comparison.bar.heights <- base::max( pairs.heights.map, na.rm = TRUE );
  # Below, minus one because the first comparison line is at height .y.location.of.adjacent.group.comparison.lines. Plus one because we give some space to see the line.
  .y.max <- .y.location.of.adjacent.group.comparison.lines + ( ( ( .num.comparison.bar.heights - 1 ) + 1 ) * ( .outcome.var.max * ( percentage.increment.per.comparison.bar / 100 ) ) );
  .y.min <- base::min( 0, .outcome.var.min );
                                                    
  .ggp <- ggplot( .the.data.frame, aes(x=group.var,y=outcome.var)) + geom_boxplot(aes(fill=group.var) ) + guides(fill=FALSE) + ylim( .y.min, .y.max ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      axis.text.x =
          element_text(
              #angle = 45, vjust = 0.95, hjust=0.95,
              colour="black",size=14,face="bold"
           ),
      axis.text.y =
          element_text(
              colour="black",size=14,face="bold"
           )
        );
  if( the.var.scale.is.discrete ) {
      .ggp <- .ggp + scale_y_discrete();
  }
  .ggp <- .ggp + ylab( the.var.name ) + xlab( "" );
  if( !is.null( the.title ) && !is.na( the.title ) ) {
    .ggp <- .ggp +
        ggtitle( the.title ) + 
        theme( plot.title = element_text( lineheight=.8, face="bold" ) );
  }

  # Line start and end locations
  if( .n.groups == 2 ) {
    .x <- c(
          # Level 1: the adjacent pair
          c( 1, 1, .n.groups )
        );
    .xend <- c(
        # Level 1: the adjacent pair
          c( 1, rep( .n.groups, 2 ) )
        );
    .y <- c(
        # Level 1: the adjacent pair
        ( .y.location.of.adjacent.group.comparison.lines + c( 0, ( .outcome.var.max * ( percentage.increment.height.of.comparison.bar / 100 ) ), 0 ) )
        );
  } else {
    # 3 or 4 groups
    .x <- c(
          # Level 1: adjacent pairs
          c( 1, 1, sapply( 2:( .n.groups - 1 ), rep, 3 ), .n.groups ),
          # Level 2: ( 1, 3 )
          c( 1, 1, 3 )
        );
    .xend <- c(
        # Level 1: adjacent pairs
        c( 1, sapply( 2:( .n.groups - 1 ), rep, 3 ), rep( .n.groups, 2 ) ),
        # Level 2: ( 1, 3 )
        c( 1, 3, 3 )
        );
    .y <- c(
        # Level 1: adjacent pairs
        rep( ( .y.location.of.adjacent.group.comparison.lines + c( 0, ( .outcome.var.max * ( percentage.increment.height.of.comparison.bar / 100 ) ), 0 ) ), .n.groups - 1 ),
        # Level 2: ( 1, 3 )
        ( .y.location.of.adjacent.group.comparison.lines + ( ( .outcome.var.max * ( percentage.increment.per.comparison.bar / 100 ) ) * 1 ) + c( 0, ( .outcome.var.max * ( percentage.increment.height.of.comparison.bar / 100 ) ), 0 ) )
        );
    if( .n.groups == 4 ) {
      .x <- c( .x,
              # Level 3: ( 2, 4 )
              c( 2, 2, 4 ),
              # Level 4: ( 1, 4 )
              c( 1, 1, 4 )
      );
      .xend <- c( .xend,
                  # Level 3: ( 2, 4 )
                  c( 2, 4, 4 ),
                  # Level 4: ( 1, 4 )
                  c( 1, 4, 4 )
      );
      .y <- c( .y, 
            # Level 3: ( 2, 4 )
            ( .y.location.of.adjacent.group.comparison.lines + ( ( .outcome.var.max * ( percentage.increment.per.comparison.bar / 100 ) ) * 2 ) + c( 0, ( .outcome.var.max * ( percentage.increment.height.of.comparison.bar / 100 ) ), 0 ) ),
            # Level 4: ( 1, 4 )
            ( .y.location.of.adjacent.group.comparison.lines + ( ( .outcome.var.max * ( percentage.increment.per.comparison.bar / 100 ) ) * 3 ) + c( 0, ( .outcome.var.max * ( percentage.increment.height.of.comparison.bar / 100 ) ), 0 ) )
       );
    } else {
       stopifnot( .n.groups == 3 );
    }
  }
  
  # The lines end at the top (the max of each triple).
  .yend <- c( sapply( seq( 1, ( length( .y ) - 2 ), by = 3 ), function ( .triple.start ) { rep( base::max( .y[ .triple.start:( .triple.start + 2 ) ] ), 3 ) } ) );
  .xtext <- ( .x + ( ( .xend - .x ) / 2 ) )[ ( .xend - .x ) > 0 ];
  .ytext <- c( sapply( seq( 1, ( length( .y ) - 2 ), by = 3 ), function ( .triple.start ) { base::max( .y[ .triple.start:( .triple.start + 2 ) ] ) + ( .outcome.var.max * ( percentage.increment.height.of.comparison.star / 100 ) ) } ) );
  
  lines.df <- data.frame( x = .x, y = .y, xend = .xend, yend = .yend );
  astpos.df <- data.frame( x = .xtext, y = .ytext );
  
  ## Now remove entries corresponding to non-significant results.
  .include.it.p.value <- c( sapply( seq( 1, ( length( .x ) - 2 ), by = 3 ), function ( .triple.start ) { .tr <- test.results.as.list.of.lists[[ levels( group )[ .x[ .triple.start ] ] ]][[ levels( group )[ .x[ .triple.start + 2 ] ] ]]; if( ( length( .tr ) == 1 ) && is.na( .tr ) ) { .rv <- NA } else if( class( .tr ) == "htest" ) { .rv <- .tr$p.value; } else { .rv <- .tr }; return( rep( ifelse( is.null( .rv ), NA, .rv ), 3 ) ); } ) );
  .include.it <- ( !is.na( .include.it.p.value ) & ( .include.it.p.value <= p.value.threshold ) );
  if( sum( .include.it ) == 0 ) {
      return( .ggp );
  }
  
  .include.it.onepertriple <- .include.it[ seq( 1, ( length( .include.it ) - 1 ), by = 3 ) ];
                                                     
  stopifnot( nrow( lines.df ) == length( .include.it ) );
  lines.df.included <- lines.df[ .include.it, , drop = FALSE ];
  
  astpos.df.included <- astpos.df[ .include.it.onepertriple, , drop = FALSE ];
                                                     
  ast.text.p.values <- .include.it.p.value[ .include.it ];
  if( length( ast.text.p.values ) == 0 ) {
      return( .ggp );
  }
  
  # But each one is there three times, so take just the first in each triple.
  if( length( ast.text.p.values ) > 0 ) {
      ast.text.p.values <- ast.text.p.values[ seq( 1, ( length( ast.text.p.values ) - 1 ), by = 3 ) ];
  }
  ast.text <- sapply( ast.text.p.values, prettyPrintPValuesTo4Digits );
  ast.text <- sapply( ast.text, function( .p.txt ) { if( substring( .p.txt, 1, 1 ) == "<" ) { paste( "P <", substring( .p.txt, 2 ) ) } else { paste( "P =", .p.txt ) } } );
  stopifnot( sum( .include.it.onepertriple ) == length( ast.text ) );

  if( nrow( lines.df.included ) > 0 ) {
      .ggp2 <- .ggp + geom_segment(data = lines.df.included, size = .5, aes(x=x, y=y, xend=xend, yend=yend))
  } else {
      .ggp2 <- .ggp;
  }
  if( length( ast.text ) > 0 ) {
      .ggp3 <- .ggp2 + geom_text(data = astpos.df.included, aes(x=x, y=y, fontface = "bold.italic"), label=ast.text, size = 4, color = "red" );
  } else {
      .ggp3 <- .ggp2;
  }
  
  return( .ggp3 );
} # createBoxplotShowingSignificance (..)

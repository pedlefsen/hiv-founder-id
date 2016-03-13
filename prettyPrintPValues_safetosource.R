createPrettyPrintPValuesToXDigits <- function( X ) {
  stopifnot( X > 0 );
  return( function( p.value, equalsText = "", lessThanText = "<"  ) { # Note that there's two chars for "0." in "0.05"...
      if( is.na( p.value ) | !is.finite( as.numeric( p.value ) ) ) { # Could be -Inf if from say max( c( NA, NA ) )
        return( "" );
      }
      p.value <- as.numeric( p.value );
      ## TODO: REMOVE
      if( p.value < 0 ) {
          print( paste( "UH OH <0", p.value ) );
      }
      stopifnot( p.value >= 0 );
      stopifnot( p.value <= 1 );
      p.value <- sprintf( paste( "%1.", X, "f", sep = "" ), p.value );
      if( p.value == sprintf( paste( "%1.", X, "f", sep = "" ), 0 ) ) {
          if( X > 1 ) {
              p.value <- paste( lessThanText, "0.", paste( rep( "0", X-1 ), collapse = "" ), "1", sep = "" );
          } else {
              p.value <- paste( lessThanText, "0.1", sep = "" );
          }
      } else if( equalsText != "" ) {
          p.value <- paste( equalsText, p.value, sep = "" );
      }
      return( p.value );
     } # <anonymous fn>( p.value )
    ); # return(..)
} # createPrettyPrintPValuesToXDigits ( X )( p.value )
prettyPrintPValuesTo4Digits <- createPrettyPrintPValuesToXDigits( 4 );
prettyPrintPValuesTo2Digits <- createPrettyPrintPValuesToXDigits( 2 );


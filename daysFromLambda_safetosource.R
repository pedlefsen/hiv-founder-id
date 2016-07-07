    ## This is originally from PFitter, but daysFromLambda has been modified to take inverse.lambda instead of lambda.  Epsilon is the per position mutation rate, per generation -- but I think that the generations / day is fixed at around 0.5, so really this is the mutation rate position per two days. Note this lambda is PFitter's (so it is the lambda paramater of the poisson process describing the intersequence HD distribution).
    default.epsilon <- 2.16e-05;
    phi <- sqrt( 1 + 4/3 );
    phi.scalar <- 1.5 * ( phi / ( 1 + phi ) );
    daysFromLambda.coefficient.of.inverse.epsilon <- function ( lambda, nb ) {
        phi.scalar * ( lambda / nb )
    }
    daysFromLambda.constant <- ( phi.scalar * ( - (1-phi)/(phi^2) ) );
    daysFromLambda <- function ( lambda, nb, inverse.epsilon = (1/default.epsilon) ) {
        daysFromLambda.constant + daysFromLambda.coefficient.of.inverse.epsilon( lambda, nb ) * inverse.epsilon;
    }


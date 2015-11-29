### R code from vignette source '/Users/pedlefsen/src/from-git/hiv-founder-id/DSPFitter.Rnw'
### Encoding: ASCII

###################################################
### code chunk number 1: DSPFitter.Rnw:18-29
###################################################
## R packages needed
library( "xtable" )
library( "coin" )
library( "ggplot2" )
library( "binom" )

## for prettyPrintPValuesTo4Digits
#source( "~/src/from-git/projects/pedlefse/rapporttemplates/rapporttemplates-util.R" )

# Setup for prettier Sweave output.
old.continue.option <- options( continue = " " )


###################################################
### code chunk number 2: DSPFitter.Rnw:32-161
###################################################
#### ERE I AM, there is still a problem with the plausibility calculation being unreasonably high when we get into the twilight zone where numeric stuff fails.  Particularly the problem is in logSubtract returning -Inf when its arguments are different but very close.  I think for now we can just sometimes give up, and that's ok.  It seems to happen when the data are very weird, as in the cases that violate the poisson model and assumed relationship between consensus and intersequence distances.  BUT TO CORRECT THIS, NEED NUMERICAL STABILITY MAYBE ONLY AVAILABLE WITH BFLOAT.

    # MAGIC N: DS.NDRAWS
    #EPSILON.NDRAWS <- 1000;
    EPSILON.NDRAWS <- 100;
    DS.NDRAWS <- 1000;
    

  # The trick for subtracting small logged doubles in log space, with logX > logY
  logSubtract <- Vectorize( function (logX, logY) { # modified from: https://facwiki.cs.byu.edu/nlp/index.php/Log_Domain_Computations
      #print( paste( "logSubtract", logX, logY ) );
       if( !is.finite( logX ) ) {
           return( logX );
       }
       negDiff <- logY - logX;
       # exp.negDiff <- exp( negDiff );
       # if( exp.negDiff == 1 ) {
       #     negDiff <- -1E-07; # Smallest value on my computer st exp( negDiff ) < 1
       #     #return( logX );
       # }
       # 4. otherwise use some nice algebra to stay in the log domain
       #    (except for negDiff)
      .rv <- logX + log1p( - exp( negDiff ) );
      #print( paste( "RV:", .rv ) );
       return( .rv );
   } ); # logSubtract (modified from https://facwiki.cs.byu.edu/nlp/index.php/Log_Domain_Computations)
PoissonDSM.calculatePlausibilityOfSingleton.one.obs <- Vectorize(
    function( x, k, log.result = FALSE ) {
        .rv <- pgamma( x, shape = k ) - pgamma( x, shape = k + 1 );
        if( ( .rv == 0 ) || log.result ) {
            .rv <- log( .rv );
            if( any( !is.finite( .rv ) ) ) {
                .rv.log.part1 <- pgamma( x, shape = k, log.p = TRUE );
                .rv.log.part2 <- pgamma( x, shape = ( k + 1 ), log.p = TRUE );
                .rv <- logSubtract( .rv.log.part1, .rv.log.part2 );
                #print( paste( "exp( ", .rv.log.part1, " ) - exp( ", .rv.log.part2, " ) = exp( ", .rv, " )", sep = "" ) );
                # NOTE this will still include non-finite elements when x is very large. 
                if( any( !is.finite( .rv ) ) ) {
                    .rv.hack <- 
                    # Very special hack:
                     # when k is 0, we get underflow that can be recovered using the special trick that logSubtract( 0, pgamma( x, shape = k, log.p = TRUE ) ) is approximately exp( -x - 2*k ).
                     ( -x - 2*k );
                }
                .rv <- ifelse( is.finite( .rv ), .rv, .rv.hack );
            }
            if( !log.result ) {
                .rv <- exp( .rv );
            }
        }
        return( .rv );
    }, vectorize.args = "k" );
PoissonDSM.calculatePlausibilityOfSingletons <- function ( xs, observed.counts, log.result = FALSE, scale.result.by.log = 0, finite.sums.coefs = NULL ) {
    .mat <-
        PoissonDSM.calculatePlausibilityOfSingleton.one.obs( xs, observed.counts, log.result = log.result );
    if( is.null( dim( .mat ) ) ) {
        .mat <- t( as.matrix( .mat ) );
    }
    if( log.result ) {
        .sums <- apply( .mat, 1, sum, na.rm = T );
        if( any( !is.finite( .sums ) ) || ( !is.null( finite.sums.coefs ) && ( length( finite.sums.coefs ) < 2 ) ) ) {
          ## Fill in the nonfinite ones by extrapolation.  Works nearly perfectly here..  Or if the coefs were figured out on a previous iteration, use them...
            .nonfinite.xs <- xs[ !is.finite( .sums ) ];
          if( is.null( finite.sums.coefs ) || ( length( finite.sums.coefs ) < 2 ) ) {
            .finite.sums <- .sums[ is.finite( .sums ) ];
            .finite.xs <- xs[ is.finite( .sums ) ];
            finite.sums.coefs <- coef( lm( .finite.sums  ~ .finite.xs ) );
          }
          if( any( !is.finite( .sums ) ) ) {
              .nonfinite.sums <- finite.sums.coefs %*% sapply( .nonfinite.xs, function( .nonfinite.x ) { c( 1, .nonfinite.x ) } );
              .sums[ !is.finite( .sums ) ] <- .nonfinite.sums;
          }
        } # Calculate slope; fix non-finite sums by extrapolation
    
        .log.plaus <- .sums;
        names( .log.plaus ) <- xs;
        if( is.numeric( scale.result.by.log ) && ( scale.result.by.log != 0 ) ) { .log.plaus <- ( .log.plaus + scale.result.by.log ); }
        attr( .log.plaus, "finite.sums.coefs" ) <- finite.sums.coefs;
        return( .log.plaus );
    } else {
        .plaus <- apply( .mat, 1, prod );
        names( .plaus ) <- xs;
        if( is.numeric( scale.result.by.log ) && ( scale.result.by.log != 0 ) ) { .log.plaus <- ( .log.plaus * exp( scale.result.by.log ) ); }
        return( .plaus );
    }
} # PoissonDSM.calculatePlausibilityOfSingletons (..)
    
# But now we want the hpd
PoissonDSM.calculatePlausibilityOfSingletons.integratedPlausibility <- function( .dta, log.result = FALSE ) {
    if( log.result ) {
        # First get a sense of what the scale factor should be; look at the peak.
        .scale.factor = 0 - suppressWarnings( optimize( f = function( .x ) { PoissonDSM.calculatePlausibilityOfSingletons( .x, .dta, log.result = TRUE, scale.result.by.log = 0 ) }, lower = 0, upper = mean( .dta ), maximum = T ) )$objective;
        print( paste( "SCALE FACTOR:", .scale.factor ) );
        # Actually to get a slope for interpolating the larger values, calc a few..
        ## finite.sums.coefs <- NULL;
        ## lower.end.of.range.of.values.to.try <- 100;
        ## while( is.null( finite.sums.coefs ) ) { 
        ##   .result.mostly.ignored <-
        ##       PoissonDSM.calculatePlausibilityOfSingletons( lower.end.of.range.of.values.to.try + 0:99, .dta, log.result = TRUE );
        ##   if( !is.null( attr( .result.mostly.ignored, "finite.sums.coefs" ) ) ) {
        ##       finite.sums.coefs <- attr( .result.mostly.ignored, "finite.sums.coefs" );
        ##       ## TODO: REMOVE
        ##       print( "FINITE SUMS COEFS" );
        ##       print( finite.sums.coefs );
        ##       break;
        ##   } else {
        ##       lower.end.of.range.of.values.to.try <- lower.end.of.range.of.values.to.try + 100;
        ##   }
        ##   #print( lower.end.of.range.of.values.to.try );
        ## } # End while searching for breakdown where finite.sums.coefs is best computed.
        .integrated.fn <- function( x ) {
            .log.rv <- PoissonDSM.calculatePlausibilityOfSingletons( x, .dta, log.result = TRUE, scale.result.by.log = .scale.factor );#, finite.sums.coefs = finite.sums.coefs );
            #print( paste( x[1], .dta[1], exp( .log.rv[1] ) ) );
            exp( .log.rv );
        };
        .foo <- function( y ) {
            tryCatch(
              {
                unname( log( ( integrate( f = .integrated.fn, lower = 0, upper = y, subdivisions = 1000 )$value ) ) - .scale.factor );
              }, error = function ( .e.msg ) {
                warning( paste( "Integration proglems in PoissonDSM.calculatePlausibilityOfSingletons.integratedPlausibility:", .e.msg ) );
                0;
              }
            );
       } # .foo
        Vectorize( .foo );
    } else {
        Vectorize( function( y ) { unname( integrate( f = function( x ) { PoissonDSM.calculatePlausibilityOfSingletons( x, .dta ) }, lower = 0, upper = y )$value ) } ) 
    }
} # PoissonDSM.calculatePlausibilityOfSingletons.integratedPlausibility (..)


###################################################
### code chunk number 3: DSPFitter.Rnw:164-187
###################################################
################################################################################## ten channel data generation
# ten.channel.gen.count <- 100;
# ten.channel.gen.truth.lambdas <- 1:ten.channel.gen.count;
# ten.channel.gen.outrageouslyhighvalue <- ten.channel.gen.truth.lambdas + 4 * sqrt( ten.channel.gen.truth.lambdas ); # The error increases a lot if this is too high.  I think it should be no more than 20 times the true lambda.  TODO: Detect it somehow.
# 
# ten.channel.gen.n <- t( apply( as.array( 1:ten.channel.gen.count ), 1, function( i ) { rpois( 10, ten.channel.gen.truth.lambdas[ i ] ) } ) );
# rownames( ten.channel.gen.n ) <- ten.channel.gen.truth.lambdas;
# 
# # estimates using the mode of the plausibility transform are the same as the naive estimate (mean):
# # First 26 use lower upper bound for search or else it gets lost
# stopifnot( all.equal( apply( ten.channel.gen.n[1:26, ], 1, function( .dta ) { optimize( function( x ) { PoissonDSM.calculatePlausibilityOfSingletons( x, .dta ) }, lower = 0, upper = 40, maximum = T, tol = 1E-5 )$maximum } ), apply( ten.channel.gen.n[1:26, ], 1, function( .dta ) { mean( .dta ) } ), tolerance = 1E-5 ) )
# # Now the rest
# stopifnot( all.equal( apply( ten.channel.gen.n[27:100, ], 1, function( .dta ) { optimize( function( x ) { PoissonDSM.calculatePlausibilityOfSingletons( x, .dta ) }, lower = 0, upper = 200, maximum = T, tol = 1E-5 )$maximum } ), apply( ten.channel.gen.n[27:100, ], 1, function( .dta ) { mean( .dta ) } )  ) )

# ten.channel.gen.normalizedPlausibilities <- function ( k ) {
#     ten.channel.gen.integratedPlausibility.totalplaus <- PoissonDSM.calculatePlausibilityOfSingletons.integratedPlausibility( ten.channel.gen.n[ k, ] )( ten.channel.gen.outrageouslyhighvalue[ k ] );
#     return( Vectorize( function( y ) { PoissonDSM.calculatePlausibilityOfSingletons.integratedPlausibility( ten.channel.gen.n[ k, ] )( y ) / ten.channel.gen.integratedPlausibility.totalplaus } ) );
# }
# ##  NOTE HACK using upper = ten.channel.gen.outrageouslyhighvalue
# ten.channel.gen.quantile <- function( k ) { function( quantile ) { optimize( function( y ) { abs( ten.channel.gen.normalizedPlausibilities( k )( y ) - quantile ) }, lower = 0, upper = ten.channel.gen.outrageouslyhighvalue[ k ] )$minimum } }
# 
# ten.channel.gen.quantiles.by.lambda <- t( sapply( truth.lambdas, function( k ) { c( ten.channel.gen.quantile( k )( 0.025 ), ten.channel.gen.quantile( k )( 0.975 ) ) } ) );
# colnames( ten.channel.gen.quantiles.by.lambda ) <- c( "2.5%", "97.5%" );


###################################################
### code chunk number 4: DSPFitter.Rnw:190-226
###################################################
####==== From PFitter.R

library(tools)
args <- commandArgs(trailingOnly=TRUE);
if( length( args ) == 0 ) {
  args <- 
    #c( "newest-Abrahams-2009aa-hiv-founder-id_resultDir/1172_removeHypermutatedSequences_PoissonFitterDir/1172_removeHypermutatedSequences_pairwiseHammingDistances.txt",
      c( "rv217_1W_gold_standard-hiv-founder-id_-fs_resultDir/rv217_anon_1w_01_fixHypermutatedSequences_PoissonFitterDir/rv217_anon_1w_01_fixHypermutatedSequences_pairwiseHammingDistances.txt",
#    c( "~/src/from-git/hiv-founder-id/Abrahams-2009aa-hiv-founder-id_resultDir/0334_fixHypermutatedSequences_PoissonFitterDir/0334_fixHypermutatedSequences_pairwiseHammingDistances.txt",
  #    c( "~/src/from-git/hiv-founder-id/Abrahams-2009aa-hiv-founder-id_resultDir/0478_fixHypermutatedSequences_removeRecombinedSequences_PoissonFitterDir/0478_fixHypermutatedSequences_removeRecombinedSequences_pairwiseHammingDistances.txt",
     "2.16e-05", "8814" );#"2517" );
}


# sample <- basename(file_path_sans_ext(infile))

#dir <- paste(dirname(infile), '/', sep='')
#outfile <- paste(dir, "DS.LOG_LIKELIHOOD.results.txt", sep="")
#outfile2 <- paste(dir, "DS.CONVOLUTION.results.txt", sep="")

cat( "DS PFitter.." );

### FUNCTIONS ###
iseven <- function(c) {
	c1 <- c/2-as.integer(c/2)
	if(c1==0){ ev <- TRUE } else { ev <- FALSE }
	return(ev)
}

phi <- sqrt(1+4/3)
days <- function(l,nb,epsilon) 1.5*((phi)/(1+phi))*(l/(epsilon*nb) - (1-phi)/(phi^2))
### end FUNCTIONS ###

#if(!file.exists(outfile)){
#	write(paste("Sample", "Lambda", "St.Dev", "NSeq", "NBases", "MeanHD", "MaxHD","Days(CI)", "Chi2","DF","Goodness_of_pval", sep="\t"), file=outfile, append=FALSE)
#}


###################################################
### code chunk number 5: DSPFitter.Rnw:292-1270
###################################################
PFitter <- function (
  infile = args[1],
  dlist = read.table( file=infile, sep="\t", stringsAsFactors=F ),
  epsilon = c(as.numeric(args[2])),
  nbases = c(as.numeric(args[3])),
  be.verbose = TRUE,
  estimate.lambda.and.days = TRUE,
  do.chisq.gof = TRUE
) {
    d0 <- dlist[which(dlist[,1]==dlist[1,1]),]
    d1 <- dlist[-which(dlist[,1]==dlist[1,1]),]	
    
    consensus.distances <- d0[ , 3 ];
    intersequence.distances <- d1[ , 3 ];
    
    mult0 <- as.numeric(sub('.+_(\\d+)$', '\\1', d0[,2]));
    nseq <- sum(mult0)
    
    seqnames <- unique(c( d1[,1], d1[,2] ));
    
    yvec <- rep(0, (1+max(intersequence.distances)))
    for(i in 1:length(seqnames)){
        tmp <- d1[which(d1[,1]==seqnames[i]),,]
        if( nrow( tmp ) == 0 ) {
            next;
        }
    	m0 <- as.numeric(sub('.+_(\\d+)$', '\\1', tmp[1,1]))
    	yvec[1] <- yvec[1] + 0.5*m0*(m0-1) ## 0 bin
    	for(j in 1:dim(tmp)[1]){
    		m1 <- as.numeric(sub('.+_(\\d+)$', '\\1', tmp[j,2]))
    		val <- tmp[j,3]
    		yvec[val+1] <- yvec[val+1] + m0*m1
    	}
    }
    
    nl <- length(yvec)
    lambda <- sum((1:(nl-1))*yvec[-1])/sum(yvec);

    if( estimate.lambda.and.days ) {
      ### Fitting
      ## U STAT from PFitter.R
      calculateUSTATvarHD <- function () {
        #### U STAT ESTIMATE OF ST DEV
        #### FORMULAE
        #### Var(HD) = (N(N-1)/2)^(-1) (2(N-2)sigma1^2 + sigma2^2)
        #### sigma1^2 = (N(N-1)(N-2)/3 -1)^(-1) sum_{i<j<l} ((Dij-mu)(Dil-mu)+(Dij-mu)(Djl-mu))
        #### sigma2^2 = (N(N-1)/2-1)^(-1) sum_{i<j} (Dij-mu)^2
        
        ### construct a matrix of Dij's
        TX <- matrix(rep(NA,nseq^2), ncol=nseq)
        rownames( TX ) <- seqnames;
        colnames( TX ) <- seqnames;
        for(i in 1:(dim(d0)[1]-1)){
        	useq <- d0[i,2]
        	TX[d1[which(d1[,1]==useq),2],i] <- d1[which(d1[,1]==useq),3];
        }
        
        sigma1 <- 0
        sigma2 <- 0
        muhat <- 0
        denmu <- (sum( !is.na( TX )))^(-1)
        ## TODO: Figure out what (if any) is the right fix to the below to handle sparse distances
        den1 <- 12*(nseq*(nseq-1)*(nseq-2)*(nseq-3))^(-1)  
        den2 <- den1/4
        
        for(n in 1:(nseq-1)){
            for(m in (n+1):nseq){
                if( !is.na( TX[ m, n ] ) ) {
                    muhat <- muhat + mult0[n]*mult0[m]*denmu*TX[m,n]
        	}
            }
        }
        
        for(n in 1:nseq){
        	dnn <- 0
        	sigma1 <- sigma1 + choose(mult0[n],3)*den1*2*(dnn-muhat)^2 
        	sigma2 <- sigma2 + choose(mult0[n],2)*den2*(dnn-muhat)^2
        	if(n != nseq){
        		for(m in (n+1):nseq){
                            dnm <- TX[m,n]
                            if( !is.na( dnm ) ) {
        			dmm <- 0
        			sigma2 <- sigma2 + mult0[n]*mult0[m]*(dnm - muhat)^2
        			sigma1 <- sigma1 + (2/3)*choose(mult0[n],2)*mult0[m]*(dnm-muhat)*(dnm+2*dnn-3*muhat)
        			sigma1 <- sigma1 + (2/3)*mult0[n]*choose(mult0[m],2)*(dnm-muhat)*(dnm+2*dmm-3*muhat)
        			if(m != nseq){
        				for(l in (m+1):nseq){
        					dnl <- TX[l,n]
        					dlm <- TX[l,m]
                                                if( !is.na( dnl ) && !is.na( dlm ) ) {
                                                    sigma1 <- sigma1 + (2/3)*mult0[n]*mult0[m]*mult0[l]*((dnm-muhat)*(dnl-muhat)+(dnm-muhat)*(dlm-muhat)+(dnl-muhat)*(dlm-muhat))
                                                }
        				}
        			}
        		    } # End if !is.na( dnm )
                        } # End for( m )
         	}
        }
        
        ## varhd <- sqrt(denmu*(2*(nseq-2)*sigma1 + sigma2))
        A <- 8/(nseq*(nseq-1)*(nseq-2)*(nseq-3))
        B <- 4/(nseq*(nseq-1)*(nseq-2)*(nseq-3))
        newvarhd <- sqrt(A*sigma1 + B*sigma2)
        return( newvarhd );
      } # calculateUSTATvarHD (..)
      newvarhd <- calculateUSTATvarHD();
      lambda.low <- lambda - 1.96*newvarhd;
      lambda.high <- lambda + 1.96*newvarhd;
      
      if( be.verbose ) {
          cat(paste("PFitter Estimated Lambda is ", format(lambda, digits=4 ), " (95% CI ", format(lambda.low, digits=4 ), " to ", format(lambda.high, digits=4 ), ")", sep=""), fill = TRUE )
      }
      
      estdays <- days(lambda, nbases, epsilon)
      upplim <- days(lambda.high, nbases, epsilon)
      lowlim <- days(lambda.low, nbases, epsilon)
      
      pfitter.results <- list( lambda = c( "estimate" = lambda, "2.5%" = lambda.low, "97.5%" = lambda.high ), days = c( "estimate" = estdays, "2.5%" = lowlim, "97.5%" = upplim ) );
      if( be.verbose ) {
        uppdays <- round(upplim)
        lowdays <- round(lowlim)
        formatteddays <-
            paste(round(estdays), " (", lowdays, ", ", uppdays, ")", sep="") 
        cat(paste("PFitter Estimated Days:", formatteddays, sep=" "), fill = TRUE )
      }
    } else { # if estimate.lambda.and.days .. else ..
      pfitter.results <- list();
    } # End if estimate.lambda.and.days .. else ..

    if( do.chisq.gof ) {
      yvec0 <- rep(0, (1+max(consensus.distances)))
      for(i in 1:(1+max(consensus.distances))){ yvec0[i] <- sum(mult0[which(consensus.distances==(i-1))]) }
      nl0 <- length(yvec0);
      clambda <- sum((1:(nl0-1))*yvec0[-1])/sum(yvec0) #### THIS IS THE LAMBDA THAT FITS THE CONSENSUS ONLY DISTRIBUTION
      ### CONSTRUCT SIGMA_ij MATRIX THEN INVERT IT
      #pk <- function(x) ((nseq^2)*(2^x)*exp(-2*clambda)*(clambda^x))/factorial(x)
      pk <- function(x) (exp( ( (log(nseq)*2)+((log(2*clambda))*x)+(-2*clambda))-lfactorial(x) ) )
      mui <- function(x) nseq*dpois(x, lambda=clambda)
      SIGMA.DIM.MAX <- 170; # Beyond this value, factorial stops working in R.
      sigma.dim <- min( SIGMA.DIM.MAX, (2*nl0) );
      eyvec <- 0.5*pk(0:(sigma.dim-1))
      eyvec[ !is.finite( eyvec ) ] <- 0; # Paul added this to avoid errors.  factorial(x) sometimes refuses to work (if x is too large).
      
      dvec1 <- 0
      for(i in 1:length(yvec)){ dvec1 <- c(dvec1, rep((i-1),yvec[i])) }
      dvec1 <- dvec1[-1]
      
      meanhd <- mean(dvec1)
      maxhd <- max(dvec1)
      
      if (lambda!=0) {
      	sigmaij <- matrix(nrow=sigma.dim, ncol=sigma.dim)
      	coeff <- (nseq^3)*exp(-3*clambda)
      	
      	#### RICORDATI!!!! EYVEC[K] == E(Y_{K-1}) !!!!!!
      	
      	for(k in 0:(sigma.dim-1)){  
        
      		for(l in 0:(sigma.dim-1)){   
      			
      			if(k>=l){ 
      				c1 <- ((clambda^k)/factorial(k))*sum(choose(k,l:0)*((clambda^(0:l))/factorial(0:l)))
      				c2 <- ((clambda^l)/factorial(l))*sum(choose(l,l:0)*((clambda^((k-l):k))/factorial((k-l):k))) 
      			}
      			if(k<l){
      				c1 <- ((clambda^l)/factorial(l))*sum(choose(l,k:0)*((clambda^(0:k))/factorial(0:k)))
      				c2 <- ((clambda^k)/factorial(k))*sum(choose(k,k:0)*((clambda^((l-k):l))/factorial((l-k):l))) 
      			}
                          if( is.na( c1 ) ) {
                              c1 <- 0;
                          }
                          if( is.na( c2 ) ) {
                              c2 <- 0;
                          }
      			sigmaij[k+1,l+1] <- 0.5*coeff*(c1+c2)
      
      			
      			if(k==l){ sigmaij[k+1,l+1] <- sigmaij[k+1,l+1] + (0.5)*pk(k) }
      			if((k==l)&(iseven(k))){ sigmaij[k+1,l+1] <- sigmaij[k+1,l+1] - (0.25)*mui(k/2) }
      		}
      	}
          
              #print( "about to run La.svd" );
      	sdec <- La.svd(sigmaij)
              #print( "ran La.svd" );
      	diag <- ifelse(sdec$d>1e-4,sdec$d,0)
      	diagmat <- matrix(rep(0,sigma.dim^2), ncol=sigma.dim)
      	for(ii in 1:sigma.dim){diagmat[ii,ii]<-ifelse(diag[ii]==0,0,1/diag[ii])}
      	sigmainv <- sdec$u%*%diagmat%*%sdec$vt
      	
      	h <- hist(dvec1[ dvec1 <= sigma.dim ], breaks=seq(-1,(sigma.dim-1),1), plot=FALSE)
      	xvec <- h$breaks
      	yvec <- h$counts
      	#nl1 <- length(yvec)
      	#aplambda <- sum((1:(nl1-1))*yvec[-1])/sum(yvec)
      			
      	#pesce <- 0.5*nseq*(nseq-1)*dpois(0:(nl1-1), lambda=aplambda)
      
      	if (length(yvec)<sigma.dim) { 
      		ccvv <- sigma.dim - length(yvec) 
      		yvec <- c(yvec, rep(0,ccvv)) 
      	}
      	
      	chisq <- ( t(abs(yvec-eyvec))%*%sigmainv%*%(abs(yvec-eyvec)) )[ 1, 1 ];
              pval <- ifelse(chisq<0,2e-16,1-pchisq(chisq,df=nl0-1))
      	if(pval==0){ pval <- 2e-16 }
      	if(chisq<0){ chisq <- NA }
      } else { 
      	chisq <- NA
      	nl0 <- NA # Paul fixed to "nl0" from "nl"
      	pval <- 0
      }
      
      if( be.verbose ) {
          .message <- paste( "The stats for the PFitter Chi-squared test for whether our intersequence counts deviate from the expected values of a Poisson with lambda equal to twice the estimate from HD to the consensus sequence are chisq = ", round( chisq, digits = 4 ), ", P = ", round( pval, digits = 4 ), ".", sep = "" );
          if( pval < 0.05 ) {
              cat( paste( "DOES NOT FIT.", .message, sep = "" ), fill = TRUE );
          } else {
              cat( paste( "FITS (insufficient evidence to conclude that it does not).\n", .message, sep = "" ), fill = TRUE );
          }
      }
      pfitter.results <- c( pfitter.results, list( twox = list( "chisq.stat" = chisq, "chisq.df" = ( nl0 - 1 ), "P" = pval ) ) );    
    } # End if do.chisq.gof

    return( pfitter.results );
} # PFitter (..)

calculateSumOfDistancesAccountingForSequenceMultiplicity <- function ( intersequence.dlist ) {
    .mult.1 <- suppressWarnings( as.numeric( gsub( "^.*_(\\d+)$", "\\1", intersequence.dlist[ , 1 ] ) ) );
    .mult.1[ is.na( .mult.1 ) | ( .mult.1 == 0 ) ] <- 1; # if it is missing the _nnn, just call it 1.
    .mult.2 <- suppressWarnings( as.numeric( gsub( "^.*_(\\d+)$", "\\1", intersequence.dlist[ , 2 ] ) ) );

    .mult.2[ is.na( .mult.2 ) | ( .mult.2 == 0 ) ] <- 1; # if it is missing the _nnn, just call it 1.
    return( list( sum = sum( intersequence.dlist[ , 3 ] * .mult.1 * .mult.2, na.rm = T ), count = sum( .mult.1 * .mult.2, na.rm = T ) ) ); 
} # calculateSumOfDistancesAccountingForSequenceMultiplicity (..)

replicateDistancesForSequenceMultiplicity <- function ( any.dlist, missing.seqnames = NULL ) {
    # Note that this isn't the entire thing because we also have to add the "0" distance distances among the duplicated seqs; see below.
    #print( missing.seqnames );
    if( is.null( any.dlist ) || is.null( dim( any.dlist ) ) || ( nrow( any.dlist ) == 0 ) ) {
        maybe.longer.dlist.list <- NULL;
    } else {
      maybe.longer.dlist.list <-
      lapply( 1:nrow( any.dlist ), function( .row.i ) {
        .row <- any.dlist[ .row.i, ];
        .mult.1 <- NA;
        if( length( grep( "^.*_(\\d+)$", .row[ 1 ], perl = TRUE ) ) > 0 ) {
            .mult.1 <- as.numeric( gsub( "^.*_(\\d+)$", "\\1", .row[ 1 ], perl = TRUE ) );
        }
        if( is.na( .mult.1 ) ) {
            .mult.1 <- 1; # if it is missing the _nnn, just call it 1.
        }
        .mult.2 <- NA;
        if( length( grep( "^.*_(\\d+)$", .row[ 2 ], perl = TRUE ) ) > 0 ) {
            .mult.2 <- as.numeric( gsub( "^.*_(\\d+)$", "\\1", .row[ 2 ], perl = TRUE ) );
        }
        if( is.na( .mult.2 ) ) {
            .mult.2 <- 1; # if it is missing the _nnn, just call it 1.
        }
        if( ( .mult.1 * .mult.2 ) > 1 ) {
          .rv <- matrix( "", ncol = 3, nrow = ( .mult.1 * .mult.2 ) );
          .rv[ , 1 ] <- .row[[ 1 ]];
          .rv[ , 2 ] <- .row[[ 2 ]];
          .rv[ , 3 ] <- .row[[ 3 ]];
          return( .rv );
      } else {
          .rv <- matrix( .row, nrow = 1 );
          mode( .rv ) <- "character";
          return( .rv );
      }
    } );
    }
    if( !is.null( maybe.longer.dlist.list ) && ( length( maybe.longer.dlist.list ) > 0 ) ) {
      maybe.longer.dlist <- do.call( rbind, maybe.longer.dlist.list );
      colnames( maybe.longer.dlist ) <- colnames( any.dlist );
    } else {
      maybe.longer.dlist <- NULL;
    }
    duplicated.entries.dlist.list <-
        lapply( unique( c( any.dlist[ , 1 ], any.dlist[ , 2 ], missing.seqnames ) ), function( .seq.name ) {
            .mult <- NA;
            if( length( grep( "^.*_(\\d+)$", .seq.name, perl = TRUE ) ) > 0 ) {
                .mult <- as.numeric( gsub( "^.*_(\\d+)$", "\\1", .seq.name, perl = TRUE ) );
            }
            if( is.na( .mult ) ) {
                .mult <- 1; # if it is missing the _nnn, just call it 1.
            }
            if( .mult == 1 ) {
                return();
            }
            .rv <- matrix( "", ncol = 3, nrow = choose( .mult, 2 ) );
            .rv[ , 1 ] <- .seq.name;
            .rv[ , 2 ] <- .seq.name;
            .rv[ , 3 ] <- 0;
            return( .rv );
    } );
    if( !all( unlist( lapply( duplicated.entries.dlist.list, is.null ) ) ) ) {
        duplicated.entries.dlist <- do.call( rbind, duplicated.entries.dlist.list );
        if( is.null( dim( duplicated.entries.dlist ) ) ) {
            duplicated.entries.dlist <- matrix( duplicated.entries.dlist, nrow = 1, ncol = 3 );
            mode( duplicated.entries.dlist ) <- "character";
        }
        colnames( duplicated.entries.dlist ) <- colnames( any.dlist );
        if( is.null( maybe.longer.dlist ) ) {
          return( duplicated.entries.dlist );
        } else {
          return( rbind( duplicated.entries.dlist, maybe.longer.dlist ) );
        }
    } else {
        return( maybe.longer.dlist );
    }
    # Put zero entries first to maintian sort order jic.
} # replicateDistancesForSequenceMultiplicity (..)

BayesPFitter <- function (
  infile = args[1],
  dlist = read.table( file=infile, sep="\t", stringsAsFactors=F ),
  epsilon = c(as.numeric(args[2])),
  nbases = c(as.numeric(args[3])),
  be.verbose = TRUE
) {
    ## Bayesian, too.  Congjugate prior is Gamma, and Gamma( 1, 1 ) (expo) has log( x ) propto 1.  [However note that alpha and beta approaching zero might give the closest analogue to the DS solution.]
    # Here I'm using the posterior median; the mode or mean are other options; can be looked up for generic gamma distrns.
    # Using posterior median, just to be different.
    prior.gamma.shape <- 1;
    prior.gamma.rate <- 1;
    .id.sum.and.count <- calculateSumOfDistancesAccountingForSequenceMultiplicity( dlist[ -which(dlist[,1]==dlist[1,1]), ] );
    .id.sum <- .id.sum.and.count[[ "sum" ]];
    .id.count  <- .id.sum.and.count[[ "count" ]];
    posterior.lambda <- qgamma( .5, prior.gamma.shape + .id.sum, prior.gamma.rate + .id.count );
    posterior.lambda.low <- qgamma( .025, prior.gamma.shape + .id.sum, prior.gamma.rate + .id.count );
    posterior.lambda.high <- qgamma( .975, prior.gamma.shape + .id.sum, prior.gamma.rate + .id.count );
    posterior.estdays <- days(posterior.lambda, nbases, epsilon)
    posterior.upplim <- days(posterior.lambda.high, nbases, epsilon)
    posterior.lowlim <- days(posterior.lambda.low, nbases, epsilon)
    
    posterior.lowdays <- round(posterior.lowlim)
    posterior.uppdays <- round(posterior.upplim)
    
    if( be.verbose ) {
        cat( paste("Bayesian PFitter Estimated Lambda is ", format(posterior.lambda, digits=4), " (95% CI ", format(posterior.lambda.low, digits=4), " to ", format(posterior.lambda.high, digits=4), ")", sep=""), fill = TRUE );
        posterior.formatteddays <- paste(round(posterior.estdays), " (", posterior.lowdays, ", ", posterior.uppdays, ")", sep="") 
        cat( paste("Bayesian PFitter Estimated Days:", posterior.formatteddays, sep=" "), fill = TRUE );
    }
    bayespfitter.results <-
        list(
            lambda = c(
                "estimate" = posterior.lambda,
                "2.5%" = posterior.lambda.low,
                "97.5%" = posterior.lambda.high ), 
            days = c(
                "estimate" = posterior.estdays,
                "2.5%" = posterior.lowlim,
                "97.5%" = posterior.upplim )
           );

    return( bayespfitter.results );
} # BayesPFitter (..)

makeDList <- function ( consensus.distances, intersequence.distances ) {
    dlist <- matrix( NA, nrow = length( consensus.distances ) + length( intersequence.distances ), ncol = 3 );
    colnames( dlist ) <- c( "from.seq", "to.seq", "dist" );
    dlist[ 1:length( consensus.distances ), 1 ] <- "Consensus";
    dlist[ 1:length( consensus.distances ), 2 ] <- names( consensus.distances );
    dlist[ 1:length( consensus.distances ), 3 ] <- consensus.distances;
    fromnames <- gsub( "^(.+) to .+$", "\\1", names( intersequence.distances ) );
    tonames <- gsub( "^.+ to (.+)$", "\\1", names( intersequence.distances ) );
    dlist[ length( consensus.distances ) + 1:length( intersequence.distances ), 1 ] <- fromnames;
    dlist[ length( consensus.distances ) + 1:length( intersequence.distances ), 2 ] <- tonames;
    dlist[ length( consensus.distances ) + 1:length( intersequence.distances ), 3 ] <- intersequence.distances;
    return( as.data.frame( dlist ) );
} # makeDList (..)

# Calculate and return the approximate one-sided (upper) p-value using the given null statistics (usually computed through permutation).
calculateUpperSidedPValue <- function ( the.stat, null.stats, add.one.for.observed.test.statistic = TRUE, use.ci = TRUE ) {
    if( is.na( the.stat ) ) {
        return( NA );
    } else {
        .nulls.more.extreme.than.the.stat <-
            sum( null.stats >= the.stat, na.rm=T )
            + as.numeric( add.one.for.observed.test.statistic );
        .total.nulls <-
            sum( !is.na( null.stats ) ) +
                as.numeric( add.one.for.observed.test.statistic );
        if( use.ci ) {
            if( !require( "binom" ) ) {
                install.packages( "binom", dependencies = TRUE );
                if( !require( "binom" ) ) {
                    stop( "ERROR LOADING \"binom\" package" );
                }
            }
            return( binom.confint( .nulls.more.extreme.than.the.stat, .total.nulls, conf.level = 0.95, methods = "wilson" )$upper );
        } else {
            return( .nulls.more.extreme.than.the.stat / .total.nulls );
        }
    }
} # calculateUpperSidedPValue (..)

DSPFitter <- function (
  infile = args[1],
  dlist = read.table( file=infile, sep="\t", stringsAsFactors=F ),
  epsilon = c(as.numeric(args[2])),
  nbases = c(as.numeric(args[3])),
  be.verbose = TRUE,
  MAXIMUM.DS.SAMPLE.SIZE = 250
) {
    ### Consensus distances.
    ## Sort first for efficiency.  Then expand.
    .con.mat <- dlist[ which(dlist[,1]==dlist[1,1]), , drop = FALSE ];
    .con.mat <- .con.mat[ order( .con.mat[ , 3 ] ), , drop = FALSE ]; # Sort it by column 3, distance.
    rownames( .con.mat ) <- .con.mat[ , 2 ];
    sorted.consensus.distances <-
        as.numeric( replicateDistancesForSequenceMultiplicity( .con.mat )[ , 3 ] );
    
    .num.consensus.distances <- length( sorted.consensus.distances );
    if( .num.consensus.distances > MAXIMUM.DS.SAMPLE.SIZE ) {
        new.sorted.consensus.distances.table <- round( MAXIMUM.DS.SAMPLE.SIZE * table( sorted.consensus.distances ) / .num.consensus.distances );
        sorted.consensus.distances <-
            unlist( sapply( 1:length( new.sorted.consensus.distances.table ), function( .i ) { rep( .i - 1, new.sorted.consensus.distances.table[ .i ]  ) } ) );
    }
    consensus.distances <-
        sorted.consensus.distances;
    
    ### Intersequence distances.
    ## Sort first for efficiency.  Then expand.
    .mat <- dlist[ -which(dlist[,1]==dlist[1,1]), , drop = FALSE ];
    .mat <- .mat[ order( .mat[ , 3 ] ), , drop = FALSE ]; # Sort it by column 3, distance.
    rownames( .mat ) <- apply( .mat[ , 1:2 ], 1, paste, collapse = " to " );
    sorted.intersequence.distances <-
        as.numeric( replicateDistancesForSequenceMultiplicity( .mat, missing.seqnames = .con.mat[ , 2 ] )[ , 3 ] );

    # Do this before reducing the sample size..
    DS.lambda <- mean( sorted.intersequence.distances, na.rm = TRUE );
    DS.estdays <- days( DS.lambda, nbases, epsilon );

    .num.intersequence.distances <- length( sorted.intersequence.distances );
    if( .num.intersequence.distances > MAXIMUM.DS.SAMPLE.SIZE ) {
        new.sorted.intersequence.distances.table <- round( MAXIMUM.DS.SAMPLE.SIZE * table( sorted.intersequence.distances ) / .num.intersequence.distances )
        sorted.intersequence.distances <-
            unlist( sapply( 1:length( new.sorted.intersequence.distances.table ), function( .i ) { rep( .i - 1, new.sorted.intersequence.distances.table[ .i ]  ) } ) );
    }
    intersequence.distances <-
        sorted.intersequence.distances;

    #print( DS.lambda );
    if( DS.lambda == 0 ) {
        # Special case. All zeros.
        DS.lambda.low <- 0;
        # TODO: ENSURE THIS IS THE PROPER UPPER BOUND.
        DS.lambda.high <- qgamma( .975, 1 + sum( intersequence.distances ), length( intersequence.distances ) );
        DS.upplim <- days(DS.lambda.high, nbases, epsilon)
        DS.lowlim <- days(DS.lambda.low, nbases, epsilon)
          
        DS.lowdays <- round(DS.lowlim)
        DS.uppdays <- round(DS.upplim)
    } else { # if( DS.lambda == 0 ) .. else ..
      .f <- PoissonDSM.calculatePlausibilityOfSingletons.integratedPlausibility( intersequence.distances, log.result = T );
      # Find range for integration.
      .maybe.maximal.point <- DS.lambda;
      .maximal.area.so.far <- .f( .maybe.maximal.point );
      .maximal.x.so.far <- .maybe.maximal.point;
      for( .additional.x in 1:(2*DS.lambda) ) {
          .logtotalplaus <-
              .f( .maybe.maximal.point + .additional.x );
          if( .logtotalplaus > .maximal.area.so.far ) {
              .maximal.area.so.far <- .logtotalplaus;
              .maximal.x.so.far <- .maybe.maximal.point + .additional.x;
          } else { # Stop if it's stable; Also stop if the maximal area starts to decline; that indicates we are in the range where the integral becomes unstable.
              break;
          }
      } # look up til we find a max
      .lowest.nonzero.x <- 0;
      .f.at.lowest.nonzero.x <- .f( .lowest.nonzero.x );
      while( !is.finite( .f.at.lowest.nonzero.x ) ) {
          .lowest.nonzero.x <- .lowest.nonzero.x + ( .maximal.x.so.far / 100 );
          stopifnot( .lowest.nonzero.x < DS.lambda );
          .f.at.lowest.nonzero.x <- .f( .lowest.nonzero.x );
      }
      HD.normalizedPlausibilities <- 
          Vectorize( function( y ) { exp( .f( y ) - .maximal.area.so.far ) } );
      ##  NOTE HACK using upper = HD.outrageouslyhighvalue/2
      HD.quantile <- function( quantile ) { optimize( function( y ) { abs( HD.normalizedPlausibilities( y ) - quantile ) }, lower = .lowest.nonzero.x, upper = .maximal.x.so.far )$minimum }
      
      DS.lambda.low <- HD.quantile( 0.025 );
      if( DS.lambda.low > DS.lambda ) {
          if( be.verbose ) {
              cat( "WARNING: NUMERICAL ISSUES PREVENT GETTING AN ACCURATE DS POSTERIOR INTERVAL.", fill = TRUE );
          }
          DS.lambda.low <- 0;
          DS.lambda.high <- Inf;
          DS.lowlim <- 0;
          DS.lowdays <- 0;
          DS.upplim <- Inf;
          DS.uppdays <- Inf;
      } else {
          DS.lambda.high <- HD.quantile( 0.975 );
          DS.upplim <- days(DS.lambda.high, nbases, epsilon)
          DS.lowlim <- days(DS.lambda.low, nbases, epsilon)
          
          DS.lowdays <- round(DS.lowlim)
          DS.uppdays <- round(DS.upplim)
      }
    } # End if( DS.lambda == 0 ) .. else ..
    if( be.verbose ) {
        cat( paste("DS PFitter Estimated Lambda is ", format(DS.lambda, digits=4), " (95% CI ", format(DS.lambda.low, digits=4), " to ", format(DS.lambda.high, digits=4), ")", sep=""), fill = T );
        DS.formatteddays <- paste(round(DS.estdays), " (", DS.lowdays, ", ", DS.uppdays, ")", sep="");
        cat( paste("DS PFitter Estimated Days:", DS.formatteddays, sep=" "), fill = T );
    }
    dspfitter.results <- list( lambda = c( "estimate" = DS.lambda, "2.5%" = DS.lambda.low, "97.5%" = DS.lambda.high ), days = c( "estimate" = DS.estdays, "2.5%" = DS.lowlim, "97.5%" = DS.upplim ) );
    
    ## For a given set of data and corresponding set of pepr variates, compute the smallest epsilon and corresponding lambda such that the data CDF passes an average of epsilon units (vertically) from Poisson( lambda ), averaged over the observed data points.
    PoissonDSM.getSmallestEpsilonAndLambda <- function ( observed.data, pepr.variates = runif( length( observed.data ) ), observed.data.are.already.in.ascending.order = FALSE, pepr.variates.are.already.in.ascending.order = FALSE, pepr.variates.are.validity.checked = FALSE, compare.to.specific.pois.with.lambda = NULL ) {
        stopifnot( length( observed.data ) > 0 );
        stopifnot( length( pepr.variates ) == length( observed.data ) );
        if( !observed.data.are.already.in.ascending.order ) {
            observed.data <- sort( observed.data );
            observed.data.are.already.in.ascending.order <- TRUE;
        }
        if( !pepr.variates.are.already.in.ascending.order ) {
            pepr.variates <- sort( pepr.variates );
            pepr.variates.are.already.in.ascending.order <- TRUE;
        }
        if( !pepr.variates.are.validity.checked ) {
            stopifnot( ( pepr.variates[ 1 ] >= 0 ) && ( pepr.variates[ 1 ] <= pepr.variates[ length( pepr.variates ) ] ) && ( pepr.variates[ length( pepr.variates ) ] <= 1 ) );
            pepr.variates.are.validity.checked <- TRUE;
        }
        if( !is.null( compare.to.specific.pois.with.lambda ) ) {
            # Force lambda for comparison.
            lambda.low <- compare.to.specific.pois.with.lambda;
            lambda.high <- compare.to.specific.pois.with.lambda;
        } else if( length( observed.data ) == 1 ) {
            lambda.low <- observed.data;
            lambda.high <- observed.data;
        } else {
            lambda.low <- max( min( observed.data ), ( mean( observed.data ) - 4 * sd( observed.data ) ) );
            lambda.high <- min( max( observed.data ), ( mean( observed.data ) + 4 * sd( observed.data ) ) );
        }
        .epsilon.from.lambda.fn <- 
          function( lambda ) {
            # Here we use the assumption that the observed data are in ascending order, which implies that the quantiles are, too.
            .quantiles <- ppois( observed.data, lambda );
            .abs.diffs <- abs( pepr.variates - .quantiles );
            # DON'T Return the largest difference.
            #return( max( .abs.diffs ) );
            # Return the average difference.
            return( mean( .abs.diffs ) );
          };
        if( lambda.low == lambda.high ) {
            the.results <- list( "lambda" = lambda.low,
                             "epsilon" = .epsilon.from.lambda.fn( lambda.low ) );
        } else {
          the.results <- 
          optimize( f = .epsilon.from.lambda.fn, lower = lambda.low, upper = lambda.high, maximum = FALSE, tol = 1E-5 );
          names( the.results ) <- c( "lambda", "epsilon" );
        }
        return( the.results );
    } # PoissonDSM.getSmallestEpsilonAndLambda ( observed.data, pepr.variates )
    .sampled.epsilons.and.lambdas <- sapply( 1:EPSILON.NDRAWS, function( .x ) {
        #print( .x ); 
        return( PoissonDSM.getSmallestEpsilonAndLambda( sorted.intersequence.distances, observed.data.are.already.in.ascending.order = TRUE ) ); } );
    ds.posterior.epsilons <-
        as.numeric( .sampled.epsilons.and.lambdas[ "epsilon", ] );
    
    if( be.verbose ) {
        cat( paste( "It seems that the CDF of the closest Poisson distribution is roughly ", sprintf( "%2.1f", ( 100 * median( ds.posterior.epsilons ) ) ), "% away from the pepr-sampled empirical CDFs (middle 95% ", sprintf( "%2.1f", ( 100 * quantile( ds.posterior.epsilons, .025 ) ) ), " to ", sprintf( "%2.1f", ( 100 * quantile( ds.posterior.epsilons, .975 ) ) ), ").", sep = "" ), fill = TRUE );
    }
    
    dspfitter.results <- c( dspfitter.results, list( epsilon = list( "estimate" = median( ds.posterior.epsilons ), "2.5%" = quantile( ds.posterior.epsilons, .025 ), "97.5%" = quantile( ds.posterior.epsilons, .975 ) ) ) );

    create.intersequence.distances.by.convolving.consensus.distances <- function ( table.of.consensus.distances ) {
        nl0 <- 1 + max( as.numeric( names( table.of.consensus.distances ) ) );
        yvec0 <- rep( 0, nl0 );
        yvec0[ 1 + as.numeric( names( table.of.consensus.distances ) ) ] <- table.of.consensus.distances;
        #print( yvec0 );
        xvec1 <- c(yvec0, rep(0,nl0))
        
        yvec1 <- rep(0,2*nl0)
        for( m in 1:length( yvec1 ) ) {
            yvec1[ m ] <- ( 1 / 2 ) * sum( unlist( sapply( 1:m, function( k ) {
                    xvec1[ k ] * xvec1[ m - k + 1 ]
                } ) ) );
            if( m %% 2 == 1 ) {
                # Odd m.
                yvec1[ m ] <- yvec1[ m ] - ( 1 / 2 ) * xvec1[ (((m-1)/2)+1) ];
            }
        }
        # print( "yvec1" );
        # print( yvec1 );
        # 
        # yvec1 <- rep(0,2*nl0)
        # yvec1[1] <- 1/2*yvec0[1]*(yvec0[1]-1)  ### freq at zero 
        # for( i in 1:nl0 ) {
        #     m <- i * 2;
        #     for(k in 1:m){
        #         .xvec1.k <- xvec1[k];
        #         .xvec1.m.minus.k.plus.1 <- xvec1[m-k+1];
        #         #yvec1[m] <- yvec1[m] + 1/2*xvec1[k]*xvec1[m-k+1]
        #         yvec1[m] <- yvec1[m] + 1/2*.xvec1.k*.xvec1.m.minus.k.plus.1;
        #         if(i != nl0) {
        #             #yvec1[m+1] <- yvec1[m+1] + 1/2*xvec1[k]*(xvec1[m-k+2] - ifelse( k == (i+1), 1, 0 ));
        #             .xvec1.m.minus.k.plus.2 <- xvec1[m-k+2];
        #             .xvec1.m.minus.k.plus.2.maybeminus1 <- .xvec1.m.minus.k.plus.2 - ifelse( k == (i+1), 1, 0 );
        #             yvec1[m+1] <- yvec1[m+1] + 1/2*.xvec1.k*.xvec1.m.minus.k.plus.2.maybeminus1;
        #         } 
        #     }
        #     if(i != nl0) { yvec1[m+1] <- yvec1[m+1] + 1/2*xvec1[m+1]*xvec1[1] }
        # }
        # print( "yvec1" );
        # print( yvec1 );
        return( yvec1 );
    } # create.intersequence.distances.by.convolving.consensus.distances ( table.of.consensus.distances )
    
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
    
    if( TRUE ) {
      n.seqs <- length( consensus.distances );
      ## The strategy here is to try 10, if at least half are higher than the test stat, do another 90, and if of those 100 the crude p-value is <= 0.10, get a precise p-value by doing a total of 500.  We threshold for further refinement based on a per-test probabilty of failing to refine (when we should refine ie p = 0.05) of 1E-5.
      POISSON.DRAW.REPS.DABBLE <- 10;
      POISSON.DRAW.REPS.INITIAL <- 100;
      POISSON.DRAW.REPS.TOTAL <- 500;
      prior.gamma.shape <- 1;
      prior.gamma.rate <- 1;
      .lambdas <- rgamma( POISSON.DRAW.REPS.TOTAL, prior.gamma.shape + sum( consensus.distances ), prior.gamma.rate + length( consensus.distances ) );
      .ds.rep.results <- lapply( 1:POISSON.DRAW.REPS.DABBLE, function( .rep.i ) {
        ## Draw the null lambda from the posterior distribution of
        ## lambdas after seeing the consensus data.
        if( be.verbose ) {
            cat( ".rep.i: ", .rep.i, " (.lambda: ", .lambdas[ .rep.i ], ")", sep = "", fill = TRUE );
        }

        .consensus.distances <- rpois( n.seqs, .lambdas[ .rep.i ] );
        .table.of.intersequence.distances <- create.intersequence.distances.by.convolving.consensus.distances( table( .consensus.distances ) );
          
        # enforce maximum
        .num.intersequence.distances <- sum( .table.of.intersequence.distances );
        if( .num.intersequence.distances > MAXIMUM.DS.SAMPLE.SIZE ) {
            .table.of.intersequence.distances <- round( MAXIMUM.DS.SAMPLE.SIZE * .table.of.intersequence.distances / .num.intersequence.distances );
        }
          
        .sorted.intersequence.distances <-
            unlist( sapply( 1:length( .table.of.intersequence.distances ), function( .i ) { rep( .i - 1, .table.of.intersequence.distances[ .i ]  ) } ) );
        .sampled.epsilons.and.lambdas <- sapply( 1:EPSILON.NDRAWS, function( .x ) {
            return( PoissonDSM.getSmallestEpsilonAndLambda( .sorted.intersequence.distances, observed.data.are.already.in.ascending.order = TRUE ) );
        } );
        .ds.posterior.epsilons <- as.numeric( .sampled.epsilons.and.lambdas[ "epsilon", ] );
        .quantiles.of.interest <- quantile( .ds.posterior.epsilons, probs = c( .025, .5, .975 ) );
        .dspfitter.results <- list( epsilon = list( "estimate" = .quantiles.of.interest[ 2 ], "2.5%" = .quantiles.of.interest[ 1 ], "97.5%" = .quantiles.of.interest[ 3 ] ) );
       # print( .dspfitter.results );
      
        return( .dspfitter.results );
      } );
      
      #print( ".ds.rep.results" );
      #print( .ds.rep.results );
      .epsilon.medians <- sapply( .ds.rep.results, function( .lst ) { unname( .lst[[ "epsilon" ]][[ "estimate" ]] ) } );
      .epsilon.lowers <- sapply( .ds.rep.results, function( .lst ) { unname( .lst[[ "epsilon" ]][[ "2.5%" ]] ) } );
      .epsilon.uppers <- sapply( .ds.rep.results, function( .lst ) { unname( .lst[[ "epsilon" ]][[ "97.5%" ]] ) } );
  
      .p.values <- list( "estimate.p.value" = calculateUpperSidedPValue( dspfitter.results[[ "epsilon" ]][[ "estimate" ]], .epsilon.medians ), "lower.p.value" = calculateUpperSidedPValue( dspfitter.results[[ "epsilon" ]][[ "2.5%" ]], .epsilon.lowers ), "upper.p.value" = calculateUpperSidedPValue( dspfitter.results[[ "epsilon" ]][[ "97.5%" ]], .epsilon.uppers ) );
      .thresh.dabble <-
          ( qbinom( 1-1E-5, POISSON.DRAW.REPS.DABBLE, .05 ) / POISSON.DRAW.REPS.DABBLE );
      if( ( any( .p.values <= .thresh.dabble ) ) && !( all( .p.values <= 0.05 ) ) ) {
        if( be.verbose ) {
          cat( "Running ", POISSON.DRAW.REPS.INITIAL - POISSON.DRAW.REPS.DABBLE, " more reps, because p-value after ", POISSON.DRAW.REPS.DABBLE, " reps the p-values are (at least one <= threshold ", sprintf( "%0.4f", .thresh.dabble ), "): ", sep = "", fill = TRUE );
          print( .p.values );
          print( paste( "Estimate:", dspfitter.results[[ "epsilon" ]][[ "estimate" ]] ) );
          print( "summary of null epsilon medians:" );
          print( summary( .epsilon.medians ) );
        }
        .ds.rep.results <- c( .ds.rep.results, lapply( ( POISSON.DRAW.REPS.DABBLE + 1):POISSON.DRAW.REPS.INITIAL, function( .rep.i ) {
          if( be.verbose ) {
              cat( ".rep.i: ", .rep.i, fill = TRUE );
          }
          .consensus.distances <- rpois( n.seqs, .lambdas[ .rep.i ] );
          .table.of.intersequence.distances <- create.intersequence.distances.by.convolving.consensus.distances( table( .consensus.distances ) );
          # enforce maximum
          .num.intersequence.distances <- sum( .table.of.intersequence.distances );
          if( .num.intersequence.distances > MAXIMUM.DS.SAMPLE.SIZE ) {
              .table.of.intersequence.distances <- round( MAXIMUM.DS.SAMPLE.SIZE * .table.of.intersequence.distances / .num.intersequence.distances );
          }
          
          .sorted.intersequence.distances <-
              unlist( sapply( 1:length( .table.of.intersequence.distances ), function( .i ) { rep( .i - 1, .table.of.intersequence.distances[ .i ]  ) } ) );
          .sampled.epsilons.and.lambdas <- sapply( 1:EPSILON.NDRAWS, function( .x ) {
              return( PoissonDSM.getSmallestEpsilonAndLambda( .sorted.intersequence.distances, observed.data.are.already.in.ascending.order = TRUE ) );
          } );
          .ds.posterior.epsilons <- as.numeric( .sampled.epsilons.and.lambdas[ "epsilon", ] );
          .quantiles.of.interest <- quantile( .ds.posterior.epsilons, probs = c( .025, .5, .975 ) );
          .dspfitter.results <- list( epsilon = list( "estimate" = .quantiles.of.interest[ 2 ], "2.5%" = .quantiles.of.interest[ 1 ], "97.5%" = .quantiles.of.interest[ 3 ] ) );
         # print( .dspfitter.results );
        
          return( .dspfitter.results );
        } ) );
          
        #print( ".ds.rep.results" );
        #print( .ds.rep.results );
        .epsilon.medians <- sapply( .ds.rep.results, function( .lst ) { unname( .lst[[ "epsilon" ]][[ "estimate" ]] ) } );
        .epsilon.lowers <- sapply( .ds.rep.results, function( .lst ) { unname( .lst[[ "epsilon" ]][[ "2.5%" ]] ) } );
        .epsilon.uppers <- sapply( .ds.rep.results, function( .lst ) { unname( .lst[[ "epsilon" ]][[ "97.5%" ]] ) } );
        .p.values <- list( "estimate.p.value" = calculateUpperSidedPValue( dspfitter.results[[ "epsilon" ]][[ "estimate" ]], .epsilon.medians ), "lower.p.value" = calculateUpperSidedPValue( dspfitter.results[[ "epsilon" ]][[ "2.5%" ]], .epsilon.lowers ), "upper.p.value" = calculateUpperSidedPValue( dspfitter.results[[ "epsilon" ]][[ "97.5%" ]], .epsilon.uppers ) );
        .thresh <-
            qbinom( 1-1E-5, POISSON.DRAW.REPS.INITIAL, .05 ) / POISSON.DRAW.REPS.INITIAL;
        if( any( .p.values <= .thresh ) && !( all( .p.values <= 0.05 ) ) ) {
          if( be.verbose ) {
            cat( "Running ", POISSON.DRAW.REPS.TOTAL - POISSON.DRAW.REPS.INITIAL, " more reps, because p-value after ", POISSON.DRAW.REPS.INITIAL, " reps the p-values are (at least one <= threshold ", sprintf( "%0.4f", .thresh ), "): ", sep = "", fill = TRUE );
            print( .p.values );
          }
          .ds.rep.results <- c( .ds.rep.results, lapply( ( POISSON.DRAW.REPS.INITIAL + 1):POISSON.DRAW.REPS.TOTAL, function( .rep.i ) {
            if( be.verbose ) {
                cat( ".rep.i: ", .rep.i, fill = TRUE );
            }
            .consensus.distances <- rpois( n.seqs, .lambdas[ .rep.i ] );
            .table.of.intersequence.distances <- create.intersequence.distances.by.convolving.consensus.distances( table( .consensus.distances ) );
            
            # enforce maximum
            .num.intersequence.distances <- sum( .table.of.intersequence.distances );
            if( .num.intersequence.distances > MAXIMUM.DS.SAMPLE.SIZE ) {
                .table.of.intersequence.distances <- round( MAXIMUM.DS.SAMPLE.SIZE * .table.of.intersequence.distances / .num.intersequence.distances );
            }
          
            .sorted.intersequence.distances <-
                unlist( sapply( 1:length( .table.of.intersequence.distances ), function( .i ) { rep( .i - 1, .table.of.intersequence.distances[ .i ]  ) } ) );
            .sampled.epsilons.and.lambdas <- sapply( 1:EPSILON.NDRAWS, function( .x ) {
                return( PoissonDSM.getSmallestEpsilonAndLambda( .sorted.intersequence.distances ) );
            } );
            .ds.posterior.epsilons <- as.numeric( .sampled.epsilons.and.lambdas[ "epsilon", ] );
            .quantiles.of.interest <- quantile( .ds.posterior.epsilons, probs = c( .025, .5, .975 ) );
            .dspfitter.results <- list( epsilon = list( "estimate" = .quantiles.of.interest[ 2 ], "2.5%" = .quantiles.of.interest[ 1 ], "97.5%" = .quantiles.of.interest[ 3 ] ) );
           # print( .dspfitter.results );
          
            return( .dspfitter.results );
          } ) );
            
          #print( ".ds.rep.results" );
          #print( .ds.rep.results );
          .epsilon.medians <- sapply( .ds.rep.results, function( .lst ) { unname( .lst[[ "epsilon" ]][[ "estimate" ]] ) } );
          .epsilon.lowers <- sapply( .ds.rep.results, function( .lst ) { unname( .lst[[ "epsilon" ]][[ "2.5%" ]] ) } );
          .epsilon.uppers <- sapply( .ds.rep.results, function( .lst ) { unname( .lst[[ "epsilon" ]][[ "97.5%" ]] ) } );
          .p.values <- list( "estimate.p.value" = calculateUpperSidedPValue( dspfitter.results[[ "epsilon" ]][[ "estimate" ]], .epsilon.medians ), "lower.p.value" = calculateUpperSidedPValue( dspfitter.results[[ "epsilon" ]][[ "2.5%" ]], .epsilon.lowers ), "upper.p.value" = calculateUpperSidedPValue( dspfitter.results[[ "epsilon" ]][[ "97.5%" ]], .epsilon.uppers ) );
        } # End if we should refine the p values further.
      } # End if we should bother doing reps for the p values 
      if( be.verbose ) {
        if( .p.values[[ "estimate.p.value" ]] <= 0.05 ) {
            cat( "DSPFitter convolution test: DOES NOT FOLLOW A STAR-PHYLOGENY (P", prettyPrintPValuesTo4Digits( .p.values[[ "estimate.p.value" ]], equalsText = " = ", lessThanText = " < " ), ")", sep = "", fill = TRUE );
        } else {
            cat( "DSPFitter convolution test: FOLLOWS A STAR-PHYLOGENY (P", prettyPrintPValuesTo4Digits( .p.values[[ "estimate.p.value" ]], equalsText = " = ", lessThanText = " < " ), ")", sep = "", fill = TRUE );
        }
          ## TODO: REMOVE
        if( .p.values[[ "lower.p.value" ]] <= 0.05 ) {
            cat( "DSPFitter lower convolution test: DOES NOT FOLLOW A STAR-PHYLOGENY (P", prettyPrintPValuesTo4Digits( .p.values[[ "lower.p.value" ]], equalsText = " = ", lessThanText = " < " ), ")", sep = "", fill = TRUE );
        } else {
            cat( "DSPFitter lower convolution test: FOLLOWS A STAR-PHYLOGENY (P", prettyPrintPValuesTo4Digits( .p.values[[ "lower.p.value" ]], equalsText = " = ", lessThanText = " < " ), ")", sep = "", fill = TRUE );
        }
          ## TODO: REMOVE
        if( .p.values[[ "upper.p.value" ]] <= 0.05 ) {
            cat( "DSPFitter upper convolution test: DOES NOT FOLLOW A STAR-PHYLOGENY (P", prettyPrintPValuesTo4Digits( .p.values[[ "upper.p.value" ]], equalsText = " = ", lessThanText = " < " ), ")", sep = "", fill = TRUE );
        } else {
            cat( "DSPFitter upper convolution test: FOLLOWS A STAR-PHYLOGENY (P", prettyPrintPValuesTo4Digits( .p.values[[ "upper.p.value" ]], equalsText = " = ", lessThanText = " < " ), ")", sep = "", fill = TRUE );
        }
      } # End if be.verbose
      dspfitter.results[[ "epsilon" ]] <- c( dspfitter.results[[ "epsilon" ]], .p.values );
    } # END IF we should do sampling to p-values for the epsilon distances.
    
    draw.one.nonconflicted.pepr.sample <- function ( sorted.observed.data, conflict.max = 5000, method = "fast" ) {
        ## The naive way, just keep drawing until one works.
          found.it <- FALSE;
          conflict <- 0;
          while( !found.it && conflict < conflict.max ) {
            # Directly draw from the Poisson DSM
            lambda.low <-
                rgamma( n = length( sorted.observed.data ), shape = sorted.observed.data );
              lambda.width <-
                  rgamma( n = length( sorted.observed.data ), shape = 1 );
              if( any( max( lambda.low ) > ( lambda.low + lambda.width ) ) ) {
                  if( method == "naive" ) {
                    ## CONFLICT.
                    conflict <- conflict + 1;
                    cat( "." );
                    next;
                  } else {
                    # All of them need to overlap, so every value of lambda.width must be at least as long as max( lambda.low ) - lambda.low.  They are iid, so the condition that any one of them is at least that value has no effect on the rest of them.
                    # But also, the distribution of an exponential is memoryless, so we can just shift everything. Note we do need to draw them again so they are not truncated (because we've conditioned on them being a minimum, but they're memoryless..)
                      lambda.width <-
                          ( max( lambda.low ) - lambda.low ) +
                          rgamma( n = length( sorted.observed.data ), shape = 1 );
                  } # End if method == "naive" .. else ..
              } # End if the draws aren't long enough.
              stopifnot( !any( max( lambda.low ) > ( lambda.low + lambda.width ) ) );
            found.it <- 1;
          } # end while( !found.it )
          if( !found.it ) { stop( "TOO MUCH CONFLICT!" ) }
          
        return( list( lambda.low = lambda.low, lambda.width = lambda.width ) );
    } # draw.one.nonconflicted.pepr.sample 
    
    draw.once.from.each.and.evaluate.assertion <- function( assertion.is.that.1.is.x.times.2.low = 1.5, assertion.is.that.1.is.x.times.2.high = 2.5, sorted.observed.data.1 = sorted.consensus.distances, sorted.observed.data.2 = sorted.intersequence.distances, be.verbose = FALSE, weakening.type = NA ) {
        stopifnot( length( assertion.is.that.1.is.x.times.2.low ) == length( assertion.is.that.1.is.x.times.2.high ) );
        stopifnot( all( assertion.is.that.1.is.x.times.2.low <= assertion.is.that.1.is.x.times.2.high ) );
        .sample1 <- draw.one.nonconflicted.pepr.sample( sorted.observed.data.1 );
        .all.width.rotations <- sapply( 1:length( .sample1$lambda.width ), function( .start.i ) { c( .sample1$lambda.width, .sample1$lambda.width )[ .start.i + 0:( length( .sample1$lambda.width ) - 1 ) ] } );
        if( is.na( weakening.type ) ) {
            .weakening.dependent.part <- min( .sample1$lambda.low + .sample1$lambda.width );
        } else if( weakening.type == "rotate" ) {
            .weakening.dependent.part <- max( apply( .all.width.rotations, 1, function( .widths ) { min( .sample1$lambda.low + .widths ) } ) );
        } else if( weakening.type == "complete" ) {
            .weakening.dependent.part <- min( sort( .sample1$lambda.low ) + sort( .sample1$lambda.width, decreasing = TRUE ) );
        } else {
            stop( paste( "unknown weakening type:", weakening.type ) );
        }
        .sample1.max = lapply( 1:length( assertion.is.that.1.is.x.times.2.low ), function( .i ) {
            assertion.is.that.1.is.x.times.2.high[ .i ] * .weakening.dependent.part;
        } );
        .weakening.dependent.part <- max( .sample1$lambda.low ); # Not actually dependent on how we weaken.
        .sample1.min = lapply( 1:length( assertion.is.that.1.is.x.times.2.low ), function( .i ) {
            assertion.is.that.1.is.x.times.2.low[ .i ] * .weakening.dependent.part;
        } );
        stopifnot( all( unlist( lapply( 1:length( assertion.is.that.1.is.x.times.2.low ), function( .i ) { .sample1.min[[ .i ]] <= .sample1.max[[ .i ]] } ) ) ) );
        if( be.verbose ) {
            cat( paste( sapply( 1:length( .sample1.min ), function( .i ) { paste( c( ".sample1.min", .sample1.min[[ .i ]],  ".sample1.max", .sample1.max[[ .i ]] ), collapse = " " ) } ), collapse = "\n" ), fill = T );
        }
        .sample2 <- draw.one.nonconflicted.pepr.sample( sorted.observed.data.2 );
        .sample2.max = min( .sample2$lambda.low + .sample2$lambda.width );
        .sample2.min = max( .sample2$lambda.low );
        stopifnot( .sample2.min <= .sample2.max );
        if( be.verbose ) {
            cat( paste( c( ".sample2.min", .sample2.min,  ".sample2.max", .sample2.max ), collapse = " " ), fill = T );
        }
        .rv <-
            sapply( 1:length( .sample1.min ), function( .i ) {
        if( ( .sample1.max[[ .i ]] < .sample2.min ) ||
                ( .sample1.min[[ .i ]] > .sample2.max ) ) {
            # Evidence against.
            if( be.verbose ) {
                cat( "Q", fill = T );
            }
            return( list( P = 0, Q = 1, R = 0 ) );
        }
        if( ( .sample1.min[[ .i ]] >= .sample2.min ) && ( .sample1.max[[ .i ]] <= .sample2.max ) ) {
            # Evidence in support.
            if( be.verbose ) {
                cat( "P", fill = T );
            }
            return( list( P = 1, Q = 0, R = 0 ) );
        }
        # If it's not evidence for or against, then we are in a "don't know" situation.
        if( be.verbose ) {
            cat( "R", fill = T );
        }
        return( list( P = 0, Q = 0, R = 1 ) );
            } );
        if( is.null( dim( .rv ) ) ) {
            .rv <- matrix( .rv, nrow = 3 );
        }
        rownames( .rv ) <- c( "P", "Q", "R" );
        return( .rv );
    } # draw.once.from.each.and.evaluate.assertion (..)
    

    # GET A BUNCH OF DRAWS
    .result <- sapply( 1:DS.NDRAWS, function( .i ) { draw.once.from.each.and.evaluate.assertion() } );
    rownames( .result ) <- c( "P", "Q", "R" );
    DS.P <- mean( unlist( .result[ "P", ] ) );
    DS.Q <- mean( unlist( .result[ "Q", ] ) );
    DS.R <- mean( unlist( .result[ "R", ] ) );
    if( DS.P == 0 ) {
        DS.Ptext <- paste( "<=", ( 1 / DS.NDRAWS ) );
    } else {
        DS.Ptext <- sprintf( paste( "= %0.", ceiling( log10( DS.NDRAWS ) ), "f", sep = "" ), DS.P )
    }
    if( DS.Q == 0 ) {
        DS.Qtext <- paste( "<=", ( 1 / DS.NDRAWS ) );
    } else {
        DS.Qtext <- sprintf( paste( "= %0.", ceiling( log10( DS.NDRAWS ) ), "f", sep = "" ), DS.Q )
    }
    if( DS.R == 0 ) {
        DS.Rtext <- paste( "<=", ( 1 / DS.NDRAWS ) );
    } else {
        DS.Rtext <- sprintf( paste( "= %0.", ceiling( log10( DS.NDRAWS ) ), "f", sep = "" ), DS.R )
    }
    #print( paste( "The DS evidence against the assertion that the Poisson rate between sequences is twice the rate of sequences to the consensus is Q ", DS.Qtext, ". The remaining evidence (", sprintf( paste( "%0.", ceiling( log10( DS.NDRAWS ) ), "f", sep = "" ), DS.R ), ") neither supports nor contradicts the assertion.", sep = "" ) );
    if( be.verbose ) {
      if( DS.Q >= 0.95 ) {
          cat( "DSPFitter test that intersequence rate = 2 x seq-consensus rate: BAD", fill = TRUE );
          cat( paste( "\tThere is evidence against the assertion that the Poisson rate between sequences is between 1.5 and 2.5 times the rate of sequences to the consensus (R ", DS.Rtext, ").", sep = "" ), fill = TRUE );
      } else {
          cat( "DSPFitter test that intersequence rate = 2 x seq-consensus rate: OK", fill = TRUE );
          cat( paste( "\tThere is not sufficient evidence against the assertion that the Poisson rate between sequences is between 1.5 and 2.5 times the rate of sequences to the consensus (R ", DS.Rtext, ").", sep = "" ), fill = TRUE );
      }
    }
    dspfitter.results <- c( dspfitter.results, list( twox = list( "R" = DS.R ) ) );

    ## TODO: REMOVE. TESTING.
    # POISSON.DRAW.REPS <- 100;#1000;
    # #DS.DECISION.Q.THRESH <- 0.95;
    # DS.DECISION.Q.THRESH <- NA;
    # #n.seqs <- c( 10, 100, 500 );
    # n.seqs <- length( consensus.distances );
    # weakening.type <- NA;
    # #weakening.type <- "complete";
    # results.by.n.seqs <- 
    # lapply( n.seqs, function( .n.seqs ) {
    #     cat( "\n.n.seqs: ", .n.seqs, fill = TRUE );
    #     .diffs <- 2;#c( 1.5, 2, 2.5 );
    #     .rv <- 
    #       lapply( .diffs, function( .diff ) {
    #           cat( ".diff: ", .diff, fill = TRUE );
    #           .assertion.widths <- 1;#c( seq( 0, .05, by = .01 ), seq( .10, .25, by = .05 ), seq( .5, 2, by = .25 ) );
    #           .ds.rep.results <- lapply( 1:POISSON.DRAW.REPS, function( .rep.i ) {
    #             cat( ".rep.i: ", .rep.i, fill = TRUE );
    #             .consensus.distances <- rpois( .n.seqs, 1 );
    #             .table.of.intersequence.distances <- create.intersequence.distances.by.convolving.consensus.distances( table( .consensus.distances ) );
    #             .sorted.intersequence.distances <- unlist( sapply( 1:length( .table.of.intersequence.distances ), function( .i ) { rep( .i - 1, .table.of.intersequence.distances[ .i ]  ) } ) );
    #             #.sorted.intersequence.distances <- sort( rpois( choose( .n.seqs, 2 ), .diff ) );
    #            .rv.list <- list();
    #            .sorted.consensus.distances <- sort( .consensus.distances );
    #            .result <- lapply( 1:DS.NDRAWS, function( .i ) { cat( ".", fill = F ); draw.once.from.each.and.evaluate.assertion( 2 - .assertion.widths, 2 + .assertion.widths, .sorted.consensus.distances, .sorted.intersequence.distances, weakening.type = weakening.type ) } );
    #            DS.P <- lapply( 1:length( .assertion.widths ), function( .assertion.width.i ) { mean( unlist( sapply( .result, function( .mat ) { .mat[ "P", .assertion.width.i ] } ) ) ) } );
    #            names( DS.P ) <- .assertion.widths;
    #            DS.Q <- lapply( 1:length( .assertion.widths ), function( .assertion.width.i ) { mean( unlist( sapply( .result, function( .mat ) { .mat[ "Q", .assertion.width.i ] } ) ) ) } );
    #            names( DS.Q ) <- .assertion.widths;
    #            DS.R <- lapply( 1:length( .assertion.widths ), function( .assertion.width.i ) { mean( unlist( sapply( .result, function( .mat ) { .mat[ "R", .assertion.width.i ] } ) ) ) } );
    #            names( DS.R ) <- .assertion.widths;
    #             ## ERE I AM.  Just changed so calc an allow DS.P > 0 (since assertion has width).  Next things are to a) run PFitter(..) with these same data to evaluate decision properties, just with "estimate.lambda.and.days = FALSE, do.chisq.gof = TRUE". b) try weakening to see if that addresses the overconfidence against exact assertions.  c) draw from a "convolution" distribution to simulate the non-independence.
    #             ## TODO: REMOVE.  ERE I AM. TESTING.
    #             if( any( DS.P > 0 ) ) {
    #                 print( "DS.P:" );
    #                 print( DS.P );
    #             }
    #             #stopifnot( all( unlist( DS.Q ) + unlist( DS.R ) == 1 ) );
    #             stopifnot( all( unlist( DS.Q ) + unlist( DS.P ) + unlist( DS.R ) == 1 ) );
    #             .rv.list <- list( "assertion.PQR" = list( "P" = DS.P, "Q" = DS.Q, "R" = DS.R ) );
    #             return( .rv.list );
    #           } );
    #           .rv.list.for.diff <- list();
    #           .ds.qs <- sapply( .ds.rep.results, function( .lst ) { .lst[[ "assertion.PQR" ]][[ "Q" ]] } );
    #           .ds.ps <- sapply( .ds.rep.results, function( .lst ) { .lst[[ "assertion.PQR" ]][[ "P" ]] } );
    #           .ds.rs <- sapply( .ds.rep.results, function( .lst ) { .lst[[ "assertion.PQR" ]][[ "R" ]] } );
    #           .rv.list.for.diff <- list( "DS.Ps" = .ds.ps, "DS.Qs" = .ds.qs, "DS.Rs" = .ds.rs );
    #           ## Dichotomize based on decision threshold Q > DS.DECISION.Q.THRESH.
    #           ## Return the fraction of reps in which the decision is made.
    #           if( !is.na( DS.DECISION.Q.THRESH ) ) {
    #               .decisions <- ( .ds.qs > DS.DECISION.Q.THRESH );
    #               .rv.list.for.diff <- c( .rv.list.for.diff, decision.fraction = apply( .decisions, 1, mean ) );
    #           }
    #           return( .rv.list.for.diff );
    #         } );
    #       names( .rv ) <- .diffs;
    #       return( .rv );
    # } );
    # names( results.by.n.seqs ) <- n.seqs;

    return( dspfitter.results );
} # DSPFitter (..)



###################################################
### code chunk number 6: DSPFitter.Rnw:1282-1285
###################################################
#.result.ignored <- PFitter( be.verbose = TRUE );
.result.ignored <- BayesPFitter( be.verbose = TRUE );
.result.ignored <- DSPFitter( be.verbose = TRUE );


###################################################
### code chunk number 7: DSPFitter.Rnw:1292-1294
###################################################
# (un)Setup for prettier Sweave output.
options( continue = old.continue.option$continue )



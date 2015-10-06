# Paul Edlefsen downloaded this on Oct 4, 2015 from "http://www.hiv.lanl.gov/tmp/POISSON_FITTER/pp429jC8Td/PFitter.R".  (followed link 'Download R Script' from PoissonFitter results page http://www.hiv.lanl.gov/cgi-bin/POISSON_FITTER/v2/pfitter.cgi)
#########################################################################
##########       PURPOSE: fit Poisson, calculate lambdas,      ##########
##########    U-stat standard deviation, and goodness of fit   ##########
##########      written by EEG, last modified on 5/26/15       ##########
##########          send questions to egiorgi@lanl.gov         ##########

# INPUT: 
#  pairwise hamming distances file, mutation rate, length of a sequence
#  distances file: tab-delimited 3-column file. seqname1(1st col), 
#                  seqname2(2nd) and distance between seq1 and seq2(3rd).
#                  based on large-scale formatted sequnce input, which 
#                  means every unique sequence is represented only once 
#                  and a seqname should end with _nnn wher nnn is the 
#                  the multiplicity of such sequence in the alignment.

#  example:  R CMD BATCH '--vanilla --args sample.dist  2.16e-05 2517' this_script 

# OUTPUT:
#  2 files, one with lambdas, maxhd, jackknife standard deviation 
#  and estimated days, goodness of fit p-values, and the other with the 
#  star-phylogeny estimated an dobserved numbers (if they coincide, you have a star)
#########################################################################


library(tools)
args = commandArgs(trailingOnly=TRUE)

infile <- args[1]
epsilon <- c(as.numeric(args[2]))
nbases <- c(as.numeric(args[3]))

sample <- basename(file_path_sans_ext(infile))
dir <- paste(dirname(infile), '/', sep='')
outfile <- paste(dir, "LOG_LIKELIHOOD.results.txt", sep="")
outfile2 <- paste(dir, "CONVOLUTION.results.txt", sep="")

### FUNCTIONS ###
iseven <- function(c) {
	c1 <- c/2-as.integer(c/2)
	if(c1==0){ ev <- TRUE } else { ev <- FALSE }
	return(ev)
}
		
phi <- sqrt(1+4/3)
#gens <- function(l,nb,epsilon) (l/(nb*epsilon)-(1-phi)/(phi^2))*((phi)/(1+phi))
days <- function(l,nb,epsilon) 1.5*((phi)/(1+phi))*(l/(epsilon*nb) - (1-phi)/(phi^2))
### end FUNCTIONS ###

if(!file.exists(outfile)){
	write(paste("Sample", "Lambda", "St.Dev", "NSeq", "NBases", "MeanHD", "MaxHD","Days(CI)", "Chi2","DF","Goodness_of_pval", sep="\t"), file=outfile, append=FALSE)
}
	
dlist <- read.table(file=infile, sep="\t", stringsAsFactors=F)

### calc HD with consensus
d0 <- dlist[which(dlist[,1]==dlist[1,1]),]
mult0 <- as.numeric(sub('.+_(\\d+)$', '\\1', d0[,2]))
nseq <- sum(mult0)
yvec0 <- rep(0, (1+max(d0[,3])))
for(i in 1:(1+max(d0[,3]))){ yvec0[i] <- sum(mult0[which(d0[,3]==(i-1))]) }

nl0 <- length(yvec0);
clambda <- sum((1:(nl0-1))*yvec0[-1])/sum(yvec0) #### THIS IS THE LAMBDA THAT FITS THE CONSENSUS ONLY DISTRIBUTION

### calc intersequence HD
d1 <- dlist[-which(dlist[,1]==dlist[1,1]),]	
yvec <- rep(0, (1+max(d1[,3])))
seqnames <- unique(d1[,1])
for(i in 1:length(seqnames)){
	tmp <- d1[which(d1[,1]==seqnames[i]),]
	m0 <- as.numeric(sub('.+_(\\d+)$', '\\1', tmp[1,1]))
	yvec[1] <- yvec[1] + 0.5*m0*(m0-1) ## 0 bin
	for(j in 1:dim(tmp)[1]){
		m1 <- as.numeric(sub('.+_(\\d+)$', '\\1', tmp[j,2]))
		val <- tmp[j,3]
		yvec[val+1] <- yvec[val+1] + m0*m1
	}
}

### Fitting

nl <- length(yvec)
lambda <- sum((1:(nl-1))*yvec[-1])/sum(yvec)
estdays <- days(lambda, nbases, epsilon)
		
print(paste("Estimated Lambda", format(lambda, digits=4), sep=" "))

#### U STAT ESTIMATE OF ST DEV
#### FORMULAE
#### Var(HD) = (N(N-1)/2)^(-1) (2(N-2)sigma1^2 + sigma2^2)
#### sigma1^2 = (N(N-1)(N-2)/3 -1)^(-1) sum_{i<j<l} ((Dij-mu)(Dil-mu)+(Dij-mu)(Djl-mu))
#### sigma2^2 = (N(N-1)/2-1)^(-1) sum_{i<j} (Dij-mu)^2

### construct a matrix of Dij's
### number of unique sequences
nuni <- dim(d0)[1]
TX <- matrix(rep(0,nuni^2), ncol=nuni)

for(i in 1:(dim(d0)[1]-1)){
	useq <- d0[i,2]
	TX[((i+1):dim(TX)[1]),i] <- d1[which(d1[,1]==useq),3]
}

sigma1 <- 0
sigma2 <- 0
muhat <- 0
denmu <- (nseq*(nseq-1)/2)^(-1)
den1 <- 12*(nseq*(nseq-1)*(nseq-2)*(nseq-3))^(-1)  
den2 <- den1/4   


for(n in 1:(nuni-1)){
	for(m in (n+1):nuni){
		muhat <- muhat + mult0[n]*mult0[m]*denmu*TX[m,n]
	}
}

for(n in 1:nuni){
	dnn <- 0
	sigma1 <- sigma1 + choose(mult0[n],3)*den1*2*(dnn-muhat)^2 
	sigma2 <- sigma2 + choose(mult0[n],2)*den2*(dnn-muhat)^2
	if(n != nuni){
		for(m in (n+1):nuni){
			dnm <- TX[m,n]
			dmm <- 0
			sigma2 <- sigma2 + mult0[n]*mult0[m]*(dnm - muhat)^2
			sigma1 <- sigma1 + (2/3)*choose(mult0[n],2)*mult0[m]*(dnm-muhat)*(dnm+2*dnn-3*muhat)
			sigma1 <- sigma1 + (2/3)*mult0[n]*choose(mult0[m],2)*(dnm-muhat)*(dnm+2*dmm-3*muhat)
			if(m != nuni){
				for(l in (m+1):nuni){
					dnl <- TX[l,n]
					dlm <- TX[l,m]
					sigma1 <- sigma1 + (2/3)*mult0[n]*mult0[m]*mult0[l]*((dnm-muhat)*(dnl-muhat)+(dnm-muhat)*(dlm-muhat)+(dnl-muhat)*(dlm-muhat)) 
				}
			}
		}
 	}
}

## varhd <- sqrt(denmu*(2*(nseq-2)*sigma1 + sigma2))
A <- 8/(nseq*(nseq-1)*(nseq-2)*(nseq-3))
B <- 4/(nseq*(nseq-1)*(nseq-2)*(nseq-3))
newvarhd <- sqrt(A*sigma1 + B*sigma2)
	
upplim <- days(lambda + 1.96*newvarhd, nbases, epsilon)
lowlim <- days(lambda - 1.96*newvarhd, nbases, epsilon)
uppdays <- round(upplim)
lowdays <- round(lowlim)
	
formatteddays <- paste(round(estdays), " (", lowdays, ", ", uppdays, ")", sep="") 

### output figures
dvec1 <- 0
for(i in 1:length(yvec)){ dvec1 <- c(dvec1, rep((i-1),yvec[i])) }
dvec1 <- dvec1[-1]

meanhd <- mean(dvec1)
maxhd <- max(dvec1)

figure <- paste(dir, sample, ".hd_freq_dist.jpg", sep="")
jpeg(file=figure, pointsize=14, width=800, height=800)
par(mar=c(5,5,4,2))

layout(matrix(c(1,2,3), nr=1, ncol=3), widths=c(5), heights=c(1), FALSE)
h <- hist(dvec1, breaks=seq(-1,max(dvec1),1), plot=FALSE)
hist(dvec1, breaks=seq(-0.5,max(dvec1)+0.5,1), freq=TRUE, xlab=sample, ylim = c(0,15+max(h$counts)), labels=ifelse(max(dvec1)<20,TRUE,FALSE), main=paste("Hamming Distance Frequency Distribution", sep=" "), cex.lab=2, cex.axis=2,cex.main=2)
lines(seq(0,max(dvec1)+1,1), 0.5*nseq*(nseq-1)*dpois(seq(0,1+max(dvec1),1), lambda=lambda),col="red", lwd=2)

dev.off()	
	
figure <- paste(dir, sample, ".hd_freq_dist.ps", sep="")
postscript(file=figure, width=400, height=300)
par(mar=c(5,5,4,3))

h <- hist(dvec1, breaks=seq(-1,max(dvec1),1), plot=FALSE)
hist(dvec1, breaks=seq(-0.5,max(dvec1)+0.5,1), freq=TRUE, xlab=sample, ylim = c(0,15+max(h$counts)),
labels=ifelse(max(dvec1)<20,TRUE,FALSE), main=paste("Hamming Distance Frequency Distribution", sep=" "), cex.lab=2, cex.axis=2,cex.main=2)
lines(seq(0,max(dvec1)+1,1), 0.5*nseq*(nseq-1)*dpois(seq(0,1+max(dvec1),1), lambda=lambda),col="red", lwd=2)

dev.off()


#### FIT THE CONSENSUS ONLY HD DISTRIBUTION
	
if (lambda!=0) {
			
	xvec1 <- c(yvec0, rep(0,nl0))
	yvec1 <- rep(0,2*nl0)
	yvec1[1] <- 1/2*yvec0[1]*(yvec0[1]-1)  ### freq at zero 
	mvals <- seq(2,2*nl0,2)

	for(m in mvals) {
		delta <- rep(0,m)
		delta[1+m/2] <- 1
		for(hj in 1:m){			
                        yvec1[m] <- yvec1[m] + 1/2*xvec1[hj]*xvec1[m-hj+1]
			if(m<2*nl0) { yvec1[m+1] <- yvec1[m+1] + 1/2*xvec1[hj]*(xvec1[m-hj+2] - delta[hj]) } 
		}
		if(m<2*nl0) { yvec1[m+1] <- yvec1[m+1] + 1/2*xvec1[m+1]*xvec1[1] }
	}
	
	dvec2 <- rep(0, yvec1[1])
	w <- which(yvec1>0)	
	for(hk in w[-1]) { dvec2 <- c(dvec2, rep((hk-1), yvec1[hk])) }
	
	mmax <- 1.5*(max(c(yvec0, yvec1)))			

	cfigure <- paste(dir, sample, ".conv_plot.jpg", sep="")
	jpeg(file=cfigure, pointsize=14, width=800, height=800)
	par(mar=c(5,5,4,2))

	layout(matrix(c(1,2,3), nr=1, ncol=3), widths=c(5), heights=c(1), FALSE)
	hh1 <- hist(dvec1, breaks=seq(-0.5,max(dvec1+0.5),1), freq=TRUE, xlab=sample, labels=ifelse(max(dvec1)<20,TRUE,FALSE), main=paste("Convolution Plot", sep=" "), col="blue",cex.lab=2, cex.axis=2,cex.main=2)
	lines(seq(0,(2*nl0-1),1), yvec1[1:length(seq(0,(2*nl0-1),1))], col="red", lwd=2)
	points(seq(0,(2*nl0-1),1), yvec1[1:length(seq(0,(2*nl0-1),1))], pch=23, col="red", lwd=2)
	legend("topright", legend=c("CONV","OBS"), fill=c("red","blue"), text.width=0.8)

	dev.off()
			
	figure <- paste(dir, sample, ".conv_plot.ps", sep="")
	postscript(file=figure, width=400, height=300)
	par(mar=c(5,5,4,3))
	
	hh1 <- hist(dvec1, breaks=seq(-0.5,max(dvec1+0.5),1), freq=TRUE, xlab=sample, labels=ifelse(max(dvec1)<20,TRUE,FALSE), main=paste("Convolution Plot", sep=" "), col="blue",cex.lab=2, cex.axis=2,cex.main=2)
	lines(seq(0,(2*nl0-1),1), yvec1[1:length(seq(0,(2*nl0-1),1))], col="red", lwd=2)
	points(seq(0,(2*nl0-1),1), yvec1[1:length(seq(0,(2*nl0-1),1))], pch=23, col="red", lwd=2)
	legend("topright", legend=c("CONV","OBS"), fill=c("red","blue"), text.width=0.8)

	dev.off()
		
	check <- 0 

	write(sample, file=outfile2, append=TRUE)
	write(paste("HD", "OBS", "CONV", sep="\t"), file=outfile2, append=TRUE)
	for(jj in 1:length(yvec1)){ 
		write(paste(jj-1, hh1$counts[jj], yvec1[jj], sep="\t"), file=outfile2, append=TRUE) 
		hey <- abs(hh1$counts[jj]-yvec1[jj])
		if(!is.na(hey)) { check <- check + hey }
	}

	check <- check/sum(yvec1)
	ifclause <- ifelse(check <= 0.1, "FOLLOWS", "DOES NOT FOLLOW")
	astring <- paste(sample, ifclause, "A STAR-PHYLOGENY", sep=" ")
	write(astring, file=outfile2, append=TRUE)
	write(" ", file=outfile2, append=TRUE)	
}
		

### CONSTRUCT SIGMA_ij MATRIX THEN INVERT IT
#pk <- function(x) ((nseq^2)*(2^x)*exp(-2*clambda)*(clambda^x))/factorial(x)
pk <- function(x) (exp( ( (log(nseq)*2)+(log(2)*x)+(-2*clambda)+log(clambda^x))-lfactorial(x) ) )
mui <- function(x) nseq*dpois(x, lambda=clambda)
SIGMA.DIM.MAX <- 170; # Beyond this value, factorial stops working in R.
sigma.dim <- min( SIGMA.DIM.MAX, (2*nl0) );
eyvec <- 0.5*pk(0:(sigma.dim-1))
eyvec[ !is.finite( eyvec ) ] <- 0; # Paul added this to avoid errors.  factorial(x) sometimes refuses to work (if x is too large).

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
    
	sdec <- La.svd(sigmaij)
	diag <- ifelse(sdec$d>1e-4,sdec$d,0)
	diagmat <- matrix(rep(0,sigma.dim^2), ncol=sigma.dim)
	for(ii in 1:sigma.dim){diagmat[ii,ii]<-ifelse(diag[ii]==0,0,1/diag[ii])}
	sigmainv <- sdec$u%*%diagmat%*%sdec$vt
	
	h <- hist(dvec1, breaks=seq(-1,max(dvec1),1), plot=FALSE)
	xvec <- h$breaks
	yvec <- h$counts
	nl1 <- length(yvec)
	aplambda <- sum((1:(nl1-1))*yvec[-1])/sum(yvec)
			
	pesce <- 0.5*nseq*(nseq-1)*dpois(0:(nl1-1), lambda=aplambda)

	if (length(yvec)<sigma.dim) { 
		ccvv <- sigma.dim - length(yvec) 
		yvec <- c(yvec, rep(0,ccvv)) 
	}
	
	chisq <- t(abs(yvec-eyvec))%*%sigmainv%*%(abs(yvec-eyvec))
	pval <- ifelse(chisq<0,2e-16,1-pchisq(chisq,df=nl0-1))		
	if(pval==0){ pval <- 2e-16 }
	if(chisq<0){ chisq <- NA }
} else { 
	chisq <- NA
	nl <- NA
	pval <- 0
}
       
write(paste(sample, format(lambda, digits=4), format(newvarhd, digits=4), nseq, nbases, format(meanhd, digits=2), maxhd, formatteddays, format(as.numeric(chisq), digits=4), nl0-1, format(as.numeric(pval), digits=4), sep="\t"), file=outfile, append=TRUE)

	









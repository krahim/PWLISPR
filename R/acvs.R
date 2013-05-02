# The following functions are contained in acvs.R
# These functions are loosley based on the lisp code provided with SAPA
# http://lib.stat.cmu.edu/sapaclisp/
# http://lib.stat.cmu.edu/sapaclisp/acvs.lisp

# acvs <- function(timeSeries, centreDataP=T, acsP=F)

# acvs_ws <- function(timeSeries, centreDataP=T, acsP=F)

# biasedACVStoUnbiasedACVS <- function(biasedACVS)

################################################################################

#libraries required
##library(waveslim);

################################################################################
# Requires:
##source("~/PWLisp/basicStatistics.R");
##source("~/PWLisp/utilities.R");
##source("~/PWLisp/tapers.R");
#source("nonparametric.R");  caused probs circular

################################################################################

#coded to match the lisp acvs function in the SAPA code
#Note: see Equations (191a) and (190b) of the SAPA book
acvs <- function(timeSeries, centreDataP=T, acsP=F) {
   N <- length(timeSeries);
   ntot <- nextPowerOf2(2*N);
   N_Xntot <- N*ntot;
   if(centreDataP) {
      timeSeries <- timeSeries - mean(timeSeries);
   }
   scratch <- c(timeSeries, array(0, ntot-N));
   scratch <- fft(scratch);
   scratch <- abs(scratch)^2;
   scratch <- fft(scratch);
   sampVar <- Re(scratch[1])/N_Xntot;
   scaleFactor <- if(acsP) {N_Xntot*sampVar;} else {N_Xntot;}
   result <- Re(scratch[1:N])/ scaleFactor;
   return(list(acvs=result, var=sampVar));
}

#uses my.acf provided in waveslim
acvs_ws <- function(timeSeries, centreDataP=T, acsP=F) {
   if(centreDataP) {
   timeSeries <- timeSeries - mean(timeSeries);
   }
   #requires waveslim
   res <- my.acf(timeSeries);
   if(acsP) {
      res <- res / res[1];
   }
   return(list(acvs=res, var=res[1]));
}

biasedACVStoUnbiasedACVS <- function(biasedACVS) {
   result <- biasedACVS*length(biasedACVS) / (length(biasedACVS):1);
   return(result);
}

## i'm not sure what this is for but acvs2 allows
## for matrix valued fft used in dpss.R
acvsB <- function(timeSeries, centreDataP=T, acsP=F) {
   N <- length(timeSeries);
   ntot <- nextPowerOf2(2*N);
   N_Xntot <- N*ntot;
   if(centreDataP) {
      timeSeries <- timeSeries - mean(timeSeries);
   }
   scratch <- c(timeSeries, array(0, ntot-N));
   scratch <- fft(scratch);
   scratch <- abs(scratch)^2;
   scratch <- fft(scratch);
   sampVar <- Re(scratch[1])/N_Xntot;
   scaleFactor <- if(acsP) {N_Xntot*sampVar;} else {N_Xntot;}
   result <- Re(scratch[1:N])/ scaleFactor;
   return(list(result=result, sampleVariance=sampVar));
}


# The following functions are contained in tapers.R
# These functions are loosley based on the lisp code provided with SAPA
# http://lib.stat.cmu.edu/sapaclisp/
# http://lib.stat.cmu.edu/sapaclisp/tapers.lisp

# cosineDataTaper <- function(  timeSeries,
#                              taperParameter=1.0,
#                              normalization="unity")

# HanningDataTaper <- function( timeSeries,
#                              taperParameter=NaN,
#                              normalization="unity")

# dpssDataTaper <- function( timeSeries, nw=4, normalization="unity")

# taperTimeSeries <- function( centredTimeSeries,
#                              dataTaperFn=NULL,
#                              taperParameter=NULL,
#                              recentreAfterTaperingP=T,
#                              restorePowerOptionP=T)
                              
# taperTimeSeriesWSupplied <- function( centredTimeSeries,
#                              dataTaperParameters,
#                              recentreAfterTaperingP=T,
#                              restorePowerOptionP=T)

################################################################################
# Requires:

################################################################################

cosineDataTaper <- function(  timeSeries,
                              taperParameter=1.0,
                              normalization="unity") {
   n <- length(timeSeries);
   sumOfSq <- 0.0;
   C_hTop <- 0.0;
   lowCutoff <- floor(taperParameter*n)/2;
   highCutoff <- (1+n) - lowCutoff;
   twoPiDenom <- 2*pi/(1+ floor(taperParameter*n));
   result <- array(0.0, n);
   for( i in 1:(n/2) ) {
      xfactor <- 1.0;
      if( 0 <= (i-1) && (i-1) <= lowCutoff) {
         xfactor <-  0.5 * (1- cos(twoPiDenom * i));
      }
      else if( lowCutoff < (i-1) && (i-1) < highCutoff) {
         xfactor <- 1.0;
      }
      else {
        cat("Fatal Error in Data Taper\n");
        return();
      }
      sumOfSq <- sumOfSq + (2*xfactor^2);
      C_hTop <- C_hTop + (2*xfactor^4);
      result[i] <- xfactor * timeSeries[i];
      result[n-i+1] <- xfactor * timeSeries[n-i+1];
   }
   if(n%%2!=0) {
      mid <- trunc(n/2)+1;
      xFactor <- 1.0;
      if( 0 <= (mid-1) && (mid-1) <= lowCutoff) {
         xfactor <-  0.5 * (1- cos(twoPiDenom * mid));
      }
      else if( lowCutoff < (i-1) && (i-1) < highCutoff) {
         xfactor <- 1.0;
      }
      else {
        cat("Fatal Error in Data Taper\n");
        return();
      }
      sumOfSq <- sumOfSq + xfactor^2;
      C_hTop <- C_hTop + xfactor^4;
      result[mid] <- xfactor * timeSeries[mid];
   }
   xfactor <- 1.0;
   if(normalization=="unity") {
      xfactor <- sqrt(1/sumOfSq);
   } else {
      xfactor <- sqrt(n/sumOfSq);
   }

   result <- result*xfactor;
   #   [2] C_h, the variance inflation factor
   #       computed using Equation (251b) in the SAPA book"
   C_h <-  n * C_hTop / sumOfSq^2;

   return(list(taperedTS=result,C_h=C_h));
}

#for(N in c(100,1000) ) {
#   for(p in c(0.0, 0.2, 0.5, 1.0)) {
#      res <- cosineDataTaper(array(1.0, N), taperParameter=p);
#      cat(paste(";;; p =  ",
#         format(p, nsmall=2) , ", N = ",
#         format(N, width=4), ",    ",
#         format(res$C_h, nsmall=6), ",   ",
#         format(sum(res$taperedTS^2), nsmall=6),
#         "\n"));
#   }
#}
#;;; p =   0.00 , N =   100 ,     1.020408 ,    1.000000
#;;; p =   0.20 , N =   100 ,     1.110349 ,    1.000000
#;;; p =   0.50 , N =   100 ,     1.338242 ,    1.000000
#;;; p =   1.00 , N =   100 ,     1.925193 ,    1.000000
#;;; p =   0.00 , N =  1000 ,     1.002004 ,    1.000000
#;;; p =   0.20 , N =  1000 ,     1.115727 ,    1.000000
#;;; p =   0.50 , N =  1000 ,     1.346217 ,    1.000000
#;;; p =   1.00 , N =  1000 ,     1.942502 ,    1.000000

HanningDataTaper <- function( timeSeries,
                              taperParameter=NaN,
                              normalization="unity") {
   #calls cosine-data-taper! with taper-parameter set to 1.0
   return(cosineDataTaper(timeSeries, taperParameter=1.0, normalization));
}

#for(N in c(100,1000) ) {
#   res <- HanningDataTaper(array(1.0, N));
#   cat(paste(";;; N = ",
#      format(N, width=4), ",    ",
#      format(res$C_h, nsmall=6), ",   ",
#      format(sum(res$taperedTS^2), nsmall=6),
#      "\n"));
#
#}
#;;; N =   100 ,     1.925193 ,    1.000000
#;;; N =  1000 ,     1.942502 ,    1.000000

dpssDataTaper <- function( timeSeries, taperParameter=4, normalization="unity") {
    ## (defun dpss-data-taper!
    ##   "given
    ##    [1] a-time-series (required)
    ##        <=> a vector containing a time series;
    ##            this vector is modified UNLESS keyword result
    ##            is bound to a different vector
    ##    [2] taper-parameter (keyword; 4.0)
    ##        ==> a number > 0.0 that specifies NW,
    ##            the product of the sample size and
    ##            the half-width of the concentration interval
    ##            (see SAPA, page 211)
    ##    [3] normalization (keyword; :unity)
    ##        ==> if :unity, taper normalized such that
    ##            its sum of squares is equal to unity (Equation (208a));
    ##            if :N, taper normalized such that
    ##            its sum of squares is equal to the number of points
    ##            in the time series
    ##    [4] result (keyword; a-time-series)
    ##        <=> a vector to contain time series multiplied
    ##            by cosine data taper
    ## returns
    ##    [1] a vector containing the tapered time series
    ##    [2] C_h, the variance inflation factor
    ##        computed using Equation (251b) in the SAPA book"
    
   N <- length(timeSeries);
   sumOfSq <- 0.0;
   C_hTop <- 0.0;
   W <- taperParameter/N;
   betaPi <- pi*W*(N-1);
   result <- array(0.0, N);
   for( i in 1:(N/2) ) {
      tmp <- betaPi * sqrt(1- ((as.double(1+2*(i-1))/N)-1)^2);
      xfactor <- besselI(tmp,0);
      sumOfSq <- sumOfSq + (2*xfactor^2);
      C_hTop <- C_hTop + (2*xfactor^4);
      result[i] <- xfactor * timeSeries[i];
      result[N-i+1] <- xfactor * timeSeries[N-i+1];
   }
   if(N%%2!=0) {
      mid <- trunc(N/2)+1;
      tmp <- betaPi * sqrt(1- ((as.double(1+2*(mid-1))/N)-1)^2);
      xfactor <- besselI(tmp,0);
      sumOfSq <- sumOfSq + xfactor^2;
      C_hTop <- C_hTop + xfactor^4;
      result[mid] <- xfactor * timeSeries[mid];
   }

   xfactor <- 1.0;
   if(normalization=="unity") {
      xfactor <- sqrt(1/sumOfSq);
   } else {
      xfactor <- sqrt(N/sumOfSq);
   }

   result <- result*xfactor;
   #   [2] C_h, the variance inflation factor
   #       computed using Equation (251b) in the SAPA book"
   C_h <-  N * C_hTop / sumOfSq^2;

   return(list(taperedTS=result,C_h=C_h));
}

## dpssDataTaper(rep(1, 100), 1)$C_h
## sum(dpssDataTaper(rep(1, 100), 1)$taperedTS**2)

## dpssDataTaper(rep(1, 100), 2)$C_h
## sum(dpssDataTaper(rep(1, 100), 2)$taperedTS**2)

## dpssDataTaper(rep(1, 100), 4)$C_h
## sum(dpssDataTaper(rep(1, 100), 4)$taperedTS**2)

## dpssDataTaper(rep(1, 100), 8)$C_h
## sum(dpssDataTaper(rep(1, 100), 8)$taperedTS**2)

## dpssDataTaper(rep(1, 1000), 1)$C_h
## sum(dpssDataTaper(rep(1, 1000), 1)$taperedTS**2)

## dpssDataTaper(rep(1, 1000), 2)$C_h
## sum(dpssDataTaper(rep(1, 1000), 2)$taperedTS**2)

## dpssDataTaper(rep(1, 1000), 4)$C_h
## sum(dpssDataTaper(rep(1, 1000), 4)$taperedTS**2)

## dpssDataTaper(rep(1, 1000), 8)$C_h
## sum(dpssDataTaper(rep(1, 1000), 8)$taperedTS**2)

#note different beteween calculated values and those in
#SAPA table 248.
#According to the lisp docs...
#;;; These are slightly off from Table 248 of SAPA
#;;; (evidently computed using a less accurate approximation
#;;; to 0th order dpss):
#for(N in c(100,1000) ) {
#   for(NW in c(1.0, 2.0, 4.0, 8.0)) {
#      res <- dpssDataTaper(array(1.0, N), taperParameter=NW);
#      cat(paste(";;; NW =  ",
#         format(NW, nsmall=2) , ", N = ",
#         format(N, width=4), ",    ",
#         format(res$C_h, nsmall=6), ",   ",
#         format(sum(res$taperedTS^2), nsmall=6),
#         "\n"));
#   }
#}
#;;; NW =   1.00 , N =   100 ,     1.404251 ,    1.000000
#;;; NW =   2.00 , N =   100 ,     1.995994 ,    1.000000
#;;; NW =   4.00 , N =   100 ,     2.820333 ,    1.000000
#;;; NW =   8.00 , N =   100 ,     3.984642 ,    1.000000
#;;; NW =   1.00 , N =  1000 ,     1.410677 ,    1.000000
#;;; NW =   2.00 , N =  1000 ,     2.005059 ,    1.000000
#;;; NW =   4.00 , N =  1000 ,     2.833080 ,    1.000000
#;;; NW =   8.00 , N =  1000 ,     4.002674 ,    1.000000

#The lisp function center&taper-time-series is replaced
#with two functions for tapering, the date is assumed
#centred before the taper functions are called.
#The two taper functions are taperTimeSeries and taperTimeSeriesWSupplied.
#tapers time series using a taper function with parameter
#P&W lisp has centre & taper time series but it is trivial to
#centre the data in R
taperTimeSeries <- function( centredTimeSeries,
                              dataTaperFn=NULL,
                              taperParameter=NULL,
                              recentreAfterTaperingP=TRUE,
                              restorePowerOptionP=TRUE) {
   if(!is.function(dataTaperFn)) {
      return();
   }
   if(is.null(taperParameter) && !identical(dataTaperFn, HanningDataTaper)) {
      return();
   }

   N <- length(centredTimeSeries);
   sumOfSquaresCentredTS <- t(centredTimeSeries)%*%centredTimeSeries;

   #currently normalize for N...

   #we normalize such that Sum h_t^2 = N
   #to compensate we now use a (delta T / N) fiddle factor for the sdf
   #like in the periodogram--see page 207 SAPA for how when a direct
   #spectral est is the same as the periodogram.
   tapers <- dataTaperFn(  centredTimeSeries,
                           taperParameter=taperParameter,
                           normalization="N");

   if(recentreAfterTaperingP) {
      tapers$taperedTS <- tapers$taperedTS - mean(tapers$taperedTS);
   }
   if(restorePowerOptionP) {
      multFactor <- as.double(sqrt(sumOfSquaresCentredTS /
            t(tapers$taperedTS)%*%tapers$taperedTS));
      tapers$taperedTS <- tapers$taperedTS * multFactor;
   }
   return(list(result=tapers$taperedTS, C_h=tapers$C_h));
}

#cts <- c(71, 63, 70, 88, 99, 90, 110)
#cts <- cts - mean(cts)
#taperTimeSeries(cts, dataTaper=dpssDataTaper, taperParameter=1.0);
#tapers <- taperTimeSeries(cts,
#               dataTaper=dpssDataTaper,
#               taperParameter=1.0,
#               recentreAfterTaperingP=F);
#tapers;
#mu <- mean(tapers$result);
#mu;
#v <- var(tapers$result)*(length(tapers$result)-1)/length(tapers$result);
#v;
#v+(mu)^2;
#
#tapers <- taperTimeSeries(cts,
#               dataTaper=dpssDataTaper,
#               taperParameter=1.0,
#               recentreAfterTaperingP=F,
#               restorePowerOptionP=F);
#tapers;
#mu <- mean(tapers$result);
#mu;
#v <- var(tapers$result)*(length(tapers$result)-1)/length(tapers$result);
#v;
#v+(mu)^2;

#The lisp function center&taper-time-series is replaced
#with two functions for tapering, the date is assumed
#centred before the taper functions are called.
#The two taper functions are taperTimeSeries and taperTimeSeriesWSupplied.
#tapers time series using a supplied data taper.
#P&W lisp has centre & taper time series but it is trivial to
#centre the data in R
taperTimeSeriesWSupplied <- function( centredTimeSeries,
                              dataTaperParameters,
                              recentreAfterTaperingP=TRUE,
                              restorePowerOptionP=TRUE) {
   if(length(centredTimeSeries) != length(dataTaperParameters)) {
      return();
   }
   N <- length(centredTimeSeries);
   sumOfSquaresCentredTS <- t(centredTimeSeries)%*%centredTimeSeries;
   result <- centredTimeSeries * dataTaperParameters;

   #the following is done from the lisp function supplied-data-taper!

   #currently normalize for N...
   sumOfSquaresTapers <- t(dataTaperParameters)%*%dataTaperParameters;
   #sumOfSquaresTapers is 1 for dpss tapers
   nfactor <- sqrt(N/as.double(sumOfSquaresTapers));

   #we normalize such that Sum h_t^2 = N
   #to compensate we now use a (delta T / N) fiddle factor for the sdf
   #like in the periodogram--see page 207 SAPA for how when a direct
   #spectral est is the same as the periodogram

   result <- nfactor * result;
   #C_h is the variance inflation factor discussed on page 248 of SAPA
   #it is computed from eqn 251b as per the lisp fn
   #supplied-data-taper in tapers.lisp

   C_h <- N*sum(dataTaperParameters^4)/(sumOfSquaresTapers^2);

   #end of supplied-data-taper! call

   if(recentreAfterTaperingP) {
      result <- result - mean(result);
   }
   if(restorePowerOptionP) {
      multFactor <- sqrt(sumOfSquaresCentredTS / t(result)%*%result);
      result <- result*multFactor;
   }
   return(list(result=result, C_h=C_h));
}

## # The following functions are contained in harmonic.R
## # These functions are loosley based on the lisp code provided with SAPA
## # http://lib.stat.cmu.edu/sapaclisp/
## # http://lib.stat.cmu.edu/sapaclisp/harmonic.lisp

## #FishergStat <- function(timeSeries,
## #                     centreData=T,
## #                     alpha=.05)

## # periodogramAtOneFrequency <- function(timeSeries,
## #                       freq,
## #                       centreData=T,
## #                       samplingTime=1.0)

## ################################################################################

## # Requires:
## ##source("utilities.R");
## ##source("nonparametric.R");

## ################################################################################

FishergStat <- function(timeSeries,
                     centreData=T,
                     alpha=.05) {
    
    ##         "given
    ##    [1] time-series (required)
    ##        ==> a vector of time series values
    ##    [2] center-data (keyword; t)
    ##        ==> if t, subtract sample mean from time series;
    ##            if a number, subtracts that number from time series;
    ##            if nil, time-series is not centered
    ##    [3] start (keyword; 0)
    ##        ==> start index of vector to be used
    ##    [4] end (keyword; length of time-series)
    ##        ==> 1 + end index of vector to be used
    ##    [5] alpha (keyword; 0.05)
    ##        ==> critical level at which Fisher's g test
    ##            is performed
    ## returns
    ##    [1] Fisher's g statistic
    ##    [2] approximation to critical level of Fisher's g test
    ##        from Equation (491c)
    ##    [3] either :reject or :fail-to-reject, depending on whether
    ##        or not we reject or fail to reject the null hypothesis
    ##        of white noise
    ## ---
    ## Note: see Section 10.9 of the SAPA book"
    result <- NULL;
    N <- length(timeSeries);
    M <- length(timeSeries);
    m <- 0;
    if(N%%2==0) {
        m <- (N-2)/2;
    } else {
        m <- (N-1)/2;
    }
    
    if(m<= 1) {
        return(list(gStat=1.0, critical=1.0, hyp="Fail To Reject"));
    }
    gF <- 1.0 - ( (alpha/m)^(1/(m-1)));
    thePeriodogram <- periodogram(timeSeries,
                                  centreData,
                                  nNonZeroFreqs="Fourier",
                                  returnEstFor0FreqP=F,
                                  sdfTransformation=F);
    maxShp <- max(thePeriodogram$sdf[1:m]);
    sumShp <- sum(thePeriodogram$sdf[1:m]);
    gStatistic <- maxShp/as.double(sumShp);
    if(gStatistic > gF) {
        result <- "Reject";
    } else {
        result <- "Fail to Reject";
    }
    return(list(gStatistic=gStatistic,
                gF=gF,
                result=result));
}


periodogramAtOneFrequency <- function(timeSeries, 
                       freq, 
                       centreData=T,
                       samplingTime=1.0) {
    ## #   "given
    ## #    [1] time-series (required)
    ## #        ==> a vector of time series values
    ## #    [2] center-data (keyword; t)
    ## #        ==> if t, subtract sample mean from time series;
    ## #            if a number, subtracts that number from time series;
    ## #            if nil, time-series is not centered
    ## #    [3] start (keyword; 0)
    ## #        ==> start index of vector to be used
    ## #    [4] end (keyword; length of time-series)
    ## #        ==> 1 + end index of vector to be used
    ## #    [5] sampling-time (keyword; 1.0)
    ## #        ==> sampling time (called delta t in the SAPA book)
    ## # returns
    ## #    [1] value of periodogram at freq
    ## #    [2] approximate conditional least squares estimate for A,
    ## #        the amplitude of cosine term in Equation (461a)
    ## #        of the SAPA book
    ## #    [3] approximate conditional least squares estimate for B
    ## #        the amplitude of sine term
    ## # ---
    ## # Note: see Section 10.2 of the SAPA book"
    
   N <- length(timeSeries);
   j <- 0;
   centreFactor <- 
      if(is.numeric(centreData)) centreData 
         else if(centreData) mean(centreData)
         else 0;
   cosSum <- 0;
   sinSum <- 0;
   twoPiFDeltT <- 2*pi*freq*samplingTime;
   tIndex <- 1;
   for(i in 0:(N -1)) {
      cosSum <- cosSum + (timeSeries[j +1] - centreFactor) *
         cos(twoPiFDeltT*tIndex);
      sinSum <- sinSum + (timeSeries[j +1] - centreFactor) *
         sin(twoPiFDeltT*tIndex);
      j <- j +1;
      tIndex <- tIndex +1;   
   }
   return(list(valueAtFreq=(samplingTime/N)*(cosSum^2+sinSum^2),
      AEst=(2.0/N)*cosSum,
      BEst=(2.0/N)*sinSum));
}
##
## #X <- c(71.0,  63.0,  70.0,  88.0,  99.0,  90.0, 110.0, 135.0, 128.0, 154.0, 156.0, 141.0, 131.0, 132.0, 141.0, 104.0, 136.0, 146.0, 124.0, 129.0)
## #periodogramAtOneFrequency(X, .2, samplingTime = .25)

## moved from analysis.R Not sure if thes belong in the PWLISPR package
## may move to multitaper Development


##### Sine multitaper Functions #####

## # sineMTM <- function(  timeSeries,
## #                  k,
## #                  deltaT=1.0,
## #                  frequencies="nextPowerOf2",
## #                  sdfTransformation=convertTodB)
                  
## # sineMultitaperSpectralEstimate <- function(timeSeries,
## #                     mDataTapers,
## #                     centreData=T,
## #                     nTapers=length(mDataTapers[1,]),
## #                     nNonZeroFreqs="halfNextPowerOf2",
## #                     returnEstFor0FreqP=F,
## #                     samplingTime=1.0,
## #                     recentreAfterTaperingP=T,
## #                     restorePowerOptionP=T,
## #                     sdfTransformation=convertTodB)

## # sine.taper.2 <- function(n,k)


sine.taper.2 <- function(n,k) {
   tapers <- matrix(nrow=n, ncol=k);
   for( i in 1:k) {
      tapers[,i] <- sqrt(2/(n+1))*sin((pi * i * 1:n)/(n + 1)) ;
   }
   return(tapers);
}

sineMultitaperSpectralEstimate <- function(timeSeries,
                     mDataTapers,
                     centreData=T,
                     nTapers=length(mDataTapers[1,]),
                     nNonZeroFreqs="halfNextPowerOf2",
                     returnEstFor0FreqP=F,
                     samplingTime=1.0,
                     recentreAfterTaperingP=T,
                     restorePowerOptionP=T,
                     sdfTransformation=convertTodB)  {
   sampleSize <- length(timeSeries);
   nFreqs <- getNFreqs(nNonZeroFreqs, sampleSize, returnEstFor0FreqP);
   nDFT <- getNDFT(nNonZeroFreqs, sampleSize);
   offSet <- if(returnEstFor0FreqP) 0 else 1;
   fiddleFactorSDF <- samplingTime / as.double(sampleSize);
   fiddleFactorFreq <- 1 / as.double(nDFT * samplingTime);
   if(centreData) {
      timeSeries <- timeSeries - mean(timeSeries);
   }
   mTaperedSpectra <- matrix(nrow=nFreqs, ncol=nTapers);
   for (k in 1:nTapers) {
      tSDF <- taperTimeSeriesWSupplied(timeSeries,
                     mDataTapers[,k],
                     recentreAfterTaperingP = recentreAfterTaperingP,
                     restorePowerOptionP = restorePowerOptionP)$result;
      tSDF <- c(tSDF, rep(0.0, nDFT-sampleSize));
      tSDF <- fft(tSDF);
      mTaperedSpectra[,k] <-
         abs(tSDF[(1+offSet):(nFreqs+offSet)])^2*fiddleFactorSDF;

   }
   parabolicWeights <- 1-((0:(nTapers-1))^2/nTapers^2);
   resultSDF <- as.vector(mTaperedSpectra %*% parabolicWeights) /
      sum(parabolicWeights);

   if(is.function(sdfTransformation)) {
      resultSDF <- sdfTransformation(resultSDF);
   }
   resultFreqs <- ((0+offSet):(nFreqs+offSet-1))*fiddleFactorFreq;
   return(list(   resultSDF=resultSDF,
                  resultFreqs=resultFreqs,
                  nFreqs=nFreqs,
                  mTaperedSpectra=mTaperedSpectra));

}

#x <- scan("http://faculty.washington.edu/dbp/DATA/Willamette.dat", skip=8);
#Smtm <- sineMultitaperSpectralEstimate(x, sine.taper.2(395,4), samplingTime=1/12);

sineMTM <- function(  timeSeries,
                    k,
                    deltaT=1.0,
                    frequencies="nextPowerOf2",
                    restorePower=F,
                    sdfTransformation=convertTodB) {
    N <- length(timeSeries);
    mDataTapers <- sine.taper.2(N, k);
    res <- sineMultitaperSpectralEstimate(
        timeSeries,
        mDataTapers,
        samplingTime=deltaT,
        nNonZeroFreqs=frequencies,
        restorePowerOptionP=restorePower,
        sdfTransformation=sdfTransformation);
    
    
    return( list(  resultSDF=res$resultSDF,
                 resultFreqs=res$resultFreqs,
                 mTaperedSpectra=res$mTaperedSpectra))##,
                 ## seems like the following was an error May 13, 2013
                 ##eigenValues=tapers$eigen));
    
}

#Smtm4 <- sineMTM(x, 4);
#plot(Smtm4$resultFreqs, Smtm4$resultSDF, type="l") ;

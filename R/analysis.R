# The following functions are contained in analysis.R
# These functions are designed to work with the
# functions based on the lisp code provided with SAPA
# http://lib.stat.cmu.edu/sapaclisp/

#these functions were not included in
#and do not contain similar functions in the SAPA code

# ThomsonFTest <- function(mEigenCoef,
#                     mDataTapers,
#                     deltaT=1.0,
#                     nFreqs=length(mEigenCoef[,1]),
#                     nTapers=length(mDataTapers[1,]))

# mtm <- function(  timeSeries,
#                  k,
#                  deltaT=1.0,
#                  nw=4,
#                  frequencies="nextPowerOf2",
#                  adaptive=T,
#                  returnAdaptiveCIs=F,
#                  recentreAfterTaper=F,
#                  restorePower=F,
#                  sdfTransformation=convertTodB,
#                  ThomsonFTest=T,
#                  returnCoherenceData=F) 



################################################################################
# Requires:
## source("utilities.R");
## source("utilities2.R");
## source("basicStatistics.R");
## source("tapers.R");
## source("nonparametric.R");
## source("multitaper.R");

################################################################################

#Thomson F-test based on SAPA eqn 499a and 499c
#details on SAPA pages 496-499
ThomsonFTest <- function(mEigenCoef,
                     mDataTapers,
                     deltaT=1.0,
                     nFreqs=length(mEigenCoef[,1]),
                     nTapers=length(mDataTapers[1,])) {

   if( !is.matrix(mEigenCoef) ||
       !is.matrix(mDataTapers)) {
      return();
   }

   H0 <- apply(mDataTapers, 2, sum);
   ##zero H0 where theoretically zero
   H0[1:(nTapers/2)*2] <- 0.0;
   ##H0 <- H0 

   #page 499 eqns can be manipulated to eliminate the need for deltaT
   #except in H0

   #manually iterate for k=1,3,5,... (note: equivalent to lisp SAPA
   # based k=0,2,4,...  page 499 eqn 499a
   #for loops are equivalent to the following two statements
   #SumH0Sq <- t(H0) %*% H0;
   #C_hatXSumH0SqDivSqrtDeltaT <- mEigenCoef %*% H0;

   SumH0Sq <- 0.0;
   C_hatXSumH0SqDivSqrtDeltaT <- array(0.0, nFreqs);
   for( k in 1:((nTapers+1)/2)*2-1 ) {
      SumH0Sq <- SumH0Sq + H0[k]^2;
      for( f in 1:nFreqs ) {
         C_hatXSumH0SqDivSqrtDeltaT[f] <-
            C_hatXSumH0SqDivSqrtDeltaT[f] + (mEigenCoef[f,k] * H0[k]);
      }
   }

   C_hatDivSqrtDeltaT <- C_hatXSumH0SqDivSqrtDeltaT / SumH0Sq;
   ##mEigenCoefEst <- t(H0%*%t(C_hatDivSqrtDeltaT));
   ##mEigenCoefEst <- C_hatDivSqrtDeltaT %*% t(H0)
   mEigenCoefEst <- tcrossprod(C_hatDivSqrtDeltaT, H0)
   FdenDivDeltaT <- apply(abs(mEigenCoef - mEigenCoefEst)^2, 1, sum);
   Fstat <- (nTapers-1) * abs(C_hatXSumH0SqDivSqrtDeltaT)^2 /
            (FdenDivDeltaT * SumH0Sq);

   return(list(FTest=Fstat,
               numDf=2,
               denDf=2*nTapers-2,
               mC_hatDivSqrtDeltaT=C_hatDivSqrtDeltaT));
}

#equation 500 S_n for simple multitaper

mtm <- function(  timeSeries,
                  k,
                  deltaT=1.0,
                  nw=4,
                  frequencies="nextPowerOf2",
                  adaptive=T,
                  returnAdaptiveCIs=F,
                  recentreAfterTaper=F,
                  restorePower=F,
                  sdfTransformation=convertTodB,
                  ThomsonFTest=T,
                  returnCoherenceData=F) {
   cis <- list(lowerCI=NULL, upperCI=NULL);
   sigma2 = sampleVariance.bias(timeSeries);
   tapers <- dpss.taper.2(length(timeSeries), k, nw=nw);
   result <- multitaperSpectralEstimate(timeSeries,
                     tapers$v,
                     nNonZeroFreqs=frequencies,
                     samplingTime=deltaT,
                     recentreAfterTaperingP=recentreAfterTaper,
                     restorePowerOptionP=restorePower,
                     sdfTransformation=F,
                     returnEigenCoef=ThomsonFTest);

   FdistrRes <- ThomsonFTest(  result$mEigenCoef,
                           tapers$v,
                           deltaT);
   if(adaptive) {
      aResult <-  eigenSpectraToAdaptiveMultitaperSpectralEstimate(
                           result$mEigenSpectra,
                           tapers$eigen,
                           sigma2,
                           sdfTransformation=F,
                           returnWeights=returnCoherenceData);

      if(is.function(sdfTransformation)) {
         aResult$resultSDF <- sdfTransformation(aResult$resultSDF);
         if(returnAdaptiveCIs) {
            cis <- createCIforAMTsdfEst( aResult$resultSDF,
                     aResult$resultDF);
         }
      }

      coherenceData = NULL;
      if(returnCoherenceData) {
         coherenceData <- list(  mEigenCoef=result$mEigenCoef,
                              mWeights=aResult$mWeights,
                              evals=tapers$eigen);
      }
      return( list(  resultSDF=aResult$resultSDF,
                     resultFreqs=result$resultFreqs,
                     df=aResult$resultDF,
                     upperCI=cis$upperCI,
                     lowerCI=cis$lowerCI,
                     coherenceData=coherenceData,
                     mEigenSpectra=result$mEigenSpectra,
                     eigenValues=tapers$eigen,
                     FdistrRes=FdistrRes));
   }

   resSDF <-  eigenspectraToMultitaperSpectralEstimate(
                        result$mEigenSpectra,
                        sdfTransformation=F);
   if(is.function(sdfTransformation)) {
      resSDF <- sdfTransformation(resSDF);
   }

   return( list(  resultSDF=resSDF,
                  resultFreqs=result$resultFreqs,
                  mEigenSpectra=result$mEigenSpectra,
                  eigenValues=tapers$eigen,
                  FdistrRes=FdistrRes));

}

#Willamette figure 512 test
#x <- scan("http://faculty.washington.edu/dbp/DATA/Willamette.dat", skip=8);
#split.screen(figs=c(2,1));
#screen(1);
#Sxx.mtm5 <- mtm(x, 5, adaptive=F, deltaT=1/12, frequencies=512);
#plot(Sxx.mtm5$resultFreqs, Sxx.mtm5$resultSDF, ylim=c(-40,20), type="l", xlab="f (cycles/year) ", ylab="dB");
#title("Willamette Multitaper K=5 NW=4 number of sdf freqs=512");
#screen(2);
#plot(Sxx.mtm5$resultFreqs, Sxx.mtm5$FdistrRes$FTest, ylim=c(0, 50), type="l", xlab="f (cycles/year) ", ylab="F Test");
#lines(Sxx.mtm5$resultFreqs, array(18.5, 512));
#lines(Sxx.mtm5$resultFreqs, array(8.6, 512));
#title("Willamette Multitaper F-Test K=5 NW=4 number of sdf freqs=512");
#close.screen(all=T);

#abs(aResult$mWeights * result$mEigenCoef)^2 %*% (sqrt(dpssNW4$eigen) * sqrt(dpssNW4$eigen)) / aResult$mWeights^2 %*% dpssNW4$eigen
#((aResult$mWeights * result$mEigenCoef) * Conj(aResult$mWeights * result$mEigenCoef))  %*% (sqrt(dpssNW4$eigen) * sqrt(dpssNW4$eigen)) / aResult$mWeights^2 %*% dpssNW4$eigen



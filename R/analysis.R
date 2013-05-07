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

##### Sine multitaper Functions #####

# sineMTM <- function(  timeSeries,
#                  k,
#                  deltaT=1.0,
#                  frequencies="nextPowerOf2",
#                  sdfTransformation=convertTodB)
                  
# sineMultitaperSpectralEstimate <- function(timeSeries,
#                     mDataTapers,
#                     centreData=T,
#                     nTapers=length(mDataTapers[1,]),
#                     nNonZeroFreqs="halfNextPowerOf2",
#                     returnEstFor0FreqP=F,
#                     samplingTime=1.0,
#                     recentreAfterTaperingP=T,
#                     restorePowerOptionP=T,
#                     sdfTransformation=convertTodB)

# sine.taper.2 <- function(n,k)

#coherenceAdaptiveJackKnifeMatrix <- function(mSpectra1,
#                              mSpectra2,
#                              mEigenCoef1,
#                              mEigenCoef2,
#                              mWeights1,
#                              mWeights2,
#                              eigenVals)

# coherenceAdaptive <- function(Samt1,
#                              Samt2,
#                              mEigenCoef1,
#                              mEigenCoef2,
#                              mWeights1,
#                              mWeights2,
#                              eigenVals1,
#                              eigenVals2)

#coherence <- function(  Smt1,
#                        Smt2)

#mtmCohere2 <- function( timeSeriesA,
#                        timeSeriesB,
#                        k,
#                        deltaT=1.0,
#                        nw=4,
#                        frequencies="nextPowerOf2",
#                        sdfTransformation=convertTodB,
#                        N=length(timeSeriesA),
#                        returnJackKnifeEst=T)

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
                 mTaperedSpectra=res$mTaperedSpectra,
                 eigenValues=tapers$eigen));
    
}

#Smtm4 <- sineMTM(x, 4);
#plot(Smtm4$resultFreqs, Smtm4$resultSDF, type="l") ;

coherenceAdaptive <- function(Samt1,
                              Samt2,
                              mEigenCoef1,
                              mEigenCoef2,
                              mWeights1,
                              mWeights2,
                              eigenVals) {
   #TODO: checks required
   f <- length(Samt1);
   k <- length(mEigenCoef1[1,]);

   num <- (((mWeights1 * mEigenCoef1) * Conj(mWeights2 * mEigenCoef2))   %*%
         eigenVals) /
         ((mWeights1*mWeights2) %*% eigenVals);
   den <- Samt1*Samt2;

   return(list(magSq=abs(as.vector(num))^2 /den, phase=Arg(as.vector(num)/sqrt(den))));
}


#returns k jacknifed mag Sq coh matrix.
#the k colmuns in this matrix is to be averaged with
#the original coherence est to obtain a jack knifed est
coherenceAdaptiveJackKnifeMatrix <- function(mSpectra1,
                              mSpectra2,
                              mEigenCoef1,
                              mEigenCoef2,
                              mWeights1,
                              mWeights2,
                              eigenVals) {
   #TODO: checks required
   f <- length(mEigenCoef1[,1]);
   nTapers <- length(mEigenCoef1[1,]);

   mJackKnifeCoherenceEst <- matrix(0.0, f, nTapers+1);
   mJackKnifePhaseEst <- matrix(0.0, f, nTapers+1);

   mWeights1Sq <- mWeights1^2;
   mWeights2Sq <- mWeights2^2;
   mWeight1XWeight2 <- mWeights1*mWeights2;
   weightedSpectra1 <- mWeights1Sq * mSpectra1;
   weightedSpectra2 <- mWeights2Sq * mSpectra1;
   weightedCoef1 <- mWeights1 * mEigenCoef1;
   conjWeightedCoef2 <- Conj(mWeights2 * mEigenCoef2);

   for( k in 1:nTapers) {
      eval_k <- eigenVals[k];
      eigenVals[k] <- 0.0;

      num <- as.vector((weightedCoef1 * conjWeightedCoef2)   %*%
            eigenVals /
            (mWeight1XWeight2 %*% eigenVals));

      Samt1 <- (weightedSpectra1 %*% eigenVals) / (mWeights1Sq %*% eigenVals);
      Samt2 <- (weightedSpectra2 %*% eigenVals) / (mWeights2Sq %*% eigenVals);
      den <- as.vector(Samt1)*as.vector(Samt2);
      mJackKnifeCoherenceEst[,k] <-  abs(num)^2 /den;
      mJackKnifePhaseEst[,k] <- Arg(num/sqrt(den));
      eigenVals[k] <- eval_k;
   }
   return(list(mJackKnifeCoherence=mJackKnifeCoherenceEst,
               mJackKnifePhase=mJackKnifePhaseEst));
}


coherence <- function(  Smt1,
                        Smt2) {
   #TODO: checks required

   num <- apply(Smt1*Conj(Smt2), 1, mean);
   return(num / (Samt1*Samt2));
}

mtmCohere2 <- function( timeSeriesA,
                        timeSeriesB,
                        k,
                        deltaT=1.0,
                        nw=4,
                        frequencies="nextPowerOf2",
                        sdfTransformation=convertTodB,
                        N=length(timeSeriesA),
                        returnJackKnifeEst=T) {
   sigma2A = sampleVariance.bias(timeSeriesA);
   sigma2B = sampleVariance.bias(timeSeriesB);
   tapers <- dpss.taper.2(N, k, nw=nw);

   resultA <- multitaperSpectralEstimate(timeSeriesA,
                     tapers$v,
                     nNonZeroFreqs=frequencies,
                     samplingTime=deltaT,
                     sdfTransformation=F,
                     returnEigenCoef=T);

   FdistrResA <- ThomsonFTest(  resultA$mEigenCoef,
                              tapers$v,
                              deltaT);

   aresultA <-  eigenSpectraToAdaptiveMultitaperSpectralEstimate(
                        resultA$mEigenSpectra,
                        tapers$eigen,
                        sigma2A,
                        sdfTransformation=F,
                        returnWeights=T);

   resultB <- multitaperSpectralEstimate(timeSeriesB,
                     tapers$v,
                     nNonZeroFreqs=frequencies,
                     samplingTime=deltaT,
                     sdfTransformation=F,
                     returnEigenCoef=T);

   FdistrResB <- ThomsonFTest(  resultA$mEigenCoef,
                              tapers$v,
                              deltaT);

   aresultB <-  eigenSpectraToAdaptiveMultitaperSpectralEstimate(
                        resultB$mEigenSpectra,
                        tapers$eigen,
                        sigma2B,
                        sdfTransformation=F,
                        returnWeights=T);

   coh <- coherenceAdaptive(aresultA$resultSDF,
                              aresultB$resultSDF,
                              resultA$mEigenCoef,
                              resultB$mEigenCoef,
                              aresultA$mWeights,
                              aresultB$mWeights,
                              tapers$eigen);
   jackKnifeCoherence <- NULL;
   jackKnifePhase <- NULL;

   if(returnJackKnifeEst) {
      jackKnifeEst <-  coherenceAdaptiveJackKnifeMatrix(
                              resultA$mEigenSpectra,
                              resultB$mEigenSpectra,
                              resultA$mEigenCoef,
                              resultB$mEigenCoef,
                              aresultA$mWeights,
                              aresultB$mWeights,
                              tapers$eigen);
      jackKnifeEst$mJackKnifeCoherence[,(k+1)] <- coh$magSq;
      jackKnifeCoherence <-
         apply(jackKnifeEst$mJackKnifeCoherence, 1, mean);
      jackKnifeEst$mJackKnifePhase[,(k+1)] <- coh$phase;
      jackKnifePhase <-
         apply(jackKnifeEst$mJackKnifePhase, 1, mean);
   }

   if(is.function(sdfTransformation)) {
      aresultA$resultSDF <- sdfTransformation(aresultA$resultSDF);
   }

   if(is.function(sdfTransformation)) {
      aresultB$resultSDF <- sdfTransformation(aresultB$resultSDF);
   }

   return(list(coherence=coh$magSq,
               phase=coh$phase,
               jackKnifeCoherence=jackKnifeCoherence,
               jackKnifePhase=jackKnifePhase,
               resultFreqs=resultA$resultFreqs,
               sdfA=aresultA$resultSDF,
               sdfB=aresultB$resultSDF,
               FTestA=FdistrResA$FTest,
               FTestB=FdistrResB$FTest,
               FdfNum=FdistrResA$numDf,
               FdfDen=FdistrResA$denDf));

}


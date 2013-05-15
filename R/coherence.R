## moved from analysis.R Not sure if thes belong in the PWLISPR package
## may move to multitaper Development


## ##coherenceAdaptiveJackKnifeMatrix <- function(mSpectra1,
## #                              mSpectra2,
## #                              mEigenCoef1,
## #                              mEigenCoef2,
## #                              mWeights1,
## #                              mWeights2,
## #                              eigenVals)

## # coherenceAdaptive <- function(Samt1,
## #                              Samt2,
## #                              mEigenCoef1,
## #                              mEigenCoef2,
## #                              mWeights1,
## #                              mWeights2,
## #                              eigenVals1,
## #                              eigenVals2)

## #coherence <- function(  Smt1,
## #                        Smt2)

## #mtmCohere2 <- function( timeSeriesA,
## #                        timeSeriesB,
## #                        k,
## #                        deltaT=1.0,
## #                        nw=4,
## #                        frequencies="nextPowerOf2",
## #                        sdfTransformation=convertTodB,
## #                        N=length(timeSeriesA),
## #                        returnJackKnifeEst=T)

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
   return(num / (Smt1*Smt2));
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


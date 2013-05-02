# The following functions are contained in multitaper.R
# These functions are loosley based on the lisp code provided with SAPA
# http://lib.stat.cmu.edu/sapaclisp/
# http://lib.stat.cmu.edu/sapaclisp/multitaper.lisp

# checkOrthonormality <- function(X)

# dpss.taper.2 <- function (n, k, nw = 4, nmax = 2^(ceiling(log(n, 2))))

# dpssTapersTriDiag <- function(   N,
#                                 numberOfTapers,
#                                 taperParameter=4,
#                                 printProgressP=F,
#                                 computeTrueEigenvaluesP=F)
                                 
# multitaperSpectralEstimate <- function( timeSeries,
#                     mDataTapers,
#                     centreData=T,
#                     nTapers=length(mDataTapers[1,]),
#                     nNonZeroFreqs="halfNextPowerOf2",
#                     samplingTime=1.0,
#                     recentreAfterTaperingP=T,
#                     restorePowerOptionP=T,
#                     returnEstFor0FreqP=F,
#                     sdfTransformation=convertTodB,
#                     returnEigenCoef=F)

# eigenspectraToMultitaperSpectralEstimate <- function (
#                        mEigenSpectra,
#                        nEigenSpectra=length(mEigenSpectra[1,]),
#                        sdfTransformation=convertTodB)

# eigenSpectraToAdaptiveMultitaperSpectralEstimate <- function(
#         mEigenSpectra,
#         eigenValues,
#         varianceOfTimeSeries,   #use sampleVariance.bias()
#         samplingTime=1.0,
#         nEigenSpectra=length(mEigenSpectra[1,]),
#         maximumNumberOfIterations=100,
#         sdfTransformation=convertTodB,
#         returnWeights=F)

# createCIforAMTsdfEst <- function( sdfDB, dfs, confidenceLevel=0.95)

################################################################################

#libraries required
library(waveslim);

################################################################################
# Requires:
source("basicStatistics.R");
source("utilities.R");
source("tapers.R");
source("nonparametric.R");
source("acvs.R");

################################################################################

checkOrthonormality <- function(X) {
   if(!is.matrix(X) && !is.array(X)) {
      return();
   }
   maxdev <- 0.0;
   for(k in 1:length(X[1,])) {
      cat(paste("k = ", k, ": sum of squares = ", t(X[,k])%*%X[,k], "\n"));
   }
   for(k in 1:(length(X[1,])-1)) {
      cat(paste("k = ", k, "\n"));
      for(j in (k+1):length(X[1,])) {
         dotprod = t(X[,k])%*%X[,j];
         if(abs(dotprod)  > maxdev) {
            maxdev <- abs(dotprod);
         }
         cat(paste("     dot product  = ", dotprod, "\n"));
      }
  }
  return(maxdev);
}

#modified dpss function
#uses inverse iteration
#http://lib.stat.cmu.edu/jcgs/bell-p-w
dpss.taper.2 <- function (n, k, nw = 4, nmax = 2^(ceiling(log(n, 2))))
{
    if (n > nmax)
        stop("length of taper is greater than nmax")
    w <- nw/n
    if (w > 0.5)
        stop("half-bandwidth parameter (w) is greater than 1/2")
    if (k <= 0)
        stop("positive dpss order (k) required")
    v <- matrix(0, nrow = nmax, ncol = (k + 1))

    storage.mode(v) <- "double"
    out <- .Fortran("dpss", nmax = as.integer(nmax), kmax = as.integer(k),
        n = as.integer(n), w = as.double(w), v = v, sig = double(k +
            1), totit = integer(1), sines = double(n), vold = double(n),
        u = double(n), scr1 = double(n), ifault = integer(1),
        PACKAGE = "waveslim")
    return(list(  v=out$v[1:n, 1:k],
                  eigen=1+out$sig[1:k],
                  iter=out$totit,
                  n=n,
                  w=w,
                  ifault=out$ifault) );
}

#format(dpssToEigenvalue(
#   c(0.016070174357819406, 0.08186425557890079, 0.21480029092555733, 0.36208333168638207, 0.40188579349948383, 0.24234798504319702, -0.059754642833675535, -0.3032839179323881, -0.3032839179323881, -0.059754642833675535, 0.24234798504319702, 0.40188579349948383, 0.36208333168638207, 0.21480029092555733, 0.08186425557890079, 0.016070174357819406),
#   4), nsmall=16);

#;;;  dpssTapersTriDiag  <- function
#;;;  computes a set of orthonormal dpss data tapers for a specified
#;;;  sample size N and taper parameter NW (the product of the duration N
#;;;  and half bandwidth W -- note that W is taken to be in standardized units
#;;;  so that 0 < W < 1/2).  In general, dpss-tapers-tri-diag should be used
#;;;  because it is generally  fast and it is generally quite accurate.

#given
#[1] N (required)
#    the sample size
#[2] numberOfTapers (required)
#    number of orthonormal dpss data tapers
#    to be computed
#[3] taperParameter (keyword; 4.0)
#    NW, the duration-half-bandwidth product
#    (must be such that 0 < NW/N < 1/2)
#[4] printProgressP (keyword; nil)
#    if t, prints a dot after each eigenvalue
#    and eigenvector has been computed
#[5] computeTrueEigenvaluesP (keyword; F)
#    if t, returns eigenvalues for eigenproblem
#    of Equation (378) of the SAPA book;
#    if nil, returns eigenvalues for tridiagonal
#    formulation
#returns
#[1] results a matrix with numberOfTapers columns and N rows
#    with orthonormal dpss's of orders 0, 1, ..., number-of-tapers - 1
#    in each column
#[2] eigenvalues an array of length number-of-tapers with
#    eigenvalues as specified by compute-true-eigenvalues-p

#Note: computes the dpss tapers using the tridiagonal
#      formulation (see Section 8.3 of the SAPA book)

dpssTapersTriDiag <- function(   N,
                                 numberOfTapers,
                                 taperParameter=4,
                                 printProgressP=F,
                                 computeTrueEigenvaluesP=F) {

   w <- taperParameter / N;
   if(0 >= w || w > 0.5) {
      return();
   }
   Nm1 <- N -1;
   Nm1o2 <- as.double(Nm1/2.0);
   ctpiW <- cos(2*pi*w);
   results <- matrix(NA, N, numberOfTapers);
   #;;; generate diagonal elements of symmetric tridiagonal system ...
   diagEl <- ctpiW * (Nm1o2- 0:Nm1)^2;
   #;;; generate off-diagonal elements ...
   offDiagEl <- 0.5 * 1:Nm1 * Nm1:1

   if(printProgressP) {
      cat("finding eigenvalues  ");
   }

   #;;; get eigenvalues ...
   eigenvalues <- largestEigenvaluesOfTridiagonalMatrix(
      diagEl,
      offDiagEl,
      numberOfTapers,
      printProgressP=printProgressP);
   if(printProgressP) {
      cat("finding eigenvectors ");
   }
   #;;; get eigenvectors (the dpss's) ...
   for( k in 1:numberOfTapers ) {
      results[,k] <- fastTridiagEigenvalueToDPSS(eigenvalues[k],
         k -1, diagEl, offDiagEl)$resultDPSS;
      if(computeTrueEigenvaluesP) {
         eigenvalues[k] <- dpssToEigenvalue( results[,k],
            taperParameter);
      }
      if(printProgressP) {
         cat(".");
      }
   }
   if(printProgressP) {
      cat("\n");
   }
   return(list(results=results, eigenvalues=eigenvalues));
}


multitaperSpectralEstimate <- function( timeSeries,
                     mDataTapers,
                     centreData=T,
                     nTapers=length(mDataTapers[1,]),
                     nNonZeroFreqs="halfNextPowerOf2",
                     samplingTime=1.0,
                     recentreAfterTaperingP=T,
                     restorePowerOptionP=T,
                     returnEstFor0FreqP=F,
                     sdfTransformation=convertTodB,
                     returnEigenCoef=F) {
   if(   !is.matrix(mDataTapers) ||
         length(timeSeries) != length(mDataTapers[,1])) {
      return();
   }
   sampleSize <- length(timeSeries);
   nFreqs <- getNFreqs(nNonZeroFreqs, sampleSize, returnEstFor0FreqP);
   nDFT <- getNDFT(nNonZeroFreqs, sampleSize);
   offSet <- if(returnEstFor0FreqP) 0 else 1;
   fiddleFactorSDF <- samplingTime / as.double(sampleSize);
   fiddleFactorFreq <- 1 / as.double(nDFT * samplingTime);
   if(centreData) {
      timeSeries <- timeSeries - mean(timeSeries);
   }
   mEigenCoef <- NULL;
   if(returnEigenCoef) {
      mEigenCoef <- matrix(0, nFreqs, nTapers);
   }
   mEigenSpectra <- matrix(0, nFreqs, nTapers);

   for (k in 1:nTapers) {
      tSDF <- taperTimeSeriesWSupplied(timeSeries,
                     mDataTapers[,k],
                     recentreAfterTaperingP = recentreAfterTaperingP,
                     restorePowerOptionP = restorePowerOptionP)$result;
      tSDF <- c(tSDF, rep(0.0, nDFT-sampleSize));
      tSDF <- fft(tSDF);
      if(returnEigenCoef) {
         mEigenCoef[,k] <-
            tSDF[(1+offSet):(nFreqs+offSet)]*sqrt(fiddleFactorSDF);
      }
      mEigenSpectra[,k] <-
         abs(tSDF[(1+offSet):(nFreqs+offSet)])^2*fiddleFactorSDF;

   }
   resultSDF = apply(mEigenSpectra, 1, mean);
   if(is.function(sdfTransformation)) {
      resultSDF <- sdfTransformation(resultSDF);
   }
   resultFreqs <- ((0+offSet):(nFreqs+offSet-1))*fiddleFactorFreq;
   return(list(   resultSDF=resultSDF,
                  resultFreqs=resultFreqs,
                  nFreqs=nFreqs,
                  mEigenCoef=mEigenCoef,
                  mEigenSpectra=mEigenSpectra));
}

#Multitaper test using Percival Sum test from page 325 6.6c
#this is also used to test the lisp code provided with
#SAPA.
#for (N in c(22, 32, 63, 64, 65)) {
#   timeSeries <- rnorm(N,mean=100);
#   #note Percival do not calculate the sample variance
#   #correction factor (n-1)/n
#   sigma2 <- sampleVariance.bias(timeSeries);
#   samplingTime <- .25;
#   #here is where we place the dpss
#   theMTMSpec <- multitaperSpectralEstimate(timeSeries,
#                     dpss.taper.2(N, 4)$v,
#                     nNonZeroFreqs="Fourier",
#                     sdfTransformation=F,
#                     samplingTime=samplingTime);
#   Nminusf <- length(theMTMSpec$resultSDF);
#   ParsevalSum <- 0;
#   if((N%%2) ==0 ) {
#      ParsevalSum <- ((2* sum(theMTMSpec$resultSDF[1:(Nminusf-1)]) +
#                              theMTMSpec$resultSDF[Nminusf]) /
#                                 (N * samplingTime))
#   }
#   else {
#      ParsevalSum <- ((2* sum(theMTMSpec$resultSDF)) / (N * samplingTime));
#   }
#   cat(paste(  " n = ",
#         format(N, width=2),
#         ", N-f = ",
#         format(Nminusf, width=2) ,
#         ": ", round(sigma2, digits=12),
#         " ",
#         round(ParsevalSum, digits=12),
#         " ",
#         format(ParsevalSum/sigma2, nsmall=12), "\n"));
#}

eigenspectraToMultitaperSpectralEstimate <- function (
                        mEigenSpectra,
                        nEigenSpectra=length(mEigenSpectra[1,]),
                        sdfTransformation=convertTodB) {
   resultSDF = apply(mEigenSpectra, 1, mean);
   if(is.function(sdfTransformation)) {
      resultSDF <- sdfTransformation(resultSDF);
   }
   return(resultSDF);
}

# Here we create two multitaper spectral estimates, one using 4 tapers
# and the other using 3 tapers.  We check the Parseval result in
# both cases.
#testMTM <- function() {
#   timeSeries <- rnorm(32,100);
#   sigma2 <- sampleVariance.bias(timeSeries);
#   samplingTime <- .25;
#   result <- multitaperSpectralEstimate(timeSeries,
#                        dpss.taper.2(32, 4)$v,
#                        nNonZeroFreqs="Fourier",
#                        sdfTransformation=F,
#                        samplingTime=samplingTime);
#   mt4 <- result$resultSDF;
#   nMinusf <- result$nFreqs;
#   mEigenSpectra <- result$mEigenSpectra;
#   cat(paste(  format(((2* sum(mt4[1:(nMinusf-1)]) +
#                     mt4[nMinusf]) /
#                     (32 * samplingTime * sigma2)), nsmall=12),
#               "\n"));
#
#   mt3 <-  eigenspectraToMultitaperSpectralEstimate(
#               mEigenSpectra,
#               nEigenSpectra=3,
#               sdfTransformation=F);
#   cat(paste(  format(((2* sum(mt3[1:(nMinusf-1)]) +
#                     mt3[nMinusf]) /
#                     (32 * samplingTime * sigma2)), nsmall=12),
#               "\n"));
#   return();
#}
#testMTM();

eigenSpectraToAdaptiveMultitaperSpectralEstimate <- function(
         mEigenSpectra,
         eigenValues,
         varianceOfTimeSeries,   #use sampleVariance.bias()
         samplingTime=1.0,
         nEigenSpectra=length(mEigenSpectra[1,]),
         maximumNumberOfIterations=100,
         sdfTransformation=convertTodB,
         returnWeights=F) {
   if( !is.matrix(mEigenSpectra) ||
         length(mEigenSpectra[1,]) <= 1) {
      return();
   }

   mWeights <- NULL;
   Nminusf <- length(mEigenSpectra[,1]);
   resultSDF <- array(0.0, Nminusf);
   sig2timesDeltaT = varianceOfTimeSeries*samplingTime;
   lambda0 <- eigenValues[1];
   lambda1 <- eigenValues[2];
   lamdba0plusLambda1 = lambda0+lambda1;
   firstEigenSpectrum <- mEigenSpectra[,1];
   secondEigenSpectrum <- mEigenSpectra[,2];
   eWeights <- array(0.0, nEigenSpectra);
   newEstBot <- 0;
   newEstTop <- 0;
   maxIterations <-0;
   resultDF <- array(0.0, Nminusf);
   if(returnWeights) {
      mWeights <- matrix(0.0, Nminusf,nEigenSpectra);
   }

   for( i in 1:Nminusf) {
      previousEst <- (  ((lambda0 * firstEigenSpectrum[i]) +
                          (lambda1 * secondEigenSpectrum[i])) /
                        lamdba0plusLambda1);
      localMaxIterations <- 1;
      for( j in 1:maximumNumberOfIterations) {
         newEstTop <- 0.0;
         newEstBot <- 0.0;
         for( k in 1:nEigenSpectra) {
            eWeights[k] <-  previousEst /
               ((eigenValues[k]* previousEst)
               + ( (1.0 - eigenValues[k]) * sig2timesDeltaT));
            newEstTop <- newEstTop + (eWeights[k]^2
               * eigenValues[k]
               * mEigenSpectra[i,k]);
            newEstBot <- newEstBot + (eWeights[k]^2 * eigenValues[k]);
         }

         resultSDF[i] <- newEstTop/newEstBot;
         if((abs(resultSDF[i] - previousEst) / previousEst) < 0.05)  {
            break;
         }

         previousEst <- resultSDF[i];
         localMaxIterations <- localMaxIterations +1;

      }
      if(localMaxIterations > maxIterations) {
         maxIterations <- localMaxIterations;
      }

      theSum <- 0.0;
      for( k in 1:nEigenSpectra) {
        theSum <- theSum + ( eWeights[k]^4 * eigenValues[k]^2);
      }
      resultDF[i] <- (2.0 * newEstBot^2) / theSum;
      if(returnWeights) {
         mWeights[i,] <- eWeights;
      }
   }

   if(is.function(sdfTransformation)) {
      resultSDF <- sdfTransformation(resultSDF);
   }
   return(list(   resultSDF=resultSDF,
                  resultDF=resultDF,
                  maxIterations=maxIterations,
                  mWeights=mWeights));
}

## testAMTM <- function(returnWeights=F) {
##   a20PtTS <- c(71.0, 63.0, 70.0, 88.0, 99.0, 90.0,110.0, 135.0, 128.0, 154.0, 156.0, 141.0, 131.0, 132.0, 141.0, 104.0, 136.0, 146.0, 124.0, 129.0);
##   nTS <- length(a20PtTS);
##   sigma2 <- sampleVariance.bias(a20PtTS);
##   samplingTime <- 0.25;
##   dpssNW4 <- dpss.taper.2(nTS, 7);
##   result <- multitaperSpectralEstimate(a20PtTS,
##                     dpssNW4$v,
##                     nNonZeroFreqs="Fourier",
##                     sdfTransformation=F,
##                     samplingTime=samplingTime);
##   mt4 <- result$resultSDF;
##   nMinusf <- result$nFreqs;
##   freqs <- result$resultFreqs;
##   eigenValues <- result$eigen;
##   mEigenSpectra <- result$mEigenSpectra;
##   cat(paste(  format(((2* sum(mt4[1:(nMinusf-1)]) +
##                     mt4[nMinusf]) /
##                     (nTS * samplingTime * sigma2)), nsmall=12),
##                     "\n"));

##   aResult <-  eigenSpectraToAdaptiveMultitaperSpectralEstimateF(
##                        mEigenSpectra,
##                        dpssNW4$eigen,
##                        sigma2,
##                        sdfTransformation=F,
##                        returnWeights=returnWeights);

##   amt4 <- aResult$resultSDF;
##   dfs <- aResult$resultDF;
##   cat(paste(  format(((2* sum(amt4[1:(nMinusf-1)]) +
##                     amt4[nMinusf]) /
##                     (nTS * samplingTime * sigma2)), nsmall=12),
##                     "\n"));

##   amt4 <- convertTodB(amt4);
##   cat(paste(  " ", format(freqs, nsmall=4), " ",
##               format(amt4, nsmall=4) , " ",
##               format(dfs, nsmall=4), "\n"), fill=1);
##   return(aResult$mWeights);
## }
## testAMTM();

createCIforAMTsdfEst <- function( sdfDB, dfs, confidenceLevel=0.95) {
   p <- (1- confidenceLevel)/2;
   oneMinusP <- 1-p;
   Nminusf <- length(sdfDB);
   upperCI <-  sdfDB;
   lowerCI <- sdfDB;
   upperCI <- upperCI + convertTodB(dfs / qchisq(p, dfs));
   lowerCI <- lowerCI + convertTodB(dfs / qchisq(oneMinusP, dfs));
   return(list(upperCI=upperCI, lowerCI=lowerCI));

}

#testCI <- function() {
#   a20PtTS <- c(71.0, 63.0, 70.0, 88.0, 99.0, 90.0,110.0, 135.0, 128.0, 154.0, 156.0, 141.0, 131.0, 132.0, 141.0, 104.0, 136.0, 146.0, 124.0, 129.0);
#   nTS <- length(a20PtTS);
#   sigma2 <- sampleVariance.bias(a20PtTS);
#   samplingTime <- 0.25;
#   dpssNW4 <- dpss.taper.2(nTS, 7);
#   result <- multitaperSpectralEstimate(a20PtTS,
#                     dpssNW4$v,
#                     nNonZeroFreqs="Fourier",
#                     samplingTime=samplingTime);
#   mt4 <- result$resultSDF;
#   nMinusf <- result$nFreqs;
#   freqs <- result$resultFreqs;
#   eigenValues <- result$eigen;
#   mEigenSpectra <- result$mEigenSpectra;
#   aResult <-  eigenSpectraToAdaptiveMultitaperSpectralEstimate(
#                        mEigenSpectra,
#                        dpssNW4$eigen,
#                        sigma2,
#                        sdfTransformation=F);
#
#   amt4 <- aResult$resultSDF;
#   dfs <- aResult$resultDF;
#   amt4 <- convertTodB(amt4);
#   cis <- createCIforAMTsdfEst(amt4, dfs);
#   cat(paste(format(freqs, nsmall=4), " ",
#            format(mt4, nsmall=4) , "     ",
#            format(cis$lowerCI, nsmall=4), " ",
#            format(amt4, nsmall=4), " ",
#            format(cis$upperCI, nsmall=4), "\n"), fill=1);
#   return();
#}
#testCI();

sturmSequenceCount <- function(testLambda,
                        diagEl,
                        offDiagEl,
                        squaredOffDiagEl,
                        machep) {
   endVal <- length(diagEl);
   count <- 0;
   bottom <- 1.0;
   for( i in 0:(endVal-1) ) {
      ratio <- NULL;
      if(i == 0) {
         ratio <- 0.0;
      } else if(bottom == 0) {
         if(squaredOffDiagEl[i] == 0) {
            ratio <- 0.0;
         } else {
            ratio <- abs(offDiagEl[i]) / machep;
         }
      } else {
         ratio <- squaredOffDiagEl[i] / bottom;
      }
      bottom <- diagEl[i+1] - testLambda - ratio;
      if(bottom < 0) {
         count <- count +1;
      }
   }
   return(count);
}

##########################################################################
#helper functions below...

largestEigenvaluesOfTridiagonalMatrix <- function(diagEl,
                  offDiagEl,
                  numberOfTapers,
                  squaredOffDiagEl=(offDiagEl)^2,
                  macheps=.Machine$double.eps,
                  printProgressP=F) {
   n <- length(diagEl);
   nm1 <- n -1;
   #;;; Set to zero all elements of squared-off-diag that correspond to
   #;;; small elements of off-diag (see do-loop 40 of tridib) ...
   previousAbs <- diagEl[1];
   for( i in 1:nm1) {
      if(offDiagEl[i] <=  (macheps * (previousAbs + diagEl[i +1]))) {
         squaredOffDiagEl[i] <- 0.0;
      }
      previousAbs <- diagEl[i +1];
   }
   #;;; Use Equation (6) of Barth, Martin and Wilkinson to find
   #;;; upper and lower bounds for all eigenvalues ..
   absOffDiagElBehind <- abs(offDiagEl[1]);
   absOffDiagElAhead <- abs(offDiagEl[nm1]);
   sumOfAbs <- NULL;
   lowerBoundAllEigenvalues <- min(diagEl[1] - absOffDiagElBehind,
      diagEl[n] - absOffDiagElAhead);
   upperBoundAllEigenvalues <- max(diagEl[1] + absOffDiagElBehind,
      diagEl[n] + absOffDiagElAhead);
   for( i in 2:(nm1) ) {
      absOffDiagElAhead <- abs(offDiagEl[i]);
      sumOfAbs <- absOffDiagElBehind + absOffDiagElAhead;
      if( lowerBoundAllEigenvalues > (diagEl[i] - sumOfAbs) ) {
         lowerBoundAllEigenvalues <- diagEl[i] - sumOfAbs;
      }
      if( upperBoundAllEigenvalues <  (diagEl[i] + sumOfAbs) ) {
         upperBoundAllEigenvalues <- diagEl[i] + sumOfAbs;
      }
      absOffDiagElBehind <- absOffDiagElAhead;
   }

   #;;; Expand upper and lower bounds a little (evidently to avoid
   #;;; numerical problems -- see code following do-loop 40 of tridib) ...
   eigenvalues <- array(NA, numberOfTapers);
   upperBounds <- array(upperBoundAllEigenvalues, numberOfTapers);
   lowerBounds <- array(lowerBoundAllEigenvalues, numberOfTapers);
   for(k in 1:numberOfTapers ) {
      #;;; use bisection to isolate eigenvalues, keeping track
      #;;; of bisections at each k to speed up subsequent k's
      updateOtherBounds <- if(k < numberOfTapers) T else F;
      upperTargetCount <- n - k +1;
      lowerTargetCount <- upperTargetCount -1;
      currentCount <- NULL;
      U <- upperBounds[k];
      L <- lowerBounds[k];
      midPoint <- (U + L)/2.0;
      while(abs(U - L) > (macheps * (abs(U) + abs(L)))) {
         currentCount <- sturmSequenceCount(midPoint,
                  diagEl,
                  offDiagEl,
                  squaredOffDiagEl,
                  macheps);
         if(updateOtherBounds) {
            N_m_currentCount <- n - currentCount;
            jToInc <- N_m_currentCount +1;
            for(j in 1:min(N_m_currentCount, numberOfTapers) ) {
               if(midPoint > lowerBounds[j]) {
                  lowerBounds[j] <- midPoint;
               }
            }
            nTapersMN_m_current <- numberOfTapers - N_m_currentCount;
            if(nTapersMN_m_current > 0) {
               for(j in 0:(nTapersMN_m_current -1) ) {
                  if(midPoint < upperBounds[jToInc]) {
                     upperBounds[jToInc] <- midPoint;
                  }
                  jToInc <- jToInc +1;
               }
            }
            if(currentCount >= lowerTargetCount) {
               updateOtherBounds <- F;
            }
         }
         if(currentCount <= lowerTargetCount) {
            L <- midPoint;
         } else {
            U <- midPoint;
         }
         midPoint <- (U + L)/2.0;
      }
      eigenvalues[k] <- midPoint;
      if(printProgressP) {
         cat(".");
      }
   }
   if(printProgressP) {
      cat("\n");
   }
   return(eigenvalues);
}

#a <- largestEigenvaluesOfTridiagonalMatrix(
#   c(3.444205370125435E-15, 2.586980922449772E-15, 1.852217110156346E-15, 1.2399139332451574E-15, 7.500713917162062E-16, 3.8268948556949295E-16, 1.3776821480501744E-16, 1.5307579422779716E-17, 1.5307579422779716E-17, 1.3776821480501744E-16, 3.8268948556949295E-16, 7.500713917162062E-16, 1.2399139332451574E-15, 1.852217110156346E-15, 2.586980922449772E-15, 3.444205370125435E-15),
#   c(7.5, 14.0, 19.5, 24.0, 27.5, 30.0, 31.5, 32.0, 31.5, 30.0, 27.5, 24.0, 19.5, 14.0, 7.5),
#   4,
#   printProgressP=F)
#
#print(a);

#largestEigenvaluesOfTridiagonalMatrix(
#   c(3.444205370125435E-15, 2.586980922449772E-15, 1.852217110156346E-15, 1.2399139332451574E-15, 7.500713917162062E-16, 3.8268948556949295E-16, 1.3776821480501744E-16, 1.5307579422779716E-17, 1.5307579422779716E-17, 1.3776821480501744E-16, 3.8268948556949295E-16, 7.500713917162062E-16, 1.2399139332451574E-15, 1.852217110156346E-15, 2.586980922449772E-15, 3.444205370125435E-15),
#   c(7.5, 14.0, 19.5, 24.0, 27.5, 30.0, 31.5, 32.0, 31.5, 30.0, 27.5, 24.0, 19.5, 14.0, 7.5),
#   4,
#   printProgressP=F)

generateInitialGuessAtDPSS <- function(
                           N,
                           orderOfTaper,
                           halfSizeP=F) {
   result <- NULL;
   if(halfSizeP) {
      result <- 1:(N/2);
   } else {
      result <- 1:N;
   }
   result <- result / (N +1);
   product <- array(1.0, length(result));
   orderPlus1 <- orderOfTaper +1;
   if(orderOfTaper != 0) {
      for(j in 1:orderOfTaper  ) {
         product <- product * ((j / orderPlus1) - result);
      }
   }
   result <- result^2 * (result -1)^2 * product;
   result <- result / sqrt(sum(result^2));
   return(result);
}
#generateInitialGuessAtDPSS(16, 4)

symmetricTridiagonalInfinityNorm <- function(diagEl, offDiagEl) {
   n <- length(diagEl);
   nm1 <- n -1;
   nm2 <- n -2;

   infinityNorm <- max(abs(diagEl[1])+abs(offDiagEl[1]),
      abs(offDiagEl[1:nm2])+abs(offDiagEl[2:nm1])+abs(diagEl[2:nm1]),
      abs(diagEl[n])+abs(offDiagEl[nm1]));
   return(infinityNorm);
}

#symmetricTridiagonalInfinityNorm(
#   c(-58.53461317004094, -58.53461317004094, -58.53461317004094, -58.53461317004094, -58.53461317004094, -58.53461317004094, -58.53461317004094, -26.53461317004094),
#   c(7.5, 14.0, 19.5, 24.0, 27.5, 30.0, 31.5, 32.0));

#symmetricTridiagonalInfinityNorm(
#   c(-38.20630088826668, -38.20630088826668, -38.20630088826668, -38.20630088826668, -38.20630088826668, -38.20630088826668, -38.20630088826668, -6.206300888266682),
#   c(7.5, 14.0, 19.5, 24.0, 27.5, 30.0, 31.5, 32.0))


# symmetricTridiagonalSolve <- function
#  given
#   [1] DiagEl (required)
#       <=> diagonal part of symmetric tridiagonal matrix A;
#           trash on output
#   [1] ofDiagEl (required)
#       <=> off-diagonal part of symmetric tridiagonal matrix;
#           trash on output
#   [2] b (required)
#       <=> on input, the right-hand side vector;
#returns
#   [1] the solution to A X = b,
#---
#Note: this is an implementation of Algorithm 4.3.6,
#p. 156, Golub and Van Loan, 1989, with modifications
#to avoid divides by zero"
#used in testing

symmetricTridiagonalSolve <- function(diagEl, offDiagEl, b) {
   temp <- NULL;
   n <- length(diagEl);
   nm1 <- n -1;
   zeroReplacement <- .Machine$double.eps *
      symmetricTridiagonalInfinityNorm(diagEl, offDiagEl);

   if(abs(diagEl[1]) < .Machine$double.eps) {
      diagEl[1] <- zeroReplacement;
   }

   for(k in 1:nm1 ) {
      kp1 <- k +1;
      tempVal <- offDiagEl[k];
      offDiagEl[k] <- tempVal / diagEl[k];
      diagEl[kp1] <- diagEl[kp1] - tempVal * offDiagEl[k];
      if( abs(diagEl[kp1]) < .Machine$double.eps) {
         diagEl[kp1] <- zeroReplacement;
      }
   }
   for( k in 1:nm1 ) {
      kp1 <- k +1;
      b[kp1] <- b[kp1] - b[k] * offDiagEl[k];
   }
   b[n] <- b[n]/diagEl[n];

   for(k in nm1:1 ) {
      b[k] <-
         b[k]/diagEl[k] - offDiagEl[k]* b[k +1];
   }
   return(b);
}

#symmetricTridiagonalSolve(
#   c(-58.53461317004094, -58.53461317004094, -58.53461317004094, -58.53461317004094, -58.53461317004094, -58.53461317004094, -58.53461317004094, -26.53461317004094),
#   c(7.5, 14.0, 19.5, 24.0, 27.5, 30.0, 31.5, 32.0),
#   c(0.026387905030664526, 0.09276997862342998, 0.1818291581019227, 0.278722246886394, 0.37107991449371985, 0.4490066965374011, 0.5050809947275632, 0.5343550768709567))
#
#


fastTridiagEigenvalueToDPSS <- function(
   eigenvalue,
   orderOfTaper,
   diagEl,
   offDiagEl,
   eps=10*.Machine$double.eps,
   maximumNumberOfIterations=25,
   b=generateInitialGuessAtDPSS(length(diagEl),
      orderOfTaper,
      halfSizeP=T)) {
   N <- length(diagEl);
   M <- if(N%%2 ==0) N/2 else (N-1)/2;

   shorterDiag <- diagEl[1:M];
   shorterOffDiag <- offDiagEl[1:M];
   #;;; even N
   if(N%%2 ==0) {
      if(orderOfTaper%%2 == 0) {
         shorterDiag[M] <- shorterDiag[M] + offDiagEl[M];
      } else {
         shorterDiag[M] <- shorterDiag[M] - offDiagEl[M];
      }
   #;;; odd N
   } else {
      if(orderOfTaper%%2 == 0) {
         shorterDiag[M] <- shorterDiag[M] +
            2.0 * offDiag[M]^2 / eigenvalue;
      }
   }
   shorterDiag <- shorterDiag - eigenvalue;
   bOld <- b;
   endSum <- max(2, trunc(N/(2 * (orderOfTaper +1))));
   iterTotal <- NULL;
   for( i in 1:maximumNumberOfIterations ) {
      b <- symmetricTridiagonalSolve(shorterDiag, shorterOffDiag, b);
      bSum <- sum(b[1:endSum]);
      tempFactor <- if(bSum > 0) 1.0 else -1.0;
      b <-  b * tempFactor / sqrt(sum(b^2));
      avgAbsDiff <- compareSeqs(b, bOld)$avgAbsDiff;
      if(avgAbsDiff < eps) {
         iterTotal <- i;
         break;
      } else {
         bOld <- b;
      }
   }
   resultDPSS <- array(NA, N);
   resultDPSS[1:M] <- b[1:M];
   if(N%%2 == 0) {
      if(orderOfTaper%%2 == 0) {
         resultDPSS[N:(N-M+1)] <- resultDPSS[1:M];
      } else {
         resultDPSS[N:(N-M+1)] <- -1 * resultDPSS[1:M];
      }
   } else if(orderOfTaper%%2 ==0) {
      resultDPSS[N:(N-M+1)] <- resultDPSS[1:M];
      resultDPSS[M +1] <- 2 * resultDPSS[M] *
         offDiagEl[M] / eigenvalue;
   } else { #;;; N is odd, and order of taper is odd
      resultDPSS[N:(N-M+1)] <- -1 * resultDPSS[1:M];
      resultDPSS[M +1] <- 0.0;
   }
   resultDPSS <- resultDPSS / sqrt(sum(resultDPSS^2));
   return(list(resultDPSS=resultDPSS, iterTotal=iterTotal));
}

#fastTridiagEigenvalueToDPSS(
#   58.53461317004094,
#   0,
#   c(3.444205370125435E-15, 2.586980922449772E-15, 1.852217110156346E-15, 1.2399139332451574E-15, 7.500713917162062E-16, 3.8268948556949295E-16, 1.3776821480501744E-16, 1.5307579422779716E-17, 1.5307579422779716E-17, 1.3776821480501744E-16, 3.8268948556949295E-16, 7.500713917162062E-16, 1.2399139332451574E-15, 1.852217110156346E-15, 2.586980922449772E-15, 3.444205370125435E-15),
#   c(7.5, 14.0, 19.5, 24.0, 27.5, 30.0, 31.5, 32.0, 31.5, 30.0, 27.5, 24.0, 19.5, 14.0, 7.5))
#

#fastTridiagEigenvalueToDPSS(
#   29.218607313170303,
#   3,
#   c(3.444205370125435E-15, 2.586980922449772E-15, 1.852217110156346E-15, 1.2399139332451574E-15, 7.500713917162062E-16, 3.8268948556949295E-16, 1.3776821480501744E-16, 1.5307579422779716E-17, 1.5307579422779716E-17, 1.3776821480501744E-16, 3.8268948556949295E-16, 7.500713917162062E-16, 1.2399139332451574E-15, 1.852217110156346E-15, 2.586980922449772E-15, 3.444205370125435E-15),
#   c(7.5, 14.0, 19.5, 24.0, 27.5, 30.0, 31.5, 32.0, 31.5, 30.0, 27.5, 24.0, 19.5, 14.0, 7.5))


#given dpss (a vector of length N) and NW,
#computes the corresponding eigenvalue using
#the method of Exercise [8.1], page 390,
#of the SAPA book
dpssToEigenvalue <- function(dpss, NW) {
   dpssACVS <- acvs(dpss, centreDataP=F)$acvs;
   N <- length(dpss);
   W <- NW/N;

   jVector <- 1:(N-1);
   vectorOfRatios <- c(2*W, sin(2*pi*W*jVector) / (pi*jVector));

   dpssACVS <-  dpssACVS * N;

   #;;; Note: both vector of ratios and dpss-acvs
   #;;;       roughy decrease in amplitude with increasing index,
   #;;;       so we sum things up in reverse order.
   eigenvalue <- sum(dpssACVS[N:2] * vectorOfRatios[N:2]);

   return(2*eigenvalue + 2*W*dpssACVS[1]);
}

#format(dpssToEigenvalue(
#   c(8.394429999264533E-4, 0.006551529503865821, 0.026942530191531977, 0.07617175228355653, 0.16388769648613713, 0.2823629404487185, 0.4007031279784092, 0.47568617075492176, 0.47568617075492176, 0.4007031279784092, 0.2823629404487185, 0.16388769648613713, 0.07617175228355653, 0.026942530191531977, 0.006551529503865821, 8.394429999264533E-4),
#   4), nsmall=16)

#dpssToEigenvalue(
#   c(8.394429999264533E-4, 0.006551529503865821, 0.026942530191531977, 0.07617175228355653, 0.16388769648613713, 0.2823629404487185, 0.4007031279784092, 0.47568617075492176, 0.47568617075492176, 0.4007031279784092, 0.2823629404487185, 0.16388769648613713, 0.07617175228355653, 0.026942530191531977, 0.006551529503865821, 8.394429999264533E-4),
#   4)
#

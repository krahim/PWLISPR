## # The following functions are contained in dft-and-fft.R
## # These functions are loosley based on the lisp code provided with SAPA
## # http://lib.stat.cmu.edu/sapaclisp/
## # http://lib.stat.cmu.edu/sapaclisp/dft-and-fft.lisp

## #note that the rest of the code uses the fft function provided by R

## #;;; The functions dft and inverse-dft implement the discrete Fourier
## #;;; transform and its inverse using, respectively, Equations (110a)
## #;;; and (111a) of the SAPA book.  They work on vectors of real
## #;;; or complex-valued numbers.  The vectors can be of any length.
## #;;; When the length of the input vector is a power of 2, these functions
## #;;; use fft_, a fast Fourier transform algorithm based upon one in Marple's
## #;;; 1987 book; when the length is not a power of 2, they use a chirp
## #;;; transform algorithm (as described in Section 3.10 of the SAPA book).

# dft <- function(X, samplingTime =1.0)
# inverseDFT <- function(X, samplingTime =1.0)

# fft_ <- function(complexVector)
# inverseFFT_ <- function(complexVector) 
# dftChirp_ <- function(complexVector)

################################################################################

## # Requires:
## ##source("utilities.R");

## ################################################################################

## #note that the rest of the code uses the fft function provided by R
## #lisp versions,  




dft <- function(X, samplingTime =1.0) {
    ##given:
    ##[1] x (required)
    ## <=> a vector of real or complex-valued numbers
    ##[2] sampling-time (keyword; 1.0)
    ##returns:
    ##[1] the discrete Fourier transform of x, namely,
    ##                      N-1
                                        #X(n) =  sampling-time SUM x(t) exp(-i 2 pi n t/N)
    ##                      t=0
    ##---
    ##Note: see Equation (110a) of the SAPA book
    
   if(is.numeric(powerOf2(length(X)))) {
      X <- fft_(X);
   } else {
      X <- dftChirp_(X);
   }
   if(samplingTime != 1.0) {
      return(X*samplingTime);
   }
   return(X);
}

inverseDFT <- function(X, samplingTime =1.0) {
    ## #given:
    ## #[1] X (required)
    ## # <=> a vector of real or complex-valued numbers
    ## #[2] sampling-time (keyword; 1.0)
    ## #returns:
    ## #[1] X, with contents replaced by
    ## # the inverse discrete Fourier transform of X, namely,
    ## #                                N-1
    ## # x(t) =  (N sampling-time)^{-1} SUM X(n) exp(i 2 pi n t/N)
    ## #                                n=0
    ## #---
    ## #Note: see Equation (111a) of the SAPA book
    N <- length(X);
    oneOverNXSamplingTime <- 1/(N*samplingTime);
    X <- Conj(X);
    if(is.numeric(powerOf2(N))) {
        X <- fft_(X);
    } else {
        X <- dftChirp_(X);
    }
    return(oneOverNXSamplingTime*Conj(X));
}

#inverseDFTEx(dftEx(c(1.0, 6.2, pi, -7.0, - pi)))
#[1]  1.000000-0i  6.200000-0i  3.141593+0i -7.000000+0i -3.141593+0i

#;;;  The following code is based on a Lisp version (with modifications) of
#;;;  the FFT routines on pages 54--6 of ``Digital Spectral Analysis with
#;;;  Applications'' by Marple, 1987.  The function fft! only works
#;;;  for sample sizes which are powers of 2.
#;;;  If the vector used as input to fft_  has elements x(0), ..., x(N-1),
#;;;  the output vector has elements X(0), ..., X(N-1),
#;;;  where
#;;;
#;;;                 N-1
#;;;         X(k) =  SUM x(n) exp(-i 2 pi n k/N)
#;;;                 n=0
#;;;
#given:
#[1] complexVector (required)
#    <=> a vector of real or complex-valued numbers
#computes the discrete Fourier transform of complex-vector
#using a fast Fourier transform algorithm and
#returns:
#[1] complexVector, with the discrete Fourier transform of the 
#input, namely,
#             N-1
#    X(n) =   SUM x(t) exp(-i 2 pi n t/N)
#             t=0
#---
#Note: see Equation (110a) of the SAPA book
#   with the sampling time set to unity
#Note this is slower than calling fft in R
fft_ <- function(complexVector) {
   n <- length(complexVector);
   exponent <- powerOf2(n);
   W <- preFFT(n);
   MM <- 1;
   LL <- n;
   NN <- NULL;
   JJ <- NULL;
   c1 <- NULL;
   c2 <- NULL;

   for(k in 0:(exponent-1)) {
      NN <- LL/2;
      JJ <- MM +1;
      i <- 1;
      kk <- NN +1;
      #chirp error correctio n to n+1....
      while(i < (n +1)) {
         c1 <- complexVector[i] + complexVector[kk];
         complexVector[kk] <- complexVector[i] - complexVector[kk];
         complexVector[i] <- c1;
         i <- i + LL;
         kk <- i + NN;
      }
      if(NN > 1) {
         j <- 1;
         while(j < NN) {
            c2 <- W[JJ];
            i <- j + 1;
            kk <- j + NN +1;
            #Chirp error correction change to n +1
            while(i < (n +1)) {
               c1 <- complexVector[i] + complexVector[kk];
               complexVector[kk] <-
                  (complexVector[i] - complexVector[kk]) * c2;
               complexVector[i] <- c1;
               i <- i + LL;
               kk <- i + NN;
            }
            JJ <- JJ + MM;
            j <- j+1;
         }
         LL <- NN;
         MM <- MM * 2;
      }
   }
   j <- 1;
   nv2 <- n/2;
   k <- NULL;
   for( i in 1:(n-1)) {
      if(i < j) {
         c1 <- complexVector[j];
         complexVector[j] <- complexVector[i];
         complexVector[i] <- c1;
      }
      k <- nv2;
      repeat {
         #Chirp error correction change to j-1
         if(k > (j -1)) {
            break;
         }
         j <- j - k;
         k <- k / 2;
      }
      j <- j + k;
   }
   return(complexVector);
}

#given:
#[1] complexVector (required)
#    <=> a vector of real or complex-valued numbers
#computes the inverse discrete Fourier transform of complex-vector
#using a fast Fourier transform algorithm and
#returns:
#[1] complexVector, with 
#    the inverse discrete Fourier transform of
#    the input X(n), namely,
#                N-1
#    x(t) =  1/N SUM X(n) exp(i 2 pi n t/N)
#                n=0
#---
#Note: see Equation (111a) of the SAPA book
#   with the sampling time set to unity
inverseFFT_ <- function(complexVector) {
   N <- length(complexVector);
   complexVector <- Conj(complexVector);
   complexVector <- fft_(complexVector);
   complexVector <- Conj(complexVector)/N;
   return(complexVector);
}

#given:
#[1] complexVector (required)
#<=> a vector of real or complex-valued numbers
#computes the discrete Fourier transform of complex-vector
#using a chirp transform algorithm and
#returns:
#[1] complexVector, with 
#the discrete Fourier transform of the input, namely,
#         N-1
#X(n) =   SUM x(t) exp(-i 2 pi n t/N)
#         t=0
#---
#Note: see Equation (110a) of the SAPA book
#with the sampling time set to unity and
#also Section 3.10 for a discussion on
#the chirp transform algorithm
dftChirp_ <- function(complexVector){
   N <- length(complexVector);
   NPot <- nextPowerOf2(2*N -1);
   chirpData <- array(0, NPot);
   chirpFactors <- array(0, NPot);
   Nm1 <- N -1;
   piOverN <- pi/N;
   #;;; copy data into chirp-data, where it will get multiplied
   #;;; by the approriate chirp factors ...
   chirpData[1:N] <- complexVector;
   chirpFactors[1] <- 1;
   complexVector[1] <- 1;
   #;;; note that, at the end of this dotimes form, complex-vector
   #;;; contains the complex conjugate of the chirp factors
   NPotMjP2 <- NPot; 
   for(j in 2:N) {
      chirpFactors[NPotMjP2] <- 
         exp( complex(real=0, imaginary=piOverN*(j-1)^2));
      chirpFactors[j] <-  chirpFactors[NPotMjP2];
      complexVector[j] <- Conj(chirpFactors[j]);
      chirpData[j] <- chirpData[j] * complexVector[j];
      NPotMjP2 <- NPotMjP2 -1;
   }
   chirpData <- fft_(chirpData);
   chirpFactors <- fft_(chirpFactors);
   chirpData <- chirpData *  chirpFactors;
   chirpData <- inverseFFT_(chirpData);
   complexVector <- complexVector * chirpData[1:N];
   return(complexVector);    
}

#helper function -----------------------------------------------------------
preFFT <- function(n, complexExpVector=NULL) {
   if(is.null(powerOf2(n))) {
      return();
   }
   if(!is.null(complexExpVector)) {
      return(complexExpVector);
   }
   s <- 2*pi/n;
   c1 <- complex(real=cos(s), imaginary=-sin(s));
   c2 <- 1.0+0i;

   result <- array(NA,n);
   for(i in 1:n) {
      result[i] <- c2;
      c2 <- c2 * c1;
   }
   return(result);
}

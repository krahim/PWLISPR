##source("~/PWLisp/acvs.R")

##  reflectionCoeffsToARCoeffs <- function(reflectionCoeffs,
##  yuleWalkerAlgorithmGivenData <- function(timeSeries, p,
## reflectionCoeffsToARCoeffs <- function(reflectionCoeffs,
## yuleWalkerAlgorithmGivenACVS <- function(acvs, p) {
## arCoeffsToPrewhitenFilter <- function(arCoeffs) {
## arCoeffsToSdf <- function(arCoeffs,


reflectionCoeffsToARCoeffs <- function(reflectionCoeffs,
                                       p=length(reflectionCoeffs)) {
    ##       "given
    ##    [1] reflection-coeffs (required)
    ##        ==> vector of length p with real-valued 
    ##            reflection coefficients
    ##            phi_{1,1}, phi_{2,2}, ... , phi_{p,p}
    ##    [2] p (keyword; length of reflection-coeffs)
    ##        ==> AR model order
    ##
    ##   results (keyword; vector of length p)
    ##        <== vector of size p into which are placed
    ##            the AR coefficients phi_{j,p} for model
    ##            x_t = phi_{1,p} * x_{t-1} + phi_{2,p} * x_{t-2}
    ##                  + ... + phi_{p,p} * x_{t-p} + e_t
    ## returns
    ##    [1] results, i.e., the AR coefficients
    ## ---
    ##   Note: see Section 9.4 of the SAPA book"
    scratch <- array(NA, p)
    result <- array(NA, p)
    if(p <= 0 ) {
        return
    }
    result[1] <- reflectionCoeffs[1]
    if(p > 1) {
        result[1] <- reflectionCoeffs[1]
        scratch[1] <- result[1]
        for(kM1 in 2:p) {
            result[kM1] <- reflectionCoeffs[kM1]
            ##;;; Use Levinson-Durbin recursions to generate
            ##;;; the remaining (k)-th order AR coefficients
            ##print(kM1)
            ##print(result)
            scratch[1:kM1] <- result[1:kM1]
            ##print(scratch)
            for(j in 1:(kM1-1)) {
                ##print(paste("j", j))
                result[j] <- scratch[j] -
                    (result[kM1] * scratch[kM1 - j])
                ##print(paste("j res j", j, result[j]))
            }
        }
    }
    return(result)

}

yuleWalkerAlgorithmGivenData <- function(timeSeries, p,
                                         centreData=T) {
    ##   "given
    ##    [1] time-series (required)
    ##        ==> a vector containing a real-valued time series x_t
    ##    [2] p (required)
    ##        ==> autoregressive model order;
    ##            p should be an integer > 0 and < end - start
    ##    [3] start (keyword; 0)
    ##        ==> start index of vector to be used
    ##    [4] end (keyword; length of the-seq)
    ##        ==> 1 + end index of vector to be used
    ##    [5] center-data (keyword; t)
    ##        ==> if t, subtract sample mean from time series;
    ##            if a number, subtracts that number from time series;
    ##            if nil, the effect of this function is to return
    ##            a copy of the relevant portion of time-series
    ##    [6] AR-coeffs (keyword; a vector of length p)
    ##        <== estimates of AR coefficients phi_{1,p}, ..., phi_{p,p}
    ##    [7] approximate-MSEs (keyword; a vector of length p+1)
    ##        <== estimates of innovations variance (mean square errors)
    ##            for models of order 0, 1, ..., p (the innovations variance
    ##            for the zeroth order model is the sample variance)
    ##    [8] reflection-coeffs (keyword; a vector of length p)
    ##        <== estimates of reflection coefficients phi_{1,1}, ..., phi_{p,p}
    ##    [9] forward-prediction-errors (keyword; a vector of length end - start - p)
    ##        <== computed forward prediction errors
    ##   [10] backward-prediction-errors (keyword; a vector of length end - start - p)
    ##        <== computed backward prediction errors
    ## uses the Yule-Walker method to estimate the phi's in the model
    ##    x_t = phi_{1,p} x_{t-1} + ... + phi_{p,p} x_{t-p} + e_t
    ## and returns
    ##    [1] AR-coeffs, a vector containing phi_{1,p}, ..., phi_{p,p}
    ##    [2] estimate of the innovations variance, i.e.,
    ##        the variance of the e_t's; this number is also given in
    ##        (elt approximate-MSEs p)
    ##    [3] approximate-MSEs
    ##    [4] reflection-coeffs
    ##    [5] forward-prediction-errors
    ##    [6] backward-prediction-errors
    ## ---
    ## Note: see p. 420 of the SAPA book"
    approximateMSEs <-  array(NA, p +1)
    reflectionCoeffs <- array(NA, p)
    nTs <- length(timeSeries)
    if(p%%1 != 0 || p <= 0 || p > nTs) {
        return
    }
    nTsPlusp <- nTs+p
    forwardJ <- array(0, nTsPlusp)
    ##backwardJ <- array(NA, nTsPlusp)
    if(!is.null(centreData)) {
        if(is.numeric(centreData)) {
            timeSeries <- timeSeries - centreData
        }
        else if(centreData==T) {
            timeSeries <- timeSeries - mean(timeSeries)
        }
    }
    forwardJ[1:nTs] <- timeSeries
    scratch <- array(NA, nTsPlusp)
    backwardJ <- c(forwardJ[nTsPlusp], forwardJ[1:(nTsPlusp-1)])
    ##print(length(forwardJ))
    ##print(length(backwardJ))
    
    ##;;; zeros at end don't matter ...
    approximateMSEs[1] <- sum(forwardJ^2) / nTs

    for( j in 1:p ) {
        reflectionCoeffs[j] <- (t(backwardJ)%*%forwardJ) /
            (t(backwardJ)%*%backwardJ)
        approximateMSEs[j +1] <- approximateMSEs[j] *
            (1.0 - reflectionCoeffs[j]^2)
        scratch <- (- reflectionCoeffs[j] * backwardJ) + forwardJ
        backwardJ <- (- reflectionCoeffs[j] * forwardJ) + backwardJ
        backwardJ <- c(backwardJ[nTsPlusp], backwardJ[1:(nTsPlusp-1)])
        forwardJ <- scratch
        
    }
    ##;;; At this point we have the reflection coefficients and
    ##;;; approximate mse's, from which we can now generate
    ##;;; the AR coefficients.

    arCoeffs <- reflectionCoeffsToARCoeffs(reflectionCoeffs, p=p)
    return(list(arCoeffs=arCoeffs,
                extInnovVar=approximateMSEs[p+1],
                approximateMSEs=approximateMSEs,
                reflectionCoeffs=reflectionCoeffs,
                forwardPredictionErrors=forwardJ[(p+1):nTs],
                backwardPredictionErrors=backwardJ[(p+2):(nTs+1)]))
}

tmp <- yuleWalkerAlgorithmGivenData(c(71, 63, 70, 88, 99, 90, 110, 135, 128, 154,
                                156, 141, 131, 132, 141, 104, 136, 146, 124,
                                129), 3)   
##timeSeries=c(71, 63, 70, 88, 99, 90, 110, 135, 128, 154,
##                                156, 141, 131, 132, 141, 104, 136, 146, 124,
##                                129)
##p=3
##centreData=NULL

## reflectionCoeffsToARCoeffs <- function(reflectionCoeffs,
##                                        p=length(reflectionCoeffs)) {
##     ##       "given
##     ##    [1] reflection-coeffs (required)
##     ##        ==> vector of length p with real-valued 
##     ##            reflection coefficients
##     ##            phi_{1,1}, phi_{2,2}, ... , phi_{p,p}
##     ##    [2] p (keyword; length of reflection-coeffs)
##     ##        ==> AR model order
##     ##
##     ##   results (keyword; vector of length p)
##     ##        <== vector of size p into which are placed
##     ##            the AR coefficients phi_{j,p} for model
##     ##            x_t = phi_{1,p} * x_{t-1} + phi_{2,p} * x_{t-2}
##     ##                  + ... + phi_{p,p} * x_{t-p} + e_t
##     ## returns
##     ##    [1] results, i.e., the AR coefficients
##     ## ---
##     ##   Note: see Section 9.4 of the SAPA book"
##     scratch <- array(NA, p)
##     result <- array(NA, p)
##     if(p <= 0 ) {
##         return
##     }
##     result[1] <- reflectionCoeffs[1]
##     if(p > 1) {
##         result[1] <- reflectionCoeffs[1]
##         scratch[1] <- result[1]
##         for(kM1 in 2:p) {
##             result[kM1] <- reflectionCoeffs[kM1]
##             ##;;; Use Levinson-Durbin recursions to generate
##             ##;;; the remaining (k)-th order AR coefficients
##             scratch[1:kM1] <- result[1:kM1]
##             for(j in 1:(kM1-1)) {
##                 result[j] <- scratch[j] -
##                     (result[kM1] * scratch[kM1 - j])
##             }
##         }
##     }
##     return(result)

## }

reflectionCoeffsToARCoeffs(c(0.771228137382831,
                             0.10369838100465757,
                             0.16589516290564163,
                             -0.4428652096133231))
## [1]  0.74751916 -0.01583855  0.46440844 -0.44286521
##reflectionCoeffs <- arrayuleWalkerAlgorithmGivenACVS(acvs(c(88, 99, 90, 110, 135, 128))$acvs, 3)
yuleWalkerAlgorithmGivenACVS <- function(acvs, p) {
    ##       "given
    ##    [1] acvs (required)
    ##        ==> a vector containing values of the acvs
    ##            for a real-valued time series x_t
    ##            from lag 0 up to (and including) p
    ##    [2] p (required)
    ##        ==> autoregressive model order;
    ##            p should be an integer > 0 and < end - start
    ##    [3] AR-coeffs (keyword; a vector of length p)
    ##        <== estimates of AR coefficients phi_{1,p}, ..., phi_{p,p}
    ##    [4] approximate-MSEs (keyword; a vector of length p+1)
    ##        <== estimates of innovations variance (mean square errors)
    ##            for models of order 0, 1, ..., p (the innovations variance
    ##            for the zeroth order model is the sample variance)
    ##    [5] reflection-coeffs (keyword; a vector of length p)
    ##        <== estimates of reflection coefficients phi_{1,1}, ..., phi_{p,p}
    ## uses the Yule-Walker method to estimate the phi's in the model
    ##    x_t = phi_{1,p} x_{t-1} + ... + phi_{p,p} x_{t-p} + e_t
    ## and returns
    ##    [1] AR-coeffs, a vector containing phi_{1,p}, ..., phi_{p,p}
    ##    [2] estimate of the innovations variance, i.e.,
    ##        the variance of the e_t's; this number is also given in
    ##        (elt approximate-MSEs p)
    ##    [3] approximate-MSEs
    ##    [4] reflection-coeffs
    ## ---
    ## Note: see Sections 9.3 and 9.4 of the SAPA book"
    if(!is.integer(p) && p < 1 && p < length(acvs)) {
        return
    }
    ARcoeffs <- array(NA, p)
    approxMSEs <- array(NA, p+1)
    reflectionCoeffs <- array(NA, p)
    approxMSEs[1] <- acvs[1]

    scratch <- array(NA, p)
    top <- NULL
    k <- NULL
    ##   ;;;
    ##   ;;; Main computational loop - once for each AR order
    ##   ;;;
    for( km1 in 0:(p-1)) {
        k <- km1 +1
        top <- acvs[k+1]
        if(km1 >= 1) {
            for(j in 1:km1) {
                top <-  top - ARcoeffs[j] * acvs[k - j + 1]
            }
        }
        ##  ;;; Calculate k-th order reflection coefficient
        ARcoeffs[k] <- top / approxMSEs[k]
        reflectionCoeffs[k] <- ARcoeffs[k]
        ##  ;;; Update mean square error using approximate method.
        approxMSEs[k +1] <- approxMSEs[k] *
            (1 - abs(ARcoeffs[k])**2)
        ##  ;;; Use Levinson-Durbin recursions to generate
        ##  ;;; the remaining (k)-th order AR coefficients
        if(km1 >= 1) {
            for( j in 1:km1) {
                scratch[j] <- ARcoeffs[j]
            }
        }
        if(km1 >= 1) {
            for( j in 1:km1) {
                ARcoeffs[j] <- scratch[j] -
                    ARcoeffs[k] * scratch[k - j ]
            }
        }
    }
        
    return(list(ARcoeffs=ARcoeffs, estInnovVar=approxMSEs[p],
                approxMSEs=approxMSEs, reflectionCoeffs=reflectionCoeffs))
}


##yuleWalkerAlgorithmGivenACVS(acvs(c(71, 63, 70, 88, 99, 90, 110, 135, 128, 154, 156, 141, 131, 132, 141, 104, 136, 146, 124, 129))$acvs, 3)
##yuleWalkerAlgorithmGivenACVS(acvs(c(88, 99, 90, 110, 135, 128))$acvs, 3)
##compare with ar.yw(c(71, 63, 70, 88, 99, 90, 110, 135, 128, 154, 156, 141, 131, 132, 141, 104, 136, 146, 124, 129), aic=F, order=3)


arCoeffsToPrewhitenFilter <- function(arCoeffs) {
    return(c(1, -arCoeffs))
}

arCoeffsToSdf <- function(arCoeffs,
                          innovationsVariance,
                          deltaT,
                          nNonZeroFreqs) {
    
    prewhitenFilter <- c(1, -arCoeffs)
    scratch <- c(prewhitenFilter, rep(0.0, 2 *
                                  (nNonZeroFreqs - 1) -
                                  length(prewhitenFilter)))
    scratch <- fft(scratch)
    numerator <- innovationsVariance * deltaT
    resultSdf <- numerator / abs(scratch[1:(1+nNonZeroFreqs)])**2
    return(resultSdf)
}

## arCoeffs <- c(0.6990553262846771,
##               -0.005832588617464066,
##               0.15582353261696208)
## innovVar <-  271.77271305830277
## deltaT <-  0.25

## ## scratch <- c(c(1, -1 * arCoeffs), rep(0.0, 2*256-length(arCoeffs+1)))
## ## scratch <- fft(scratch)
## ## numerator <- innovVar * deltaT

## ## resultSdf <- numerator / abs(scratch[1:257])**2
## ## resultSdf <- 10*log10(resultSdf)
##  10*log10(arCoeffsToSdf(arCoeffs, innovVar, deltaT, 256))
## tmp2 <- c(71, 63, 70, 88, 99, 90, 110, 135, 128, 154, 156, 141, 131, 132, 141, 104, 136, 146, 124, 129)
## tmp2 <- tmp2 - mean(tmp2)

## arcoef1 <- c(0.7475191644352183, -0.01583855174510637, 0.46440844342490983, -0.4428652096133231)

## prewhiten <- c(1, -1 * arcoef1)

## nfft <- 4*2**(ceiling(log2(length(tmp2))))

## dataFFT <- c(tmp2, rep(0.0, nfft-length(tmp2)))
## coeffsFFT <- c(prewhiten, rep(0.0, nfft-length(prewhiten)))

## dataFFT <- fft(dataFFT)
## coeffsFFT <- fft(coeffsFFT)
## res <- Conj(dataFFT*coeffsFFT)/nfft
## res <- fft(res)



## burgAlgorithm <- function(timeSeries, p, centreData=T,
                          
## (defun burg-algorithm
##        (time-series
##         p
##         &key
##         (start 0)
##         (end (length time-series))
##         (center-data t)  ;t, nil or value to be subtracted off ...
##         (AR-coeffs (make-array p))
##         (approximate-MSEs (make-array (1+ p)))
##         (reflection-coeffs (make-array p))
##         (forward-prediction-errors (make-array (- end start p)))
##         (backward-prediction-errors (make-array (- end start p)))
##         (exact-MSEs (make-array (1+ p))))
##   "given
##    [1] time-series (required)
##        ==> a vector containing a real-valued or
##            complex-valued time series x_t
##    [2] p (required)
##        ==> autoregressive model order;
##            p should be an integer > 0 and < end - start
##    [3] start (keyword; 0)
##        ==> start index of vector to be used
##    [4] end (keyword; length of the-seq)
##        ==> 1 + end index of vector to be used
##    [5] center-data (keyword; t)
##        ==> if t, subtract sample mean from time series;
##            if a number, subtracts that number from time series;
##            if nil, the effect of this function is to return
##            a copy of the relevant portion of time-series
##    [6] AR-coeffs (keyword; a vector of length p)
##        <== estimates of AR coefficients phi_{1,p}, ..., phi_{p,p}
##    [7] approximate-MSEs (keyword; a vector of length p+1)
##        <== estimates of innovations variance (mean square errors)
##            for models of order 0, 1, ..., p (the innovations variance
##            for the zeroth order model is the sample variance)
##    [8] reflection-coeffs (keyword; a vector of length p)
##        <== estimates of reflection coefficients phi_{1,1}, ..., phi_{p,p}
##    [9] forward-prediction-errors (keyword; a vector of length end - start - p)
##        <== computed forward prediction errors
##   [10] backward-prediction-errors (keyword; a vector of length end - start - p)
##        <== computed backward prediction errors
##   [11] exact-MSEs (keyword; a vector of length p+1)
##        <== another set of estimates of innovations variance
##            for models of order 0, 1, ..., p; these estimates
##            are based on Equation (419a) of the SAPA book
## uses Burg's algorithm to estimate the phi's in the model
##    x_t = phi_{1,p} x_{t-1} + ... + phi_{p,p} x_{t-p} + e_t
## and returns
##    [1] AR-coeffs, a vector containing phi_{1,p}, ..., phi_{p,p}
##    [2] estimate of the innovations variance, i.e.,
##        the variance of the e_t's; this number is also given in
##        (elt approximate-MSEs p)
##    [3] approximate-MSEs
##    [4] reflection-coeffs
##    [5] forward-prediction-errors
##    [6] backward-prediction-errors
##    [7] exact-MSEs
## ---
## Note: see Section 9.5 of the SAPA book"
##   (let ((N-ts (- end start)))
##     (assert (and (integerp p) (> p 0) (< p N-ts)))
##     (let* ((forward-scratch (center&taper-time-series
##                              time-series
##                              :center-data center-data
##                              :start start
##                              :end end))
##            (backward-scratch (copy-seq forward-scratch))
##            (sh2-old (float (/ (sum-of-squares forward-scratch)
##                               N-ts)))
##            (mseb-old sh2-old)
##            (scratch (make-array p))
##            top bottom k)
##       ;(declare (dynamic-extent forward-scratch backward-scratch))
##       ;;; Set zeroth order mean square prediction error
##       (setf (aref approximate-MSEs 0)
##             (setf (aref exact-MSEs 0) sh2-old))
##       ;;;
##       ;;; Main computational loop - once for each AR order
##       ;;;
##       (dotimes (k-minus-1 p)
##         (setf k (1+ k-minus-1)
##               top 0.0
##               bottom 0.0)
##         (dotimes
##           (j (- N-ts k))
##           (incf top 
##                 (* (aref forward-scratch (1+ j))
##                    (conjugate (aref backward-scratch j))))
##           (incf bottom
##                 (+ (expt (abs (aref forward-scratch (1+ j))) 2)
##                    (expt (abs (aref backward-scratch j)) 2))))
##         ;;; Calculate k-th order reflection coefficient
##         (setf (aref reflection-coeffs k-minus-1)
##               (setf (aref AR-coeffs k-minus-1)
##                     (/ (* 2.0 top) bottom)))
##         ;;; Update mean square error using:
##         ;;; (1) approximate method
##         (setf (aref approximate-MSEs k)
##               (setf sh2-old
##                     (* sh2-old
##                        (- 1.0
##                           (expt (abs (aref AR-coeffs k-minus-1)) 2)))))
##         ;;; (2) exact method
##         (setf (aref exact-MSEs k)
##               (setf mseb-old
##                     (* (/ (- (* mseb-old (1+ (- N-ts k)))
##                              (/ (+ (expt (abs (aref forward-scratch 0)) 2)
##                                    (expt (abs (aref backward-scratch (- N-ts k))) 2))
##                                 2.0))
##                           (- N-ts k))
##                        (- 1.0
##                           (expt (abs (aref AR-coeffs k-minus-1)) 2)))))
##         ;;; Update prediction errors
##         (dotimes
##           (j (- N-ts k))
##           (setf (aref forward-scratch j)
##                 (- (aref forward-scratch (1+ j))
##                    (* (aref AR-coeffs k-minus-1)
##                       (aref backward-scratch j))))
##           (decf (aref backward-scratch j)
##                 (* (conjugate (aref AR-coeffs k-minus-1))
##                    (aref forward-scratch (1+ j)))))
##         ;;; Use Levinson-Durbin recursions to generate
##         ;;; the remaining (k)-th order AR coefficients
##         (cond
##          ((> k 1)
##           (dotimes
##             (j k-minus-1)
##             (setf (aref scratch j) (aref AR-coeffs j)))
##           (dotimes
##             (j k-minus-1)
##             (setf (aref AR-coeffs j)
##                   (- (aref scratch j)
##                      (* (aref AR-coeffs k-minus-1)
##                         (conjugate (aref scratch (- k-minus-1 j 1))))))))))
##       (values AR-coeffs
##               (elt approximate-MSEs p)
##               approximate-MSEs
##               reflection-coeffs
##               (copy-vector forward-scratch forward-prediction-errors
##                            :end (- N-ts p))
##               (copy-vector backward-scratch backward-prediction-errors
##                            :end (- N-ts p))
##               exact-MSEs))))

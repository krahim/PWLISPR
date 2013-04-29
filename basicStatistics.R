# The following functions are contained in basicStatistics.R
# These functions are loosley based on the lisp code provided with SAPA
# http://lib.stat.cmu.edu/sapaclisp/
# http://lib.stat.cmu.edu/sapaclisp/basic-statistics.lisp

# sampleVariance.bias <-function(timeSeries)

################################################################################



#S_{N} is what is used not S_{N-1} in the SAPA Lisp code
sampleVariance.bias <-function(timeSeries) {
   N <- length(timeSeries);
   return((as.double(t(timeSeries)%*%timeSeries) - sum(timeSeries)^2/N)/N);
}
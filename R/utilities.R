# The following functions are contained in utilities.R
# These functions are loosley based on the lisp code provided with SAPA
# http://lib.stat.cmu.edu/sapaclisp/
# http://lib.stat.cmu.edu/sapaclisp/utilities.lisp

# convertTodB <- function(X)

# convertFromdB <- function(X)

# powerOf2 <- function(x)

# nextPowerOf2 <- function(x)

# NyquistFrequency  <- function(sampleTime)

# sampleTime <- function(NyquistFrequency)

# compareSeqs<- function(oneSeq, anotherSeq)

################################################################################

#libraries required
##library(bitops);

################################################################################


convertTodB <- function(X) {
   if(!identical(X[X<=0], numeric(0))) {
      return();
   }
   return(10*log10(X));
}

convertFromdB <- function(X) {
   return(10^(X/10));
}

#library(bitops)

powerOf2 <- function(x) {
   if(as.integer(x)!=x) {
      return();
   }
   if(bitAnd(x,x-1) == 0) {
      return(log2(x));
   }
   else {
      return();
   }
}

nextPowerOf2 <- function(x) {
   return(2^(ceiling(log2(x))));
}

NyquistFrequency  <- function(sampleTime) {
   return(1/(2*sampleTime));
}

sampleTime <- function(NyquistFrequency) {
   return(1/(2*NyquistFrequency));
}

#  "given
#   [1] one-seq (required)
#       ==> any sequence
#   [2] another-seq (required)
#       ==> any sequence (must be same length as one-seq)
#returns
#   [1] maximum absolute difference between corresponding
#       elements of the two sequences
#   [2] average absolute difference
#---
#Note: useful for comparing results that in theory
#should be identical"
compareSeqs<- function(oneSeq, anotherSeq) {
   if(length(oneSeq) != length(anotherSeq)) {
      return();
   }
   absDiff <- abs(oneSeq - anotherSeq);
   maxAbsDiff <- max(absDiff);
   avgAbsDiff <- mean(absDiff);
   return(list(avgAbsDiff=avgAbsDiff,maxAbsDiff=maxAbsDiff));
}
#compareSequence( c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
#   c(1, 2, 3, 4, 5, 6.1, 7, 8, 9, 10))

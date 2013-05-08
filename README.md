PWLISPR
===========

R code implementing tools from Percival and Walden (1993). This is based on their LISP code. This code can make use of the filtering repo which uses FFTW and requires compiling. 

I have included the sine multitaper spectral estimate coded in R in here. I must have added this code some time ago and it is included in this package as it may depend on underlying routines. 

This is currently an R package that I do not plan to put on CRAN, and all documentation is in the form of source code comments.

Build instructions for general Linux machine.

1) download and unzip to a folder called PWLISPR-master 2) from the parent folder: a) R CMD build PWLISPR-master/ b) R CMD check PWLISPR_0.1-1.tar.gz (optional) c) R CMD INSTALL PWLISPR_0.1-1.tar.gz

Note: The version number may change, and you will likely have to set the PATH variable for other operating systems.
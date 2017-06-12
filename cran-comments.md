## Test environments

* local windows 7 install, R 3.3.1
* local ubuntu 16.04 install, R 3.4.0
* ubuntu 12.04 (on travis-ci), R 3.3.2
* osx 10.11.6 (on travis-ci), R 3.3.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

Win-builder devel and release note:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Joris Chau <j.chau@uclouvain.be>'

Ubuntu 16.04 note:

* checking installed package size ... NOTE
  installed size is  7.2Mb
  sub-directories of 1Mb or more:
    libs   6.9Mb
      
My understanding is that the inflation of the libs subdirectory is due to the use of `Rcpp` and `RcppArmadillo`. 



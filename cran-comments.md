## Test environments

* Windows 7, R 3.3.1, pass
* Linux Ubuntu 17.04 install, R 3.4.2, 1 note
* Linux CentOS 6.9, R 3.3.1, pass
* MacOS Sierra 10.12.6 (on travis-ci), R 3.4.3, pass
* Win-builder (R-oldrelease, release and devel), 1 note

## R CMD check results

0 errors | 0 warnings | 1 note

**Win-builder note:**

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Joris Chau <j.chau@uclouvain.be>'

Possibly mis-spelled words in DESCRIPTION:
  Hermitian (3:32, 7:3, 10:6, 11:30)
  covariance (7:64)

**Ubuntu 17.04 note:**

* checking installed package size ... NOTE
  installed size is  5.3Mb
  sub-directories of 1Mb or more:
    libs  4.2Mb
      
My understanding is that the inflation of the libs subdirectory is due to the use of `Rcpp` and `RcppArmadillo`. 



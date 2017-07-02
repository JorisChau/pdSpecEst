## Resubmission
This is a resubmission. In this version I have: 

* changed all C++ calls 'std::sqrt(x)' to 'std::sqrt((double)x)' removing the ERROR `call of overloaded ‘std::sqrt(int)’ is ambiguous` on Solarix-x86. 

* removed link to `shapes::distcov` in .Rd file (`pdDist.Rd`) so that the `shapes` package no longer needs to be imported. This removes the NOTE on Linux-fedora.

## Test environments

* local windows 7 install, R 3.3.1, pass
* local ubuntu 16.04 install, R 3.4.0, 1 note
* ubuntu 12.04 (on travis-ci), R 3.3.2, pass
* osx 10.11.6 (on travis-ci), R 3.3.2, pass
* win-builder (devel and release), 1 note

## R CMD check results

0 errors | 0 warnings | 1 note

**Win-builder devel and release note:**

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Joris Chau <j.chau@uclouvain.be>'

Possibly mis-spelled words in DESCRIPTION:
  Hermitian (3:29, 7:3, 9:60, 11:3)
  covariance (7:64)

**Ubuntu 16.04 (local) and Ubuntu 12.04 (travis-ci) note:**

* checking installed package size ... NOTE
  installed size is  8.2Mb
  sub-directories of 1Mb or more:
    libs  7.8Mb
      
My understanding is that the inflation of the libs subdirectory is due to the use of `Rcpp` and `RcppArmadillo`. 



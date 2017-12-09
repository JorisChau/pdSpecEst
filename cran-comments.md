## Resubmission

This is a resubmission. In this version I have resolved the note:
* checking dependencies in R code ... NOTE
Namespaces in Imports field not imported from:
  ‘R.utils’ ‘grid’
  All declared Imports should be used.

## Test environments

* Windows 7, R 3.3.1, pass
* Linux Ubuntu 17.04, R 3.4.3, 1 note
* Linux CentOS 6.9, R 3.3.1, pass
* Linux Ubuntu 14.04 (on travis-ci), R 3.4.2, 1 note
* MacOS Sierra 10.12.6 (on travis-ci), R 3.4.3, pass
* Win-builder (oldrelease, release and devel), 1 note

## CRAN package version 1.1.1

Since RcppArmadillo 0.8.300.1.0 (2017-12-06), the calls to the armadillo function `arma::chol` break on Linux (cannot find the LAPACK routine `zpbtrf`). I removed all calls to `arma::chol` to fix the errors in the CRAN package check results.

## R CMD check results

0 errors | 0 warnings | 1 note

**Win-builder oldrelease and release note:**

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Joris Chau <j.chau@uclouvain.be>'

Possibly mis-spelled words in DESCRIPTION:
  Hermitian (3:32, 7:3, 10:6, 11:30)
  covariance (7:64)

**Ubuntu 17.04 and 14.04 note:**

* checking installed package size ... NOTE
  installed size is  5.3Mb
  sub-directories of 1Mb or more:
    libs  4.2Mb
      
My understanding is that the inflation of the libs subdirectory is due to the use of `Rcpp` and `RcppArmadillo`. 



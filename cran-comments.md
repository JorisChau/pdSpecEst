## Test environments
* local windows 7 install, R 3.3.1
* Ubuntu 12.04 (on travis-ci), R 3.3.2
* OS X 10.11.6 (on travis-ci), R 3.3.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

Win-builder devel and release note:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Joris Chau <j.chau@uclouvain.be>'

Travis-ci note under Ubuntu 12.04:

* checking installed package size ... NOTE
    installed size is  8.1Mb
    sub-directories of 1Mb or more:
      libs   7.9Mb
- My understanding is that the inflation of the libs subdirectory is due to the use of Rcpp and RcppArmadillo. 


## Reverse dependencies

This is a new release, so there are no reverse dependencies.


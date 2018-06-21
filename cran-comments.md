## CRAN package version 1.2.2

### Test environments

* Ubuntu Linux 17.10 (local install) R 3.4.2, and 14.04 (travis-ci) R-release, 1 note (see below)
* Debian Linux, R-release, gcc (rhub), 1 note (see below)
* CentOS 6 with Redhat Developer Toolset, R from EPEL (rhub), 1 note (see below)
* Oracle Solaris 10, x86, 32 bit, R-patched (rhub), pass
* macOS 10.11 El Capitan, R-release (rhub and travis-ci), pass

* Win-builder R-release, R-devel and R-oldrel, 1 note

checking CRAN incoming feasibility ... NOTE
Maintainer: 'Joris Chau <j.chau@uclouvain.be>'

Possibly mis-spelled words in DESCRIPTION:
  Hermitian (3:32, 7:3, 9:68, 11:12, 13:47)
  denoising (10:63)

* Ubuntu, Debian, Fedora and CentOS Linux note:

checking installed package size ... NOTE
  installed size is  7.4Mb
  sub-directories of 1Mb or more:
    libs   6.5Mb
R CMD check results
0 errors | 0 warnings | 1 note 

My understanding is that the inflation of the libs subdirectory is due to the use of `Rcpp` and `RcppArmadillo`. 

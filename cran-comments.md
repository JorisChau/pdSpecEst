## CRAN package version 1.2.3

### Comments

- Fixed WARNING:
     checkRd: (-1) ...: Non-ASCII contents without declared encoding
  by removing non-ASCII characters from REFERENCES.bib and setting 
  UTF-8 encoding in DESCRIPTION file 
- Fixed NOTE:
     src/Makevars.win: SHLIB_OPENMP_CXXFLAGS is included in PKG_CXXFLAGS but not in PKG_LIBS
     src/Makevars.win: SHLIB_OPENMP_CFLAGS is included in PKG_LIBS but not in PKG_CFLAGS


### Test environments

* Ubuntu Linux 18.04 (local install) R-release, and 16.04 (rhub) R-release + R-devel, 2 notes (see below)
* Debian Linux, R-release + R-devel, gcc + gcc ASAN/UBSAN (rhub), 2 notes (see below)
* CentOS 6 with Redhat Developer Toolset, R from EPEL (rhub), pass
* Oracle Solaris 10, x86, 32 bit, R-patched (rhub), pass
* macOS 10.11 El Capitan + 10.9 Mavericks, R-release + R-oldrel (travis-ci, rhub), pass
* Win-builder R-release, R-devel and R-oldrel, 1 note

### Notes

* checking CRAN incoming feasibility ... NOTE
New maintainer:
  Joris Chau <joris.chau@openanalytics.eu>
Old maintainer(s):
  Joris Chau <j.chau@uclouvain.be>

* Possibly mis-spelled words in DESCRIPTION:
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

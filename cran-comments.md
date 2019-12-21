## CRAN package version 1.2.4

### Comments

* Fixed ERROR:
    Running 'testthat.R' [5s/6s]
    Running the tests in 'tests/testthat.R' failed.
  by editing test that checked against R's messages  
    
* Updated references    
    
### Test environments

* Ubuntu Linux 18.10 (local install) R-release, and 16.04 (rhub) R-release + R-devel, 1 note (see below)
* Fedora Linux R-devel (rhub), 2 notes (see below)
* Debian Linux R-release + R-devel (rhub), 1 note (see below)
* Win-builder R-release, R-devel, pass

### Notes

Fedora, Debian Linux note:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Joris Chau <joris.chau@openanalytics.eu>’

Possibly mis-spelled words in DESCRIPTION:
  denoising (10:63)
  
Ubuntu, Fedora Linux note:

* checking installed package size ... NOTE
  installed size is  7.4Mb
  sub-directories of 1Mb or more:
    libs   6.5Mb


My understanding is that the inflation of the libs subdirectory is due to the use of `Rcpp` and `RcppArmadillo`. 

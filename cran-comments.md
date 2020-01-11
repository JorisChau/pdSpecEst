## CRAN package version 1.2.5

### Comments

CRAN checks produce an error on all Windows flavors:
  
  Error in library.dynam(lib, package, package.lib) : 
     DLL 'stringi' not found: maybe not installed for this architecture?   
     
Cannot reproduce error in Windows (local install) R-release or Win-builder R-release, R-devel, R-old-release
R-pkg-devel mailing list suggests this might be due to a race condition and it is advised to resubmit the package.


### Test environments

* Ubuntu Linux 18.10 (local install) R-release, and 16.04 (rhub) R-release + R-devel, 1 note (see below)
* Fedora Linux R-devel (rhub), 2 notes (see below)
* Debian Linux R-release + R-devel (rhub), 1 note (see below)
* Windows (local install) + Win-builder R-release, R-devel, R-old-release pass

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

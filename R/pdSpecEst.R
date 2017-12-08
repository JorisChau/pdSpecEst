#' pdSpecEst: An Analysis Toolbox for Hermitian Positive Definite Matrices
#'
#' The \code{pdSpecEst} (\strong{p}ositive \strong{d}efinite \strong{Spec}tral \strong{Est}imation)
#' package provides data analysis tools for samples of symmetric or Hermitian positive definite matrices,
#' such as collections of (non-degenerate) covariance matrices or spectral density matrices.
#'
#' The tools in this package can be used to perform:
#' \itemize{
#'    \item \emph{Intrinsic wavelet regression} and \emph{clustering} for curves of Hermitian positive definite matrices.
#'    These implementations are based on the paper (Chau and von Sachs, 2017a).
#'    \item \emph{Exploratory data analysis} and \emph{inference} for samples of Hermitian positive definite matrices by
#'    means of intrinsic data depth and intrinsic rank-based hypothesis tests. These implementations are based on the paper
#'    (Chau, Ombao and von Sachs, 2017b).
#'  }
#' For more details and examples on how to use the package see the accompanying vignettes in the vignettes folder.
#' A demo Shiny app to test out the implemented functions in the package is available
#' \href{https://jchau.shinyapps.io/pdSpecEst/}{here}.
#'
#' Author and maintainer: \strong{Joris Chau} (\email{j.chau@@uclouvain.be}).
#'
#' Install the current development version via \code{devtools::install_github("JorisChau/pdSpecEst")}.
#'
#' @references Chau, J. and von Sachs, R. (2017a). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Chau, J., Ombao, H., and von Sachs, R. (2017b). \emph{Data depth and rank-based
#' tests for covariance and spectral density matrices}. Available at \url{http://arxiv.org/abs/1706.08289}.
#'
#' @docType package
#' @useDynLib pdSpecEst, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import utils
#' @import stats
#' @name pdSpecEst
NULL
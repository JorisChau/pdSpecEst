% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{pdParTrans}
\alias{pdParTrans}
\title{Riemannian HPD parallel transport}
\usage{
pdParTrans(P, V, W)
}
\arguments{
\item{P}{a \eqn{(d,d)}-dimensional HPD matrix.}

\item{V}{a \eqn{(d,d)}-dimensional Hermitian matrix corresponding to a vector in the tangent space of \code{P}.}

\item{W}{a \eqn{(d,d)}-dimensional Hermitian matrix corresponding to a vector in the tangent space of \code{P}.}
}
\value{
a \eqn{(d,d)}-dimensional Hermitian matrix corresponding to the parallel transportation of \code{W} in
the direction of \code{V} along a geodesic curve for a unit time step.
}
\description{
\code{pdParTrans} computes the parallel transport on the manifold of HPD matrices
equipped with the affine-invariant Riemannian metric as described in e.g., Chapter 2 of \insertCite{C18}{pdSpecEst}. That is,
the function computes the parallel transport of a Hermitian matrix \code{W} in the tangent space
at the HPD matrix \code{P} along a geodesic curve in the direction of the Hermitian matrix \code{V}
in the tangent space at \code{P} for a unit time step.
}
\examples{
## Transport the vector W to the tangent space at the identity
W <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
diag(W) <- rnorm(3)
W[lower.tri(W)] <- t(Conj(W))[lower.tri(W)]
p <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
P <- t(Conj(p)) \%*\% p

pdParTrans(P, Logm(P, diag(3)), W) ## whitening transport

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{Expm}, \link{Logm}}
}

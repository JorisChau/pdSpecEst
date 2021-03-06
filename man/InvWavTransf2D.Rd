% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/invwavtrans.R
\name{InvWavTransf2D}
\alias{InvWavTransf2D}
\title{Inverse AI wavelet transform for surface of HPD matrices}
\usage{
InvWavTransf2D(D, M0, order = c(3, 3), jmax, metric = "Riemannian",
  ...)
}
\arguments{
\item{D}{a list of arrays containing the 2D pyramid of wavelet coefficients, where each array contains the
(\eqn{d,d})-dimensional wavelet coefficients from the coarsest wavelet scale \code{j = 0} up to the finest
wavelet scale \code{j = jmax}. This is the same format as the \code{$D} component given as output by
 \code{\link{WavTransf2D}}.}

\item{M0}{a numeric array containing the midpoint(s) at the coarsest scale \code{j = 0} in the 2D midpoint pyramid.
This is the same format as the \code{$M0} component given as output by \code{\link{WavTransf2D}}.}

\item{order}{a 2-dimensional numeric vector \eqn{(1,1) \le} \code{order} \eqn{\le (9,9)} corresponding to the marginal
orders of the intrinsic 2D AI refinement scheme, defaults to \code{order = c(3, 3)}.}

\item{jmax}{the maximum scale (resolution) up to which the 2D surface of HPD midpoints (i.e. scaling coefficients) are
reconstructed. If \code{jmax} is not specified it is set equal to the resolution in the finest wavelet scale
\code{jmax = length(D)}.}

\item{metric}{the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
\code{"Euclidean"}. See also the Details section below.}

\item{...}{additional arguments for internal use.}
}
\value{
Returns a (\eqn{d, d, n_1, n_2})-dimensional array corresponding to a rectangular surface of size \eqn{n_1} by
\eqn{n_2} of (\eqn{d,d})-dimensional HPD matrices.
}
\description{
\code{InvWavTransf2D} computes the inverse intrinsic average-interpolation (AI) wavelet
transform mapping an array of coarsest-scale HPD midpoints combined with a 2D pyramid of Hermitian
wavelet coefficients to a surface in the manifold of HPD matrices equipped with a metric specified by the
user, as described in Chapter 5 of \insertCite{C18}{pdSpecEst}. This is the inverse operation of the
function \code{\link{WavTransf2D}}.
}
\details{
The input list of arrays \code{D} and array \code{M0} correspond to a 2D pyramid of wavelet coefficients and
the coarsest-scale HPD midpoints respectively, both are structured in the same way as in the output of
\code{WavTransf2D}. As in the forward AI wavelet transform, the marginal refinement orders should be smaller
or equal to 9, and the function computes the wavelet transform using a fast wavelet refinement scheme based on weighted
intrinsic averages with pre-determined weights as explained in Chapter 5 of \insertCite{C18}{pdSpecEst}. By default
\code{WavTransf2D} computes the inverse intrinsic 2D AI wavelet transform equipping the space of HPD matrices with (i)
the affine-invariant Riemannian metric as detailed in e.g., \insertCite{B09}{pdSpecEst}[Chapter 6] or \insertCite{PFA05}{pdSpecEst}.
Instead, the space of HPD matrices can also be equipped with one of the following metrics; (ii) the Log-Euclidean metric, the
Euclidean inner product between matrix logarithms; (iii) the Cholesky metric, the Euclidean inner product between Cholesky
decompositions; (iv) the Euclidean metric and (v) the root-Euclidean metric. The default choice of metric (affine-invariant Riemannian)
satisfies several useful properties not shared by the other metrics, see \insertCite{C18}{pdSpecEst} for more details. Note that this
comes at the cost of increased computation time in comparison to one of the other metrics.
}
\examples{
P <- rExamples2D(c(2^4, 2^4), 2, example = "tvar")
P.wt <- WavTransf2D(P$f) ## forward transform
P.f <- InvWavTransf2D(P.wt$D, P.wt$M0) ## backward transform
all.equal(P.f, P$f)

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{WavTransf2D}}, \code{\link{pdSpecEst2D}}, \code{\link{pdNeville}}
}

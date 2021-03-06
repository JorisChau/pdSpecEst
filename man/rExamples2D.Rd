% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/examples.R
\name{rExamples2D}
\alias{rExamples2D}
\title{Several example surfaces of HPD matrices}
\usage{
rExamples2D(n, d = 2, example = c("smiley", "tvar", "facets", "peak"),
  replicates = 1, noise = "riem-gaussian", noise.level = 1,
  df.wishart = NULL)
}
\arguments{
\item{n}{integer vector \code{c(n1, n2)} specifying the number of sampled matrices to be generated on a rectangular surface.}

\item{d}{row- (resp. column-)dimension of the generated matrices. Defaults to \code{d = 2}.}

\item{example}{the example target HPD matrix surface, one of \code{'smiley'}, \code{'tvar'}, \code{'facets'} or \code{'peak'}.}

\item{replicates}{a positive integer specifying the number of replications of noisy HPD matrix surfaces to be generated based on the
target surface of HPD matrices. Defaults to \code{replicates = 1}}

\item{noise}{noise distribution for the generated noisy surfaces of HPD matrices, one of \code{'riem-gaussian'},
\code{'log-gaussian'}, \code{'wishart'}, \code{'log-wishart'} or \code{'periodogram'}, defaults to \code{'riem-gaussian'}.
Additional details are given below.}

\item{noise.level}{parameter to tune the signal-to-noise ratio for the generated noisy HPD matrix observations.
If \code{noise.level = 0}, the noise distributions are degenerate and the noisy HPD matrix observations coincide with the target HPD matrices.
Defaults to \code{noise.level = 1}.}

\item{df.wishart}{optional parameter to specify the degrees of freedom in the case of a Wishart noise distribution (\code{noise = 'wishart'} or
\code{noise = 'log-wishart'}). By default \code{df.wishart} is equal to the dimension \code{d} to guarantee positive definiteness of the
generated noisy matrices.}
}
\value{
Returns a list with two components:
  \item{\code{f} }{ a (\code{d,d,n[1],n[2]})-dimensional array, corresponding to the \eqn{(n_1 \times n_2)}-sized example target surface of
  \eqn{(d,d)}-dimensional HPD matrices.}
  \item{\code{P} }{ a (\code{d,d,n[1],n[2]})-dimensional array, corresponding to the \eqn{(n_1 \times n_2)}-sized surface of noisy \eqn{(d,d)}-dimensional
  HPD matrices centered around the smooth target HPD matrix surface \code{f}. If \code{replicates > 1}, \code{P} is a
  \code{(d,d,n[1],n[2],length(replicates))}-dimensional array, corresponding to a collection of replicated \eqn{(n_1 \times n_2)}-sized surfaces
  of noisy \eqn{(d,d)}-dimensional HPD matrices centered around the smooth target HPD matrix surface \code{f}.}
}
\description{
\code{rExamples2D()} generates several example (locally) smooth target \emph{surfaces} of HPD matrices corrupted by
noise in a manifold of HPD matrices for testing and simulation purposes. For more details, see also Chapter 2 and 5 in
\insertCite{C18}{pdSpecEst}.
}
\details{
The examples include: (i) a \eqn{(d,d)}-dimensional \code{'smiley'} HPD matrix surface consisting of constant surfaces of random HPD matrices in
the shape of a smiley face; (ii) a \eqn{(d,d)}-dimensional \code{'tvar'} HPD matrix surface generated from a time-varying vector-auto-
regressive process of order 1 with random time-varying coefficient matrix (\eqn{\Phi}); (iii) a \eqn{(d,d)}-dimensional \code{'facets'} HPD matrix
surface consisting of several facets generated from random geodesic surfaces; and (iv) a \eqn{(d,d)}-dimensional \code{'peak'} HPD matrix surface
containing a pronounced peak in the center of its 2-d (e.g., time-frequency) domain.\cr
In addition to the (locally) smooth target surface of HPD matrices, the function also returns a noisy version of the target surface of HPD matrices, corrupted
by a user-specified noise distribution. By default, the noisy HPD matrix observations follow an intrinsic signal plus i.i.d. noise model with
respect to the affine-invariant Riemannian metric, with a matrix log-Gaussian noise distribution (\code{noise = 'riem-gaussian'}), such that the
Riemannian Karcher means of the observations coincide with the target surface of HPD matrices. Additional details can be found in Chapters 2, 3,
and 5 of \insertCite{C18}{pdSpecEst}. Other available signal-noise models include: (ii) a Log-Euclidean signal plus i.i.d. noise model, with
a matrix log-Gaussian noise distribution (\code{noise = 'log-gaussian'}); (iii) a Riemannian signal plus i.i.d. noise model, with a complex
Wishart noise distribution (\code{noise = 'wishart'}); (iv) a Log-Euclidean signal plus i.i.d. noise model, with a complex Wishart noise
distribution (\code{noise = 'log-wishart'}).
}
\examples{
example <- rExamples2D(n = c(32, 32), example = "smiley")

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{rExamples1D}}, \code{\link{pdPgram2D}}
}

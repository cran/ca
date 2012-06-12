\name{iterate.mjca}
\alias{iterate.mjca}
\title{Updating a Burt matrix in Joint Correspondence Analysis}
\description{Updating a Burt matrix in Joint Correspondence Analysis based on iteratively weighted least squares.}
\usage{iterate.mjca(B, lev.n, nd = 2, maxit = 50, epsilon = 0.0001)}
\arguments{
  \item{B      }{A Burt matrix.}
  \item{lev.n  }{The number of levels for each factor from the original response pattern matrix.}
  \item{nd     }{The required dimensionality of the solution.}
  \item{maxit  }{The maximum number of iterations.}
  \item{epsilon}{A convergence criterion for the maximum absolute difference of updated values compared to the previous values. The iteration is completed when all differences are smaller than \code{epsilon}.}
          }
\details{The function \code{iterate.mjca} computes the updated Burt matrix. This function is called from the function \code{\link{mjca}} when the option \kbd{lambda="JCA"}, i.e. when a Joint Correspondence Analysis is performed.}
\value{
  \item{B.star}{The updated Burt matrix}
  \item{crit  }{Vector of length 2 containing the number of iterations and epsilon}
      }

\seealso{\code{\link{mjca}}}
\keyword{multivariate}

\name{mjca}
\alias{mjca}
\title{Multiple and joint correspondence analysis}
\description{Computation of multiple and joint correspondence analysis.}
\usage{mjca(obj, nd = 2, lambda = "adjusted", supcol = NA, maxit = 50, epsilon = 0.0001)}
\arguments{
  \item{   obj}{A response pattern matrix containing factors.}
  \item{    nd}{Number of dimensions to be included in the output; if NA the maximum possible dimensions are included.}
  \item{lambda}{Gives the scaling method. Possible values include "indicator", "Burt", "adjusted" and "JCA".
                Using \code{lambda = "JCA"} results in a joint correspondence analysis using iterative adjusment of the Burt matrix in the solution space.}
  \item{   supcol}{Indices of supplementary columns.}
  \item{    maxit}{The maximum number of iterations (Joint Correspondence Analysis).}
  \item{  epsilon}{A convergence criterion (Joint Correspondence Analysis).}
          }
\details{The function \code{mjca} computes a multiple or joint correspondence analysis based on the eigenvalue decomposition of the Burt matrix.}
\value{
  \item{sv         }{Eigenvalues (lambda = "indicator") or singular values (lambda = "Burt", "adjusted" or "JCA") }
  \item{lambda     }{Scaling method}
  \item{inertia.e  }{Percentages of explained inertia}
  \item{inertia.t  }{Total inertia}
  \item{inertia.et }{Total percentage of explained inertia with the \code{nd}-dimensional solution}
  \item{levelnames }{Names of the factor/level combinations}
  \item{levels.n   }{Number of levels in each factor}
  \item{nd         }{User-specified dimensionality of the solution}
  \item{nd.max     }{Maximum possible dimensionality of the solution}
  \item{rownames   }{Row names}
  \item{rowmass    }{Row masses}
  \item{rowdist    }{Row chi-square distances to centroid}
  \item{rowinertia }{Row inertias}
  \item{rowcoord   }{Row standard coordinates}
  \item{colnames   }{Column names}
  \item{colmass    }{Column masses}
  \item{coldist    }{Column chi-square distances to centroid}
  \item{colinertia }{Column inertias}
  \item{colcoord   }{Column standard coordinates}
  \item{colsup     }{Indices of column supplementary points (of the Burt and Indicator matrix)}
  \item{Burt       }{Burt matrix}
  \item{Burt.upd   }{The updated Burt matrix (JCA only)}
  \item{subinertia }{Inertias of sub-matrices}
  \item{JCA.iter   }{Vector of length two containing the number of iterations and the epsilon (JCA only)}
  \item{call       }{Return of \code{match.call}}
      }

\seealso{\code{\link{eigen}}, \code{\link{plot.mjca}}, \code{\link{summary.mjca}}, \code{\link{print.mjca}} }
\examples{ 
library(MASS)
data(farms)
mjca(farms)

# Joint correspondence analysis:
mjca(farms, lambda = "JCA")

 }
\keyword{multivariate}

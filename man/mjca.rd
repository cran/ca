\name{mjca}
\alias{mjca}
\title{Multiple and joint correspondence analysis}
\description{Computation of multiple and joint correspondence analysis.}
\usage{mjca(obj, nd = 2, lambda = "adjusted", supcol = NA, subsetcol = NA, 
     ps = "", maxit = 50, epsilon = 0.0001)}
\arguments{
  \item{obj      }{A response pattern matrix (data frame containing factors).}
  \item{nd       }{Number of dimensions to be included in the output; if NA the maximum possible dimensions are included.}
  \item{lambda   }{Gives the scaling method. Possible values include \kbd{"indicator"}, \kbd{"Burt"}, \kbd{"adjusted"} and \kbd{"JCA"}.
                Using \kbd{lambda = "JCA"} results in a joint correspondence analysis using iterative adjusment of the Burt matrix in the solution space.}
  \item{supcol   }{Indices of supplementary columns.}
  \item{subsetcol}{Indices of subset categories.}
  \item{ps       }{Separator used for combining variable and category names.}
  \item{maxit    }{The maximum number of iterations (Joint Correspondence Analysis).}
  \item{epsilon  }{A convergence criterion (Joint Correspondence Analysis).}
          }
\details{The function \code{mjca} computes a multiple or joint correspondence analysis based on the eigenvalue decomposition of the Burt matrix.}
\value{
  \item{sv         }{Eigenvalues (\kbd{lambda = "indicator"}) or singular values (\kbd{lambda = "Burt"}, \kbd{"adjusted"} or \kbd{"JCA"}) }
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
  \item{subsetcol  }{Indices of subset columns}
  \item{Burt       }{Burt matrix}
  \item{Burt.upd   }{The updated Burt matrix (JCA only)}
  \item{subinertia }{Inertias of sub-matrices}
  \item{JCA.iter   }{Vector of length two containing the number of iterations and the epsilon (JCA only)}
  \item{call       }{Return of \code{match.call}}
      }

\references{Nenadic, O. and Greenacre, M. (2007), Correspondence analysis in R, with two- and three-dimensional graphics: The ca package. Journal of Statistical Software, 20 (3), available at \url{http://www.jstatsoft.org/v20/i03/}\cr
            Nenadic, O. and Greenacre, M. (2007), Computation of Multiple Correspondence Analysis, with Code in R, in Multiple Correspondence Analysis and Related Methods (eds. M. Greenacre and J. Blasius), Chapmann & Hall / CRC, Boca Raton, London, New York, pp. 523-551.\cr
            Greenacre, M.J. and Pardo, R. (2006), Subset correspondence analysis: visualizing relationships among a selected set of response categories from a questionnaire survey. Sociological Methods and Research, 35, pp. 193-218.}
\seealso{\code{\link{eigen}}, \code{\link{plot.mjca}}, \code{\link{summary.mjca}}, \code{\link{print.mjca}} }
\examples{ 
library(MASS)
data(farms)
mjca(farms)

# Joint correspondence analysis:
mjca(farms, lambda = "JCA")

 }
\keyword{multivariate}

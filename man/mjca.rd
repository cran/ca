\name{mjca}
\alias{mjca}
\alias{mjca.data.frame}
\alias{mjca.table}
\alias{mjca.array}
\alias{mjca.default}

\title{Multiple and joint correspondence analysis}
\description{Computation of multiple and joint correspondence analysis.}


\usage{

mjca(obj, ...)

\method{mjca}{data.frame}(obj, ...)
\method{mjca}{table}(obj, ...)
\method{mjca}{array}(obj, ...)

\method{mjca}{default}(obj, nd = 2, lambda = c("adjusted", "indicator", "Burt", "JCA"), 
     supcol = NA, subsetcat = NA, 
     ps = ":", maxit = 50, epsilon = 0.0001, reti = FALSE, ...)

}

\arguments{
  \item{obj      }{A response pattern matrix (data frame containing factors), or a frequency table (a \dQuote{table} object)
                   or an integer array.}
  \item{nd       }{Number of dimensions to be included in the output; if NA the maximum possible dimensions are included.}
  \item{lambda   }{Gives the scaling method. Possible values include \kbd{"indicator"}, \kbd{"Burt"}, \kbd{"adjusted"} and \kbd{"JCA"}.
                Using \kbd{lambda = "JCA"} results in a joint correspondence analysis using iterative adjusment of the Burt matrix in the solution space. See Details for descriptions of these options.}
  \item{supcol   }{Indices of supplementary columns.}
  \item{subsetcat}{Indices of subset categories (previously \kbd{subsetcol}).}
  \item{ps       }{Separator used for combining variable and category names.}
  \item{maxit    }{The maximum number of iterations (Joint Correspondence Analysis).}
  \item{epsilon  }{A convergence criterion (Joint Correspondence Analysis).}
  \item{reti     }{Logical indicating whether the indicator matrix should be included in the output.}
  \item{...      }{Arguments passed to \code{mjca.default}}
}

\details{
The function \code{mjca} computes a multiple or joint correspondence analysis based on the eigenvalue decomposition of the Burt matrix. The \code{lambda} option selects the scaling variant desired for
reporting inertias.

\itemize{
  \item \code{lambda="indicator"} gives multiple correspondence analysis based on the correspondence analysis of the indicator matrix, with corresponding inertias (eigenvalues).  
  \item \code{lambda="Burt"} gives the version of multiple correspondence analysis based on the correspondence analysis of the Burt matrix, the inertias of which are the squares of those for the indicator option.  
  \item \code{lambda="adjusted"} is the default option, giving improved percentages of inertia based on fitting the off-diagonal submatrices of the Burt matrix by rescaling the multiple correspondence analysis solution.  All these first three options give the same standard coordinates of the categories. 
  \item \code{lambda="JCA"} gives a joint correspondence analysis, which uses an iterative algorithm that optimally fits the off-diagonal submatrices of the Burt matrix.  The JCA solution does not have strictly nested  dimensions, so the percentage of inertia explained is given for the whole solution of chosen dimensionality, not for each dimension, but this percentage is optimal.
}

}

\value{
  \item{sv         }{Eigenvalues (\kbd{lambda = "indicator"}) or singular values (\kbd{lambda = "Burt"}, \kbd{"adjusted"} or \kbd{"JCA"}) }
  \item{lambda     }{Scaling method}
  \item{inertia.e  }{Percentages of explained inertia}
  \item{inertia.t  }{Total inertia}
  \item{inertia.et }{Total percentage of explained inertia with the \code{nd}-dimensional solution}
  \item{levelnames }{Names of the factor/level combinations, joined using \code{ps}}
  \item{factors    }{A matrix containing the names of the factors and the names of the factor levels}
  \item{levels.n   }{Number of levels in each factor}
  \item{nd         }{User-specified dimensionality of the solution}
  \item{nd.max     }{Maximum possible dimensionality of the solution}
  \item{rownames   }{Row names}
  \item{rowmass    }{Row masses}
  \item{rowdist    }{Row chi-square distances to centroid}
  \item{rowinertia }{Row inertias}
  \item{rowcoord   }{Row standard coordinates}
  \item{rowpcoord  }{Row principal coordinates}
  \item{rowctr     }{Row contributions}
  \item{rowcor     }{Row squared correlations}
  \item{colnames   }{Column names}
  \item{colmass    }{Column masses}
  \item{coldist    }{Column chi-square distances to centroid}
  \item{colinertia }{Column inertias}
  \item{colcoord   }{Column standard coordinates}
  \item{colpcoord  }{Column principal coordinates}
  \item{colctr     }{column contributions}
  \item{colcor     }{Column squared correlations}
  \item{colsup     }{Indices of column supplementary points (of the Burt and Indicator matrix)}
  \item{subsetcol  }{Indices of subset columns (\kbd{subsetcat})}
  \item{Burt       }{Burt matrix}
  \item{Burt.upd   }{The updated Burt matrix (JCA only)}
  \item{subinertia }{Inertias of sub-matrices}
  \item{JCA.iter   }{Vector of length two containing the number of iterations and the epsilon (JCA only)}
  \item{indmat     }{Indicator matrix if \code{reti} was set to \code{TRUE}}
  \item{call       }{Return of \code{match.call}}
      }

\references{Nenadic, O. and Greenacre, M. (2007), Correspondence analysis in R, with two- and three-dimensional graphics: The ca package. \emph{Journal of Statistical Software}, \bold{20 (3)}, \url{http://www.jstatsoft.org/v20/i03/}\cr
            Nenadic, O. and Greenacre, M. (2007), Computation of Multiple Correspondence Analysis, with Code in R, in \emph{Multiple Correspondence Analysis and Related Methods} (eds. M. Greenacre and J. Blasius), Boca Raton: Chapmann & Hall / CRC, pp. 523-551.\cr
            Greenacre, M.J. and Pardo, R. (2006), Subset correspondence analysis: visualizing relationships among a selected set of response categories from a questionnaire survey. \emph{Sociological Methods and Research}, \bold{35}, pp. 193-218.}

\seealso{\code{\link{eigen}}, \code{\link{plot.mjca}}, \code{\link{summary.mjca}}, \code{\link{print.mjca}} }

\examples{ 
data("wg93")
mjca(wg93[,1:4])

# table input
data(UCBAdmissions)
mjca(UCBAdmissions)
\dontrun{plot(mjca(UCBAdmissions))}

### Different approaches to multiple correspondence analysis:
# Multiple correspondence analysis based on the indicator matrix:
\dontrun{mjca(wg93[,1:4], lambda = "indicator")}

# Multiple correspondence analysis based on the Burt matrix:
\dontrun{mjca(wg93[,1:4], lambda = "Burt")}

# "Adjusted" multiple correspondence analysis (default setting):
\dontrun{mjca(wg93[,1:4], lambda = "adjusted")}

# Joint correspondence analysis:
\dontrun{mjca(wg93[,1:4], lambda = "JCA")}


### Subset analysis and supplementary variables:
# Subset analysis:
\dontrun{mjca(wg93[,1:4], subsetcat = (1:20)[-seq(3,18,5)])}

# Supplementary variables:
\dontrun{mjca(wg93, supcol = 5:7)}

 }
\keyword{multivariate}

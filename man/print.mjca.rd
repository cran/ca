\name{print.mjca}
\alias{print.mjca}
\title{Printing mjca objects}
\description{Printing method for multiple and joint correspondence analysis objects}
\usage{\method{print}{mjca}(x, ...) }
\arguments{
  \item{x}{Multiple or joint correspondence analysis object returned by \code{\link{mjca}}}
  \item{...}{Further arguments are ignored}
          }
\details{
The function \code{print.mjca} gives the basic statistics of the \code{mjca} object.  First the eigenvalues (that is, principal inertias) and their percentages with respect to total inertia are printed.   Then for the rows and columns respectively, the following are printed: the masses, chi-square distances of the points to the centroid (i.e., centroid of the active points), point inertias (for active points only) and principal coordinates on the first \code{nd} dimensions requested (default = 2 dimensions).  The function \code{\link{summary.mjca}} gives more detailed results about the inertia contributions of each point on each principal axis.\cr
For supplementary points, masses and inertias are not applicable.
}
\seealso{\code{\link{mjca}}}
\examples{
data("wg93")
print(mjca(wg93[,1:4]))
# equivalent to:
mjca(wg93[,1:4])
}

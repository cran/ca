\name{print.ca}
\alias{print.ca}
\title{Printing ca objects}
\description{Printing method for correspondence analysis objects}
\usage{\method{print}{ca}(x, ...) }
\arguments{
  \item{x}{Simple correspondence analysis object returned by \code{\link{ca}}}
  \item{...}{Further arguments are ignored}
          }
\details{
The function \code{print.ca} gives the basic statistics of the \code{ca} object.  First the eigenvalues (that is, principal inertias) and their percentages with respect to total inertia are printed.   Then for the rows and columns respectively, the following are printed: the masses, chi-square distances of the points to the centroid (i.e., centroid of the active points), point inertias (for active points only) and principal coordinates on the first \code{nd} dimensions requested (default = 2 dimensions).  The function \code{\link{summary.ca}} gives more detailed results about the inertia contributions of each point on each principal axis.\cr
For supplementary points, masses and inertias are not applicable.
}
\seealso{\code{\link{ca}}}
\examples{
data("smoke")
print(ca(smoke))
}

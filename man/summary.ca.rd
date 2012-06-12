\name{summary.ca}
\alias{summary.ca}
\title{Summarizing simple correspondence analysis}
\description{Textual output summarizing the results of \code{\link{ca}}, including a scree-plot of the principal inertias and row and column contributions.}
\usage{\method{summary}{ca}(object, scree = TRUE, ...)}
\arguments{
  \item{object}{Simple correspondence analysis object returned by \code{\link{ca}}.}
  \item{scree}{Logical flag specifying if a scree-plot should be included in the output.}
  \item{...}{Further arguments (ignored)}
          }
\details{
The function \code{summary.ca} gives the detailed numerical results of the \code{\link{ca}} function. All the eigenvalues (principal inertias) are listed, their percentages with respect to total inertia, and a bar chart (also known as a scree plot). Then for the set of rows and columns a table of results is given in a standard format, where quantities are either multiplied by 1000 or expressed in permills (thousandths): the mass of each point (x1000), the quality of display in the solution subspace of \code{nd} dimensions, the inertia of the point (in permills of the total inertia), and then for each dimension of the solution the principal coordinate (x1000), the (relative) contribution COR of the principal axis to the point inertia (x1000) and the (absolute) contribution CTR of the point to the inertia of the axis (in permills of the principal inertia). \cr
For supplementary points, masses, inertias and absolute contributions (CTR) are not applicable, but the relative contributions (COR) are valid as well as their sum over the set of chosen \code{nd} dimensions (QLT).
}
\examples{
data("smoke")
summary(ca(smoke))
}

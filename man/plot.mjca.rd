\name{plot.mjca}
\alias{plot.mjca}
\title{Plotting 2D maps in multiple and joint correspondence analysis}
\description{Graphical display of multiple and joint correspondence analysis results in two dimensions}
\usage{\method{plot}{mjca}(x, dim = c(1,2), map = "symmetric", centroids = FALSE, what = c("all", "all"), 
               mass = c(FALSE, FALSE), contrib = c("none", "none"), 
               col = c("#000000", "#FF0000"), pch = c(16, 1, 17, 24), 
               labels = c(2, 2), arrows = c(FALSE, FALSE), ...) }
\arguments{
  \item{x}{Multiple or joint correspondence analysis object returned by \code{\link{mjca}}}
  \item{dim}{Numerical vector of length 2 indicating the dimensions to plot on horizontal and vertical axes respectively; default is first dimension horizontal and second dimension vertical.}
  \item{map}{Character string specifying the map type. Allowed options include \cr
              \code{"symmetric"} (default) \cr
              \code{"rowprincipal"} \cr
              \code{"colprincipal"} \cr
              \code{"symbiplot"} \cr
              \code{"rowgab"} \cr
              \code{"colgab"} \cr
              \code{"rowgreen"} \cr
              \code{"colgreen"}
            }
  \item{centroids}{Logical indicating if column centroids should be added to the plot}
  \item{what}{Vector of two character strings specifying the contents of the plot. First entry sets the rows and the second entry the columns. Allowed values are \cr
              \code{"all"} (all available points, default) \cr
              \code{"active"} (only active points are displayed) \cr
              \code{"passive"} (only supplementary points are displayed) \cr
              \code{"none"} (no points are displayed) \cr
              The status (active or supplementary) of rows and columns is set in \code{\link{mjca}} using the options \code{suprow} and \code{supcol}.}
  \item{mass}{Vector of two logicals specifying if the mass should be represented by the area of the point symbols (first entry for rows, second one for columns)}
  \item{contrib}{Vector of two character strings specifying if contributions (relative or absolute) should be represented by different colour intensities. Available options are\cr
                 \code{"none"} (contributions are not indicated in the plot).\cr
                 \code{"absolute"} (absolute contributions are indicated by colour intensities).\cr
                 \code{"relative"} (relative contributions are indicated by colour intensities).\cr
                 If set to \code{"absolute"} or \code{"relative"}, points with zero contribution are displayed in white. The higher the contribution of a point, the closer the corresponding colour to the one specified by the \code{col} option.}
  \item{col}{Vector of length 2 specifying the colours of row and column point symbols, by default black for rows and red for columns. Colours can be entered in hexadecimal (e.g. \code{"#FF0000"}), rgb (e.g. \code{rgb(1,0,0)}) values or by R-name (e.g. \code{"red"}). }
  \item{pch}{Vector of length 4 giving the type of points to be used for row active and supplementary, column active and supplementary points. See \code{\link{pchlist}} for a list of symbols.}
  \item{labels}{Vector of length two specifying if the plot should contain symbols only (\code{0}), labels only (\code{1}) or both symbols and labels (\code{2}). Setting \code{labels} to \code{2} results in the symbols being plotted at the coordinates and the labels with an offset.}
  \item{arrows}{Vector of two logicals specifying if the plot should contain points (FALSE, default) or arrows (TRUE). First value sets the rows and the second value sets the columns.}
  \item{...}{Further arguments passed to \code{\link{plot}} and \code{\link{points}}.}
          }
\details{
The function \code{plot.mjca} makes a two-dimensional map of the object created by \code{mjca} with respect to two selected dimensions.  By default the scaling option of the map is \code{"symmetric"}, that is the so-called \emph{symmetric map}. In this map both the row and column points are scaled to have inertias (weighted variances) equal to the principal inertia (eigenvalue) along the principal axes, that is both rows and columns are in pricipal coordinates. Other options are as follows:  

\item{-}{\code{"rowprincipal"} or \code{"colprincipal"} - these are the so-called \emph{asymmetric maps}, with either rows in principal coordinates and columns in standard coordinates, or vice versa (also known as row-metric-preserving or column-metric-preserving respectively). These maps are biplots;}

\item{-}{\code{"symbiplot"} - this scales both rows and columns to have variances equal to the singular values (square roots of eigenvalues), which gives a symmetric biplot but does not preserve row or column metrics;}

\item{-}{\code{"rowgab"} or \code{"colgab"} - these are asymmetric maps (see above) with rows (respectively, columns) in principal coordinates and columns (respectively, rows) in standard coordinates multiplied by the mass of the corresponding point. These are also biplots and were proposed by Gabriel & Odoroff (1990);}

\item{-}{\code{"rowgreen"} or \code{"colgreen"} - these are similar to \code{"rowgab"} and \code{"colgab"} except that the points in standard coordinates are multiplied by the square root of the corresponding masses, giving reconstructions of the standardized residuals.}

This function has options for sizing and shading the points.  If the option \code{mass} is \code{TRUE} for a set of points, the size of the point symbol is proportional to the relative frequency (mass) of each point.  If the option \code{contrib} is \code{"absolute"} or \code{"relative"} for a set of points, the colour intensity of the point symbol is proportional  to the absolute contribution of the points to the planar display or, respectively, the quality of representation of the points in the display.

}

\references{
Gabriel, K.R. and Odoroff, C. (1990). Biplots in biomedical research. \emph{Statistics in Medicine}, 9, pp. 469-485. \cr
Greenacre, M.J. (1993) \emph{Correspondence Analysis in Practice}.  Academic Press, London. \cr
Greenacre, M.J. (1993) Biplots in correspondence Analysis, \emph{Journal of Applied Statistics}, 20, pp. 251 - 269.
}
\seealso{\code{\link{mjca}}, \code{\link{summary.mjca}}, \code{\link{print.mjca}}, \code{\link{pchlist}}}
\examples{
library(MASS)
data(farms)

# A two-dimensional map with standard settings
plot(mjca(farms))

# Mass for columns represented by the size of the point symbols
plot(mjca(farms), mass = c(FALSE, TRUE))

}
\keyword{}

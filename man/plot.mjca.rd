\name{plot.mjca}
\alias{plot.mjca}
\title{Plotting 2D maps in multiple and joint correspondence analysis}
\description{Graphical display of multiple and joint correspondence analysis results in two dimensions}
\usage{\method{plot}{mjca}(x, dim = c(1,2), map = "symmetric", centroids = FALSE, 
     what = c("none", "all"), mass = c(FALSE, FALSE), 
     contrib = c("none", "none"), col = c("#000000", "#FF0000"), 
     pch = c(16, 1, 17, 24), labels = c(2, 2), 
     arrows = c(FALSE, FALSE), xlab = "_auto_", ylab = "_auto_", ...) }
\arguments{
  \item{x}{Multiple or joint correspondence analysis object returned by \code{\link{mjca}}}
  \item{dim}{Numerical vector of length 2 indicating the dimensions to plot on horizontal and vertical axes respectively; default is first dimension horizontal and second dimension vertical.}
  \item{map}{Character string specifying the map type. Allowed options include \cr
              \kbd{"symmetric"} (default) \cr
              \kbd{"rowprincipal"} \cr
              \kbd{"colprincipal"} \cr
              \kbd{"symbiplot"} \cr
              \kbd{"rowgab"} \cr
              \kbd{"colgab"} \cr
              \kbd{"rowgreen"} \cr
              \kbd{"colgreen"}
            }
  \item{centroids}{Logical indicating if column centroids should be added to the plot}
  \item{what}{Vector of two character strings specifying the contents of the plot. First entry sets the rows and the second entry the columns. Allowed values are \cr
              \kbd{"all"} (all available points, default) \cr
              \kbd{"active"} (only active points are displayed) \cr
              \kbd{"passive"} (only supplementary points are displayed) \cr
              \kbd{"none"} (no points are displayed) \cr
              The status (active or supplementary) of columns is set in \code{\link{mjca}} using the option \code{supcol}.}
  \item{mass}{Vector of two logicals specifying if the mass should be represented by the area of the point symbols (first entry for rows, second one for columns)}
  \item{contrib}{Vector of two character strings specifying if contributions (relative or absolute) should be represented by different colour intensities. Available options are\cr
                 \kbd{"none"} (contributions are not indicated in the plot).\cr
                 \kbd{"absolute"} (absolute contributions are indicated by colour intensities).\cr
                 \kbd{"relative"} (relative contributions are indicated by colour intensities).\cr
                 If set to \kbd{"absolute"} or \kbd{"relative"}, points with zero contribution are displayed in white. The higher the contribution of a point, the closer the corresponding colour to the one specified by the \code{col} option.}
  \item{col}{Vector of length 2 specifying the colours of row and column point symbols, by default black for rows and red for columns. Colours can be entered in hexadecimal (e.g. \kbd{"#FF0000"}), rgb (e.g. \kbd{rgb(1,0,0)}) values or by R-name (e.g. \kbd{"red"}). }
  \item{pch}{Vector of length 4 giving the type of points to be used for row active and supplementary, column active and supplementary points. See \code{\link{pchlist}} for a list of symbols.}
  \item{labels}{Vector of length two specifying if the plot should contain symbols only (\kbd{0}), labels only (\kbd{1}) or both symbols and labels (\kbd{2}). Setting \code{labels} to \kbd{2} results in the symbols being plotted at the coordinates and the labels with an offset.}
  \item{arrows}{Vector of two logicals specifying if the plot should contain points (\kbd{FALSE}, default) or arrows (\kbd{TRUE}). First value sets the rows and the second value sets the columns.}
  \item{xlab, ylab}{Labels for horizontal and vertical axes.  The default, \code{"_auto_"} means that the function auto-generates a label of the form \code{Dimension X (xx.xx \%}}
  \item{...}{Further arguments passed to \code{\link{plot}} and \code{\link{points}}.}
          }
\details{
The function \code{plot.mjca} makes a two-dimensional map of the object created by \code{mjca} with respect to two selected dimensions.  By default the scaling option of the map is \kbd{"symmetric"}, that is the so-called \emph{symmetric map}. In this map both the row and column points are scaled to have inertias (weighted variances) equal to the principal inertia (eigenvalue) along the principal axes, that is both rows and columns are in pricipal coordinates. Other options are as follows:
\itemize{  
\item{-}{\kbd{"rowprincipal"} or \kbd{"colprincipal"} - these are the so-called \emph{asymmetric maps}, with either rows in principal coordinates and columns in standard coordinates, or vice versa (also known as row-metric-preserving or column-metric-preserving respectively). These maps are biplots;}
}
\itemize{
\item{-}{\kbd{"symbiplot"} - this scales both rows and columns to have variances equal to the singular values (square roots of eigenvalues), which gives a symmetric biplot but does not preserve row or column metrics;}
}
\itemize{
\item{-}{\kbd{"rowgab"} or \kbd{"colgab"} - these are asymmetric maps (see above) with rows (respectively, columns) in principal coordinates and columns (respectively, rows) in standard coordinates multiplied by the mass of the corresponding point. These are also biplots and were proposed by Gabriel & Odoroff (1990);}
}
\itemize{
\item{-}{\kbd{"rowgreen"} or \kbd{"colgreen"} - these are similar to \kbd{"rowgab"} and \kbd{"colgab"} except that the points in standard coordinates are multiplied by the square root of the corresponding masses, giving reconstructions of the standardized residuals.}
}
This function has options for sizing and shading the points.  If the option \code{mass} is \kbd{TRUE} for a set of points, the size of the point symbol is proportional to the relative frequency (mass) of each point.  If the option \code{contrib} is \kbd{"absolute"} or \kbd{"relative"} for a set of points, the colour intensity of the point symbol is proportional  to the absolute contribution of the points to the planar display or, respectively, the quality of representation of the points in the display.
To globally resize all the points (and text labels), use \code{par("cex"=)} before the plot.
}

\value{
In addition to the side effect of producing the plot, the function invisibly returns the coordinates
of the plotted points, a list of two components, with names \code{rows} and \code{cols}.
These can be used to further annotate the plot using base R plotting functions.
}

\references{
Gabriel, K.R. and Odoroff, C. (1990). Biplots in biomedical research. \emph{Statistics in Medicine}, \bold{9}, pp. 469-485. \cr
Greenacre, M.J. (1993) \emph{Correspondence Analysis in Practice}. London: Academic Press. \cr
Greenacre, M.J. (1993) Biplots in correspondence Analysis, \emph{Journal of Applied Statistics}, \bold{20}, pp. 251 - 269.
}
\seealso{\code{\link{mjca}}, \code{\link{summary.mjca}}, \code{\link{print.mjca}}, \code{\link{pchlist}}}
\examples{
data("wg93")

# A two-dimensional map with standard settings
plot(mjca(wg93[,1:4]))

}


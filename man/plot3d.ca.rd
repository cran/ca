\name{plot3d.ca}
\alias{plot3d.ca}
\title{Plotting 3D maps in correspondence analysis}
\description{Graphical display of correspondence analysis in three dimensions}
\usage{plot3d.ca(x, dim = c(1, 2, 3), map = "symmetric", what = c("all", "all"), 
             contrib = c("none", "none"), col = c("#6666FF","#FF6666"), 
             labcol  = c("#0000FF", "#FF0000"), pch = c(16, 1, 18, 9), 
             labels = c(2, 2), sf = 0.00002, arrows  = c(FALSE, FALSE), ...) }
\arguments{
  \item{x}{Simple correspondence analysis object returned by ca}
  \item{dim}{Numerical vector of length 2 indicating the dimensions to plot}
  \item{map}{Character string specifying the map type. Allowed options include \cr 
              "symmetric" (default) \cr
              "rowprincipal" \cr
              "colprincipal" \cr
              "symbiplot" \cr 
              "rowgab" \cr
              "colgab" \cr
              "rowgreen" \cr
              "colgreen"
            }
  \item{what}{Vector of two character strings specifying the contents of the plot. First entry sets the rows and the second entry the columns. Allowed values are \cr
              "none" (no points are displayed) \cr
              "active" (only active points are displayed, default) \cr
              "supplementary" (only supplementary points are displayed) \cr
              "all" (all available points) \cr
              The status (active or supplementary) is set in ca() using the tl-option.}
  \item{contrib}{Vector of two character strings specifying if contributions (relative or absolute) should be indicated by different colour intensities. Available options are\cr
                 "none" (contributions are not indicated in the plot).\cr
                 "absolute" (absolute contributions are indicated by colour intensities).\cr
                 "relative" (relative conrributions are indicated by colour intensities).\cr
                 If set to "absolute" or "relative", points with zero contribution are displayed in white. The higher the contribution of a point, the closer the corresponding colour to the one specified by the col option.}
  \item{col}{Vector of length 2 specifying the colours of row and column profiles. Colours can be entered in hexadecimal (e.g. "\#FF0000"), rgb (e.g. rgb(1,0,0)) values or by R-name (e.g. "red"). }
  \item{labcol}{Vector of length 2 specifying the colours of row and column labels. }
  \item{pch}{Vector of length 2 giving the type of points to be used for rows and columns.}
  \item{labels}{Vector of length two specifying if the plot should contain symbols only (\code{0}), labels only (\code{1}) or both symbols and labels (\code{2}). Setting \code{labels} to \code{2} results in the symbols being plotted at the coordinates and the labels with an offset.}
  \item{sf}{A scaling factor for the volume of the 3d primitives.}
  \item{arrows}{Vector of two logicals specifying if the plot should contain points (FALSE, default) or arrows (TRUE). First value sets the rows and the second value sets the columns.}
  \item{...}{Further arguments passed to the rgl functions.}
          }
\details{ }

\references{}
\seealso{ca}
\examples{}
\keyword{}

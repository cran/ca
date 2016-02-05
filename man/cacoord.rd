\name{cacoord}
\alias{cacoord}
\title{Extracting coordinates from ca and mjca objects.}
\description{Extracting standard and principal coordinates as well as various row and column scaling configurations for visual display from \kbd{ca} and \kbd{mjca} objects.}
\usage{cacoord(obj, 
               type = c("standard", "principal", 
                        "symmetric", "rowprincipal", "colprincipal", "symbiplot", 
                        "rowgab", "colgab", "rowgreen", "colgreen"), 
               dim  = NA,
               rows = NA,
               cols = NA,
               ...)}
\arguments{
  \item{obj }{A \kbd{ca} or \kbd{mjca} object returned by \code{\link{ca}} or \code{\link{mjca}}.}
  \item{type}{The type of coordinates to extract (\kbd{"standard"} or \kbd{"principal"}). The remaining options (\kbd{"symmetric"}, ..., \kbd{"colgreen"}) return the corresponding row/column coordinate configuration for the map scaling options described in \code{\link{plot.ca}} where the corresponding argument is \kbd{map}.}
  \item{dim }{The dimensions to return. If \kbd{NA}, all available dimensions are returned.}
  \item{rows}{Logical indicating whether to return the row coordinates (see below for details).}
  \item{cols}{Logical indicating whether to return the column coordinates (see below for details).}
  \item{... }{Further arguments (ignored).}
  }
\details{The function \code{cacoord} returns the standard or principal coordinates of a CA or MCA solution. Additionally, row and column scaling configurations for plotting methods can be computed (see \code{\link{plot.ca}} for details).\cr 
Note that by default row and column coordinates are computed (i.e. for \kbd{(rows=NA&cols=NA)|(rows=TRUE&cols=TRUE)}). Using \kbd{rows=TRUE} (and \kbd{cols=NA} or \kbd{cols=FALSE}) returns a matrix with the row coordinates, and for \kbd{cols=TRUE} (and \kbd{cols=NA} or \kbd{cols=FALSE}) a matrix with the column coordinates is returned.}
\value{A list with the slots \kbd{rows} (row coordinates) and \kbd{columns} (column coordinates). When computing only row or only column coordinates, a matrix (with the corresponding row or column coordinates) is returned.}
\seealso{\code{\link{ca}},\code{\link{mjca}},\code{\link{plot.ca}},\code{\link{plot.mjca}}}
\keyword{multivariate}

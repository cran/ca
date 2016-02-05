\name{multilines}
\alias{multilines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw lines for groups distinguished by a factor
}
\description{
This is a convenience function for drawing a set of lines
distinguished by the levels of a factor.
It can be used to make more attractive plots than
available via \code{\link{plot.mjca}}.

}
\usage{
multilines(XY, group=NULL, which=1:nf, sort=1, type='l', col=palette(), lwd=1, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XY}{
A two-column data frame or matrix
}
  \item{group}{
A factor; a separate line is drawn for each level included in \code{which}
}
  \item{which}{
An integer vector used to select the factors for which lines are drawn.
By default, all lines are drawn.
}

  \item{sort}{
Column of \code{XY} to sort upon before drawing the line for each group
}
  \item{type}{
Line type: \code{"l"} for line, \code{"b"} for line and points
}
  \item{col}{
A vector of colors to be used for the various lines, in the order of the levels
in \code{group}; recycled as necessary.
}
  \item{lwd}{
A vector of line widths to be used for the various lines; recycled as necessary
}
  \item{\dots}{
Other graphic parameters passed to \code{\link[graphics]{lines}}, e.g., \code{lty}
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
none
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
Michael Friendly
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[graphics]{lines}}
}
\examples{
if (require(vcd)) {
  data(PreSex, package="vcd")
  presex.mca <- mjca(PreSex)
  res <- plot(presex.mca, labels=0, pch='.', cex.lab=1.2)
  coords <- data.frame(res$cols, presex.mca$factors)      
  nlev <- rle(as.character(coords$factor))$lengths
  fact <- unique(as.character(coords$factor))
  
  cols <- c("blue", "red", "brown", "black")
  lwd <- c(2, 2, 2, 4)
  
  plot(Dim2 ~ Dim1, type='n', data=coords)
  points(coords[,1:2], pch=rep(16:19, nlev), col=rep(cols, nlev), cex=1.2)
  text(coords[,1:2], labels=coords$level, col=rep(cols, nlev), pos=3, cex=1.2, xpd=TRUE)
  
  multilines(coords[, c("Dim1", "Dim2")], group=coords$factor, col=cols, lwd=lwd)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{aplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

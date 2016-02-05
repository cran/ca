\name{caconv}
\alias{caconv}
\title{Converting data types in CA and MCA}
\description{Conversion from and to a number of different data types commonly used in CA and MCA (frequency tables, response pattern matrices, indicator matrices and Burt matrices).}
\usage{caconv(x, from = c("freq", "rpm", "ind", "Burt"), to = c("rpm", "ind", "Burt", "freq"), 
              nlev = NA, vars = c(1,2), ...)}
\arguments{
  \item{x   }{A matrix (two-way frequency table, indicator matrix, or Burt matrix) or data frame (response pattern matrix).}
  \item{from}{The type of input data in \kbd{x}: a frequency table (\kbd{"freq"}), or a response pattern matrix (\kbd{"rpm"}), or an indicator matrix (\kbd{"ind"}), or a Burt matrix (\kbd{"Burt"}).}
  \item{to  }{The data type into which \kbd{x} should be converted.}
  \item{nlev}{A vector containing the number of levels for each categorical variable (for \kbd{from="ind"} or \kbd{from="Burt"}). If \kbd{NA}, \kbd{nlev} is computed from the data.}
  \item{vars}{A vector of length 2 specifying the index of the variables to use for converting to \kbd{"freq"} (i.e. to a regular two-way frequency table).}
  \item{... }{Further arguments (ignored).}
          }
\details{The function \code{caconv} converts between data types in CA and MCA. Note that a conversion from \kbd{from="Burt"} to \kbd{to="ind"} or \kbd{to="rpm"} is not supported.}
\value{A matrix or data frame containing the converted data (with the type specified in \kbd{to}).}
\seealso{\code{\link{ca}},\code{\link{mjca}}}
\keyword{multivariate}

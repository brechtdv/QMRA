\name{giardia}
\alias{giardia}
\docType{data}
\title{Giardia dataset}
\description{
The Giardia dataset contains concentrations of Giardia oocysts and left-censoring indicators.
}
\usage{data(giardia)}
\format{
  A data frame with 26 observations on the following 2 variables.
  \describe{
    \item{\code{x}}{observed concentrations per 100mL; a numeric vector}
    \item{\code{d}}{censoring indicator (1=left-censored; 0=observed); a numeric vector}
  }
}

\references{
\itemize{
  \item{
    Haas CN, Rose JB, Gerba CP (1999)
    \emph{Quantitative Microbial Risk Assessment.}
	John Wiley & Sons, Inc.
	}
}}

\seealso{
  \code{\link{ea_conc}}, \code{\link{bea_conc}}
}
\examples{
data(giardia)

## Exposure assessment from concentration data
ea_conc(x = x, d = d, data = giardia, model = "weibull")
}
\keyword{datasets}

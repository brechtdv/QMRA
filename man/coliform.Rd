\name{coliform}
\alias{coliform}
\docType{data}
\title{MPN coliform dilution assay dataset}
\description{
Dilution assays report the presence of absence of the hazard in a given quantity. The \code{coliform} dataset contains dilution assay results and the corresponding volumes and replicates.
}
\usage{data(coliform)}
\format{
  A data frame with 3 observations on the following 3 variables.
  \describe{
    \item{\code{volume}}{a numeric vector}
    \item{\code{replicates}}{a numeric vector}
    \item{\code{positive}}{a numeric vector}
  }
}

\source{
\itemize{
  \item{
    Haas CN, Rose JB, Gerba CP (1999)
    \emph{Quantitative Microbial Risk Assessment.}
	John Wiley & Sons, Inc.
	}
}}

\seealso{
  \code{\link{ea_presence}}, \code{\link{bea_presence}}
}

\examples{
data(coliform)

## Exposure assessment from presence/absence data
ea_presence(x = positive, q = volume, rep = replicates, data = coliform)
}
\keyword{datasets}

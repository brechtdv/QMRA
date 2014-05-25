\name{crypto}
\alias{crypto}
\docType{data}
\title{Cryptosporidium dataset}
\description{
The Cryptosporidium dataset contains oocyst counts and corresponding volumes from 52 water samples.
}
\usage{data(crypto)}
\format{
  A data frame with 52 observations on the following 2 variables.
  \describe{
    \item{\code{volume}}{sample volumes; a numeric vector}
    \item{\code{counts}}{Cryptosporidium oocyst counts; a numeric vector}
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

\examples{
data(crypto)

## Exposure assessment from count data
ea_count(x = counts, q = volume, model = "negbin", data = crypto)
}
\keyword{datasets}

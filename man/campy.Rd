\name{campy}
\alias{campy}
\docType{data}
\title{Campylobacter dose-response dataset}
\description{
The Campylobacter dataset contains dose-response data of 72 individuals enrolled in a feeding trial.
}
\usage{data(campy)}
\format{
  A data frame with 7 observations on the following 3 variables.
  \describe{
    \item{\code{dose}}{administered dose; a numeric vector}
    \item{\code{total}}{total number of exposed individuals; a numeric vector}
    \item{\code{infected}}{number of infected individuals; a numeric vector}
  }
}

\source{
\itemize{
  \item{
    Haas CN, Rose JB, Gerba CP (1999)
    \emph{Quantitative Microbial Risk Assessment.}
	John Wiley & Sons, Inc.
	}
  \item{
    Black RE, Levine MM, Clements ML, Hughes TP, Blaser MJ (1988). Experimental \emph{Campylobacter jejuni} infection in humans. \emph{Journal of Infectious Diseases}, 157(3), 472-479.
  }
}}

\seealso{
  \code{\link{drm}}, \code{\link{bdrm}}
}

\examples{
data(campy)

## Fit Beta-Poisson dose-response model to Campylobacter dataset
drm(x = infected, n = total, dose = dose, data = campy, model = "bp")
}
\keyword{datasets}

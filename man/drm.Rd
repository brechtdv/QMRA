\name{drm}
\alias{drm}

\title{Dose-response modelling}

\description{
Parametric modelling of dose-response data.
}

\usage{
drm(x, n, dose, data,
    model = c("betapoisson", "bp", 
              "exponential", "exp", 
              "loglogistic", "ll",
              "logprobit", "lp",
              "extremevalue", "ev"),
    ...)
}

\arguments{
  \item{x}{Vector of infected individuals.}
  \item{n}{Vector of exposed individuals.}
  \item{dose}{Vector of administered doses.}
  \item{data}{An optional data frame, containing the variables in the model.
              If not found in \code{data}, the variables are taken from the
              environment from which \code{drm} is called.}
  \item{model}{A character string naming the dose-response model. See details below.}
  \item{\dots}{Arguments to be passed to \code{\link{mle}}.}
}

\details{
  Available models:

\itemize{
  \item{\strong{Beta-Poisson:  }}{\code{model = "betapoisson"} or \code{"bp"}}
  \item{\strong{Exponential:  }}{\code{model = "exponential"} or \code{"exp"}}
  \item{\strong{Log-Logistic:  }}{\code{model = "loglogistic"} or \code{"ll"}}
  \item{\strong{Log-Probit:  }}{\code{model = "logprobit"} or \code{"lp"}}
  \item{\strong{Extreme Value:  }}{\code{model = "extremevalue"} or \code{"ev"}}
}}

\value{An object of class \code{"\linkS4class{drm}"}.}

\author{\email{Brecht.Devleesschauwer@UGent.be}}

\references{
\itemize{
  \item{
    Haas CN, Rose JB, Gerba CP (1999)
    \emph{Quantitative Microbial Risk Assessment.}
	John Wiley & Sons, Inc.
	}
}}

\seealso{
  \code{\link{bdrm}}, for Bayesian dose-response modelling\cr
}

\examples{
## Fit Beta-Poisson dose-response model to Campylobacter dataset
drm(x = infected, n = total, dose = dose, data = campy, model = "bp")
}

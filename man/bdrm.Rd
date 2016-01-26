\name{bdrm}
\alias{bdrm}

\title{Bayesian parametric dose-response modelling}

\description{
Bayesian parametric modelling of dose-response data.
}

\usage{
bdrm(x, n, dose, data,
     model = c("betapoisson", "bp", 
               "exponential", "exp", 
               "loglogistic", "ll",
               "logprobit", "lp",
               "extremevalue", "ev"),
     inits = NULL, nchains = 2, burnin = 10000, update = 10000,
     verbose = FALSE)
}

\arguments{
  \item{x}{Vector of infected individuals.}
  \item{n}{Vector of exposed individuals.}
  \item{dose}{Vector of administered doses.}
  \item{data}{An optional data frame, containing the variables in the model.
              If not found in \code{data}, the variables are taken from the
              environment from which \code{bdrm} is called.}
  \item{model}{A character string naming the dose-response model. See details below.}
  \item{inits}{Named list of initial values; defaults to \code{NULL}, leading \code{JAGS} to generate inits automatically (see \code{\link{jags.model}}).}
  \item{nchains}{Number of model chains, should be \eqn{\ge 2};
                 defaults to 2.}
  \item{burnin}{Number of samples to discard as burn-in; defaults to 10000.}
  \item{update}{Number of samples to retain; defaults to 10000.}
  \item{verbose}{Should JAGS process info be printed to the R console?
                 defaults to \code{FALSE}.}
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

\author{\email{brechtdv@gmail.com}}

\references{
\itemize{
  \item{
    Haas CN, Rose JB, Gerba CP (1999)
    \emph{Quantitative Microbial Risk Assessment.}
  John Wiley & Sons, Inc.
	}
}}

\seealso{
  \code{\link{drm}}, for frequentist dose-response modelling\cr
}

\examples{
## Fit Beta-Poisson dose-response model to Campylobacter dataset
bdrm(x = infected, n = total, dose = dose,
     data = campy, model = "bp",
     inits = list(alpha = 1, beta = 1))
}

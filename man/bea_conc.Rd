\name{bea_conc}
\alias{bea_conc}

\title{Bayesian exposure assessment from concentration data}

\description{
Bayesian parametric modelling of exposure concentration from
(possibly left-censored) concentration data.
}

\usage{
bea_conc(x, d, data,
         model = c("gamma", "g",
                   "lognorm", "ln",
                   "weibull", "w",
                   "invgauss", "ig"),
         inits = NULL, nchains = 2, burnin = 5000, update = 5000,
         verbose = FALSE)
}

\arguments{
  \item{x}{The observed concentrations.}
  \item{d}{A vector indicating if the observations were left-censored (1) or
           not (0).}
  \item{data}{An optional data frame, containing the variables in the model.
              If not found in \code{data}, the variables are taken from the
              environment from which \code{bea_conc} is called.}
  \item{model}{A character string naming the concentration model. See details below.}
  \item{inits}{Named list of initial values; defaults to \code{NULL}, leading \code{JAGS} to generate inits automatically (see \code{\link{jags.model}}).}
  \item{nchains}{Number of model chains, should be \eqn{\ge 2};
                 defaults to 2.}
  \item{burnin}{Number of samples to discard as burn-in; defaults to 5000.}
  \item{update}{Number of samples to retain; defaults to 5000.}
  \item{verbose}{Should JAGS process info be printed to the R console?
                 defaults to \code{FALSE}.}
}

\details{
  Distributions available:

\itemize{
  \item{\strong{Gamma:  }}{\code{model = "gamma"} or \code{"g"}}
  \item{\strong{Log-Normal:  }}{\code{model = "lognorm"} or \code{"ln"}}
  \item{\strong{Weibull:  }}{\code{model = "weibull"} or \code{"w"}}
  \item{\strong{Inverse Gaussian:  }}{\code{model = "invgauss"} or \code{"ig"}}
}}

\value{An object of class \code{"\linkS4class{bea}"}.}

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
  \code{\link{bea_count}}, for modelling exposure from count data\cr
  \code{\link{bea_presence}}, for modelling exposure from presence/absence data

}

\examples{
bea_conc(x = x, d = d, data = giardia, model = "weibull")
}

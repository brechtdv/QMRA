\name{bea_count}
\alias{bea_count}

\title{Bayesian exposure assessment from count data}

\description{
Bayesian parametric modelling of exposure concentration from count data.
}

\usage{
bea_count(x, q = 1, data,
          model = c("poisson", "p",
                    "negbin", "nb",
                    "poislognorm", "pln",
                    "poisweibull", "pw"),
          inits = NULL, nchains = 2, burnin = 5000, update = 5000,
          verbose = FALSE)
}

\arguments{
  \item{x}{The observed counts.}
  \item{q}{The quantities (e.g., volumes, masses) in which the counts were
           observed; defaults to 1.}
  \item{data}{An optional data frame, containing the variables in the model.
              If not found in \code{data}, the variables are taken from the
              environment from which \code{ea_count} is called.}
  \item{model}{A character string naming the count model to be fitted. See details below.}
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
  \item{\strong{Poisson:  }}{\code{model = "poisson"} or \code{"p"}}
  \item{\strong{Negative Binomial:  }}{\code{model = "negbin"} or \code{"nb"}}
  \item{\strong{Poisson-LogNormal:  }}{\code{model = "poislognorm"} or \code{"pln"}}
  \item{\strong{Poisson-Weibull:  }}{\code{model = "poisweibull"} or \code{"pw"}}
}}

\value{An object of class \code{"\linkS4class{bea}"}.}

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
  \code{\link{bea_presence}}, for modelling exposure from presence/absence data\cr
  \code{\link{bea_conc}}, for modelling exposure from concentration data
}

\examples{
## Bayesian exposure assessment from count data
bea_count(x = counts, q = volume, model = "negbin", data = crypto,
          inits = list(shape = 1, rate = 1))
}

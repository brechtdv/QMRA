\name{bea_presence}
\alias{bea_presence}

\title{Bayesian exposure assessment from presence/absence data}

\description{
Bayesian parametric modelling of exposure concentration from
presence/absence data.
}

\usage{
bea_presence(x, q = 1, replicates = rep(1, length(x)), data,
             model = c("poisson"),
             nchains = 2, burnin = 5000, update = 5000, verbose = FALSE)
}

\arguments{
  \item{x}{A vector indicating the number of positive samples per examined
           quantity.}
  \item{q}{The quantities (e.g., volumes, masses) in which presence/absence
           was observed; defaults to 1.}
  \item{replicates}{The number of replicates of each quantity \code{q};
                    defaults to \code{rep(1, length(x))}.}
  \item{data}{An optional data frame, containing the variables in the model.
              If not found in \code{data}, the variables are taken from the
              environment from which \code{ea_presence} is called.}
  \item{model}{A character string naming the model to be fitted;
                currently only \code{"poisson"} is supported.
                See details below.}
  \item{nchains}{Number of model chains, should be \eqn{\ge 2};
                 defaults to 2.}
  \item{burnin}{Number of samples to discard as burn-in; defaults to 5000.}
  \item{update}{Number of samples to retain; defaults to 5000.}
  \item{verbose}{Should JAGS process info be printed to the R console?
                 defaults to \code{FALSE}.}
}

\details{
  Available distributions:

\itemize{
  \item{\strong{Poisson:  }}{\code{model = "poisson"}}
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
  \code{\link{bea_count}}, for modelling exposure from count data\cr
  \code{\link{bea_conc}}, for modelling exposure from concentration data
}

\examples{
## Bayesian exposure assessment from presence/absence data
bea_presence(x = positive, q = volume, rep = replicates, data = coliform)
}

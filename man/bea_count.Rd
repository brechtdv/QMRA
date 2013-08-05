\name{bea_count}
\alias{bea_count}

\title{Bayesian exposure assessment from count data}

\description{
Bayesian parametric modelling of exposure concentration from count data.
}

\usage{
bea_count(x, q = 1, data,
          family = c("poisson", "negbin", "poislognorm", "poisinvgauss", "poisgeninvgauss"),
          nchains = 2, burnin = 1000, update = 5000, verbose = FALSE)
}

\arguments{
  \item{x}{The observed counts.}
  \item{q}{The quantities (e.g., volumes, masses) in which the counts were
           observed; defaults to 1.}
  \item{data}{An optional data frame, containing the variables in the model.
              If not found in \code{data}, the variables are taken from the
              environment from which \code{ea_count} is called.}
  \item{family}{A character string naming the family function;
                should be one of \code{"poisson"}, \code{"negbin"},
                \code{"poislognorm"}, \code{"poisinvgauss"},
                \code{"poisgeninvgauss"}. See details below.}
  \item{nchains}{Number of model chains, should be \eqn{\ge 2};
                 defaults to 2.}
  \item{burnin}{Number of samples to discard as burn-in; defaults to 1000.}
  \item{update}{Number of samples to retain; defaults to 5000.}
  \item{verbose}{Should JAGS process info be printed to the R console?
                 defaults to \code{FALSE}.}
}

\details{
  Distributions available:

\itemize{
  \item{\strong{Poisson:  }}{\code{family = "poisson"}}
  \item{\strong{Negative Binomial:  }}{\code{family = "negbin"}}
  \item{\strong{Poisson-LogNormal:  }}{\code{family = "poislognorm"}}
  \item{\strong{Poisson-Inverse Gaussian:  }}{\code{family = "poisinvgauss"}}
  \item{\strong{Poisson-Generalized Inverse Gaussian:  }}{\code{family = "poisgeninvgauss"}}
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
bea_count(x = counts, q = volume, family = "negbin", data = crypto)
}

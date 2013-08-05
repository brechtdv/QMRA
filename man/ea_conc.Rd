\name{ea_conc}
\alias{ea_conc}

\title{Exposure assessment from concentration data}

\description{
Parametric modelling of exposure concentration from (possibly left-censored)
concentration data.
}

\usage{
ea_conc(x, d, data,
        family = c("gamma", "lognormal", "weibull", "invgauss"))
}

\arguments{
  \item{x}{The observed concentrations.}
  \item{d}{A vector indicating if the observations were left-censored (1) or
           not (0).}
  \item{data}{An optional data frame, containing the variables in the model.
              If not found in \code{data}, the variables are taken from the
              environment from which \code{ea_count} is called.}
  \item{family}{A character string naming the family function;
                should be one of \code{"gamma"}, \code{"lognormal"},
                \code{"weibull"}, \code{"invgauss"}. See details below.}
}

\details{
  Distributions available:

\itemize{
  \item{\strong{Gamma:  }}{\code{family = "gamma"}}
  \item{\strong{Log-Normal:  }}{\code{family = "lognormal"}}
  \item{\strong{Weibull:  }}{\code{family = "weibull"}}
  \item{\strong{Inverse Gaussian:  }}{\code{family = "invgauss"}}
}}

\value{An object of class \code{"\linkS4class{ea}"}.}

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
  \code{\link{ea_count}}, for modelling exposure from count data\cr
  \code{\link{ea_presence}}, for modelling exposure from presence/absence data

}

\examples{
ea_conc(x = x, d = d, data = giardia, family = "weibull")
}

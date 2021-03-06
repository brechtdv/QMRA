\name{ea_conc}
\alias{ea_conc}

\title{Exposure assessment from concentration data}

\description{
Parametric modelling of exposure concentration from (possibly left-censored)
concentration data.
}

\usage{
ea_conc(x, d, data,
        model = c("gamma", "g",
                  "lognorm", "ln",
                  "weibull", "w",
                  "invgauss", "ig"),
        \dots)
}

\arguments{
  \item{x}{The observed concentrations.}
  \item{d}{A vector indicating if the observations were left-censored (1) or
           not (0).}
  \item{data}{An optional data frame, containing the variables in the model.
              If not found in \code{data}, the variables are taken from the
              environment from which \code{ea_count} is called.}
  \item{model}{A character string naming the concentration model. See details below.}
  \item{\dots}{Arguments to be passed to \code{\link{mle}}}
}

\details{
  Distributions available:

\itemize{
  \item{\strong{Gamma:  }}{\code{model = "gamma"} or \code{"g"}}
  \item{\strong{Log-Normal:  }}{\code{model = "lognorm"} or \code{"ln"}}
  \item{\strong{Weibull:  }}{\code{model = "weibull"} or \code{"w"}}
  \item{\strong{Inverse Gaussian:  }}{\code{model = "invgauss"} or \code{"ig"}}
}}

\value{An object of class \code{"\linkS4class{ea}"}.}

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
  \code{\link{ea_count}}, for modelling exposure from count data\cr
  \code{\link{ea_presence}}, for modelling exposure from presence/absence data
}

\examples{
ea_conc(x = x, d = d, data = giardia, model = "weibull")
}

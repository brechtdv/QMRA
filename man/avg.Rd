\name{avg}
\alias{avg}

\title{Model averaging}

\description{
Average model estimates according to Akaike's Information Criterion.
}

\usage{
avg(\dots)
}

\arguments{
  \item{\dots}{Possibly named objects of class \code{"\linkS4class{ea}"} or \code{"\linkS4class{bea}"}.}
}

\details{
Starting from a set of \eqn{K} plausible candidate models, the averaged model estimate \eqn{\hat{\theta}} is defined as \eqn{\sum_{i=1}^{K}{w_i \theta_i}}. The weights are defined as
\deqn{w_i = \frac{exp(-0.5 \Delta_i)}{\sum_{j=1}^{K}{exp(-0.5 \Delta_j)}}}
with \eqn{\Delta_i = AIC_i - AIC_{min}}

The variance of the averaged model estimate is given by
\deqn{var(\hat{\theta}) = [\sum_{i=1}^{K}{w_i \sqrt{var(\theta_i) + (\theta_i - \hat{\theta})^2}}]^2}
}

\value{An object of S3 class \code{"avg"}.}

\author{\email{Brecht.Devleesschauwer@UGent.be}}

\references{
\itemize{
  \item{
    Haas CN, Rose JB, Gerba CP (1999)
    \emph{Quantitative Microbial Risk Assessment.}
	John Wiley & Sons, Inc.
	}
  \item{
    Burnham KP, Anderson DR (2002)
    \emph{Model Selection and Multimodel Inference.}
  Springer-Verlag New York, Inc.
	}
}}

\examples{
## exposure assessment from concentration data
gam <- ea_conc(x = x, d = d, data = giardia, family = "gamma")
lno <- ea_conc(x = x, d = d, data = giardia, family = "lognormal")
wei <- ea_conc(x = x, d = d, data = giardia, family = "weibull")
inv <- ea_conc(x = x, d = d, data = giardia, family = "invgauss")

## model averaging
avg("Gamma" = gam,
    "Log-normal" = lno,
    "Weibull" = wei,
    "Inverse Gaussian" = inv)
}

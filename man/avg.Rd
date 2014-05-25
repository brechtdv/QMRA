\name{avg}
\alias{avg}
\alias{avg_ea}
\alias{avg_drm}

\title{Model averaging}

\description{
Average model estimates according to Akaike's Information Criterion.
}

\usage{
avg_ea(\dots)
avg_drm(\dots, dose)
}

\arguments{
  \item{\dots}{Possibly named objects of class \code{"\linkS4class{ea}"} or \code{"\linkS4class{drm}"}.}
  \item{dose}{Doses for which model-averaged responses need to be obtained.}
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
  \item{
  Namata H, Aerts M, Faes C, Teunis P (2008). Model averaging in microbial risk assessment using fractional polynomials. \emph{Risk analysis}, 28(4), 891-905.
  }
}}

\examples{
## Exposure assessment from concentration data
gam <- ea_conc(x = x, d = d, data = giardia, model = "gamma")
lno <- ea_conc(x = x, d = d, data = giardia, model = "lognormal")
wei <- ea_conc(x = x, d = d, data = giardia, model = "weibull")
inv <- ea_conc(x = x, d = d, data = giardia, model = "invgauss")

## Model averaging
avg_ea("Gamma" = gam,
       "Log-Normal" = lno,
       "Weibull" = wei,
       "Inverse Gaussian" = inv)


## Fit several dose-response models to the Campylobacter dataset
bp <- drm(infected, total, dose, campy, "betapoisson")
ll <- drm(infected, total, dose, campy, "loglogistic")
lp <- drm(infected, total, dose, campy, "logprobit")
ev <- drm(infected, total, dose, campy, "extremevalue")

## Model averaging
avg_drm("bp" = bp, "ll" = ll, "lp" = lp, "ev" = ev,
        dose = c(1, 10, 100))
}

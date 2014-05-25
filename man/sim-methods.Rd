\name{sim-methods}
\docType{methods}
\alias{sim}
\alias{sim-methods}
\alias{sim,ea-method}
\alias{sim,bea-method}
\alias{sim,drm-method}
\alias{sim,bdrm-method}

\title{Methods for Function \code{sim} in Package \pkg{QMRA}}

\description{
Simulate concentration from objects of class \code{"\linkS4class{ea}"}, \code{"\linkS4class{bea}"} and dose-reponse from objects of class \code{"\linkS4class{drm}"} and \code{"\linkS4class{bdrm}"}.
}

\usage{
\S4method{sim}{ea}(x, n, \dots)

\S4method{sim}{bea}(x, n, \dots)

\S4method{sim}{drm}(x, n, dose, \dots)

\S4method{sim}{bdrm}(x, n, dose, \dots)
}

\arguments{
  \item{x}{Object of class \code{"\linkS4class{ea}"}, \code{"\linkS4class{bea}"}, \code{"\linkS4class{drm}"} or \code{"\linkS4class{bdrm}"}.}
  \item{n}{Number of simulations.}
  \item{dose}{Doses for which the reponse needs to be simulated.}
  \item{\dots}{Other arguments to be passed to \code{\link{sim}}; currently ignored.}
}

\section{Methods}{
\describe{
  \item{\code{signature(x = "bdrm")}}{}
  \item{\code{signature(x = "bea")}}{}
  \item{\code{signature(x = "drm")}}{}
  \item{\code{signature(x = "ea")}}{}
}}

\examples{
## Exposure assessment from count data
ea_fit <- ea_count(x = counts, q = volume,
                   data = crypto, model = "negbin")
## Simulate 100 concentrations from fitted distribution
sim(ea_fit, n = 100)

## Fit Beta-Poisson model to Campylobacter dataset
drm_fit <- drm(x = infected, n = total, dose = dose,
               data = campy, model = "bp")
## Simulate 10 probabilities of response for doses 1:10
sim(drm_fit, n = 10, dose = 1:10)
}

\keyword{methods}

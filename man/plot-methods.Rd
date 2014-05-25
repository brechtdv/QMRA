\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ANY,ANY-method}
\alias{plot,bdrm,ANY-method}
\alias{plot,drm,ANY-method}

\title{Methods for Function \code{plot} in Package \pkg{QMRA}}

\description{
Plot dose-response curves from objects of class \code{"\linkS4class{drm}"} and \code{"\linkS4class{bdrm}"}.
}

\usage{
\S4method{plot}{drm,ANY}(x, y, se = TRUE, add = FALSE, n = NULL, min_log10dose = 0,
     xlab = "log10(dose)", ylab = "P(infection)",
     type = "l", lwd = 2, col = "red",
     se_pars = list(type = "l", lty = 2, lwd = 2, col = "blue"),
     sim_pars = list(type = "l", lty = 1, col = rgb(0, 0, 0, 0.1)),
     \dots)

\S4method{plot}{bdrm,ANY}(x, y, n = 100, add = FALSE, \dots)}

\arguments{
  \item{x}{Object of class \code{"\linkS4class{drm}"} or \code{"\linkS4class{bdrm}"}.}
  \item{y}{Ignored.}
  \item{se}{Logical. Should a 95\% confidence interval be plotted? Defaults to \code{TRUE}.}
  \item{add}{Logical. If \code{TRUE}, plots are added to current one, using \code{lines}.}
  \item{n}{Number of simulated dose-response curves to be plotted (generated using \code{\link{sim}}. Defaults to \code{NULL}.}
  \item{min_log10dose}{Lower bound of x axis; defaults to 0.}
  \item{xlab,ylab}{Titles for x and y axes.}
  \item{type,lwd,col}{Type, line width and color of dose-response curve.}
  \item{se_pars}{List of graphical parameters for the 95\% confidence interval curves; see \code{par}.}
  \item{sim_pars}{List of graphical parameters for the simulated dose-response curves; see \code{par}.}
  \item{\dots}{Other arguments to be passed to \code{plot}.}
}

\section{Methods}{
\describe{
  \item{\code{signature(x = "ANY", y = "ANY")}}{}
  \item{\code{signature(x = "bdrm", y = "ANY")}}{}
  \item{\code{signature(x = "drm", y = "ANY")}}{}
}}

\examples{
## Fit Beta-Poisson model to Campylobacter dataset
drm_fit <- drm(x = infected, n = total, dose = dose,
               data = campy, model = "bp")
## Plot fitted dose-response curve and 100 simulated curves
plot(drm_fit, se = TRUE, sim = 100)
           
## Fit Bayesian Beta-Poisson model to Campylobacter dataset
bdrm_fit <- bdrm(x = infected, n = total, dose = dose,
                 data = campy, model = "bp",
                 inits = list(alpha = 1, beta = 1))
## Plot 100 simulated dose-response curves
plot(bdrm_fit, sim = 100)
}

\keyword{methods}

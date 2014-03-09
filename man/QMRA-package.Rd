\name{QMRA-package}
\alias{QMRA-package}
\alias{QMRA}

\docType{package}

\title{Parametric Models for Quantitative Microbial Risk Assessment}

\description{
  The \pkg{QMRA} package provides Maximum-likelihood and Bayesian parametric
  methods for exposure assessment and dose-response assessment.
}

\details{
  \tabular{ll}{
    Package: \tab QMRA\cr
    Type: \tab Package\cr
    Version: \tab 0.0.7\cr
    Date: \tab 2014-03-09\cr
    Authors: \tab Brecht Devleesschauwer, Niko Speybroeck, Christel Faes\cr
    Maintainer: \tab \email{Brecht.Devleesschauwer@UGent.be}\cr
    BugReports: \tab \url{https://github.com/brechtdv/QMRA/issues}\cr
    Depends: \tab R (>= 3.0.0), rjags, coda, methods, lme4, stats4, SuppDists, HyperbolicDist\cr
    SystemRequirements: \tab jags (>= 3.2.0)\cr
    License: \tab GPL (>= 2)\cr
    LazyLoad: \tab yes\cr
}

Available functions in the \pkg{QMRA} package:
  \tabular{ll}{
    \code{\link{ea_count}} \tab Maximum-likelihood-based parametric modelling of exposure concentration from count data.\cr
	\code{\link{ea_conc}} \tab Maximum-likelihood-based parametric modelling of exposure concentration from (possibly left-censored) concentration data.\cr
    \code{\link{ea_presence}} \tab Maximum-likelihood-based parametric modelling of exposure concentration from presence/absence data.\cr
    \code{\link{bea_count}} \tab Bayesian parametric modelling of exposure concentration from count data.\cr
	\code{\link{bea_conc}} \tab Bayesian parametric modelling of exposure concentration from (possibly left-censored) concentration data.\cr
    \code{\link{bea_presence}} \tab Bayesian parametric modelling of exposure concentration from presence/absence data.\cr
    \code{\link{avg}} \tab Model averaging.\cr
    \code{\link{sim}} \tab Simulate from fitted concentration distributions. Method for objects of class \code{"\linkS4class{ea}"} and \code{"\linkS4class{bea}"}.\cr
    \code{\link{summarize}} \tab Summarize fitted concentration distributions. Method for objects of class \code{"\linkS4class{ea}"} and \code{"\linkS4class{bea}"}.\cr
  }
}

\note{
The Bayesian functions in the \pkg{prevalence} package call on \pkg{JAGS (Just Another Gibbs Sampler)}, through the \pkg{\link[rjags:rjags]{rjags}} package. Therefore, JAGS has to be installed on the user's system.\cr\cr
JAGS can be downloaded from \url{http://sourceforge.net/projects/mcmc-jags/}
}

\keyword{package}

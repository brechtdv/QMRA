\name{QMRA-package}
\alias{QMRA-package}
\alias{QMRA}

\docType{package}

\title{Parametric Models for Quantitative Microbial Risk Assessment}

\description{
  The \pkg{QMRA} package provides frequentist and Bayesian parametric models
  for exposure and dose-response assessment.
}

\details{
  \tabular{ll}{
    Package: \tab QMRA\cr
    Type: \tab Package\cr
    Version: \tab 0.0.1\cr
    Date: \tab 2013-08-01\cr
    Authors: \tab Brecht Devleesschauwer, Niko Speybroeck, Christel Faes\cr
    Maintainer: \tab \email{Brecht.Devleesschauwer@UGent.be}\cr
    BugReports: \tab \url{https://github.com/brechtdv/QMRA/issues}\cr
    Depends: \tab R (>= 3.0.0), rjags, coda, methods, lme4\cr
    SystemRequirements: \tab jags (>= 3.2.0)\cr
    License: \tab GPL (>= 2)\cr
    URL: \tab \url{http://users.ugent.be/~bdvleess/R/QMRA/}\cr
    LazyLoad: \tab yes\cr
}

Available functions in the \pkg{QMRA} package:
  \tabular{ll}{
    \code{\link{ea_count}} \tab Parametric modelling of exposure concentration from count data.\cr
	\code{\link{ea_conc}} \tab Parametric modelling of exposure concentration from (possibly left-censored) concentration data.\cr
    \code{\link{ea_presence}} \tab Parametric modelling of exposure concentration from presence/absence data.\cr
  }
}

\note{
The Bayesian functions in the \pkg{prevalence} package call on \pkg{JAGS (Just Another Gibbs Sampler)}, through the \pkg{\link[rjags:rjags]{rjags}} package. Therefore, JAGS has to be installed on the user's system.\cr\cr
JAGS can be downloaded from \url{http://sourceforge.net/projects/mcmc-jags/}
}

\keyword{package}

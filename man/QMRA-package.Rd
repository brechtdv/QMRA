\name{QMRA-package}
\alias{QMRA-package}
\alias{QMRA}

\docType{package}

\title{Parametric Models for Quantitative Microbial Risk Assessment}

\description{
The \pkg{QMRA} package provides a suite of functions for performing exposure assessment, dose-response assessment and risk characterisation in the context of microbial risk assessment. Different parametric models are available and can be fitted in both a maximum likelihood and Bayesian framework.
}

\details{
  \tabular{ll}{
    Package: \tab QMRA\cr
    Version: \tab 0.0.14\cr
    Date: \tab 2016-03-12\cr
    Authors: \tab Brecht Devleesschauwer, Christel Faes, Arie Havelaar, Niko Speybroeck\cr
    Maintainer: \tab \email{brechtdv@gmail.com}\cr
    URL: \tab \url{https://github.com/brechtdv/QMRA}\cr
    BugReports: \tab \url{https://github.com/brechtdv/QMRA/issues}\cr
    Depends: \tab R (>= 3.2.0)\cr
    Imports: \tab methods, rjags, lme4, stats4, MASS, SuppDists, HyperbolicDist\cr
    SystemRequirements: \tab JAGS (>= 4.1.0) (see http://mcmc-jags.sourceforge.net)\cr
    License: \tab GPL (>= 2)\cr
    LazyLoad: \tab yes\cr
}

Available functions in the \pkg{QMRA} package:\cr

\emph{Exposure Assessment:}
  \tabular{ll}{
    \code{\link{ea_count}} \tab ML parametric modelling of exposure concentration from count data. \cr
	\code{\link{ea_conc}} \tab ML parametric modelling of exposure concentration from (possibly left-censored) concentration data.\cr
    \code{\link{ea_presence}} \tab ML parametric modelling of exposure concentration from presence/absence data.\cr
    \code{\link{bea_count}} \tab Bayesian parametric modelling of exposure concentration from count data.\cr
	\code{\link{bea_conc}} \tab Bayesian parametric modelling of exposure concentration from (possibly left-censored) concentration data.\cr
    \code{\link{bea_presence}} \tab Bayesian parametric modelling of exposure concentration from presence/absence data.\cr
    \code{\link{avg_ea}} \tab Model averaging for objects of class \code{"\linkS4class{ea}"}.\cr
    \code{\link{sim}} \tab Simulate from fitted concentration distributions. Method for objects of class \code{"\linkS4class{ea}"} and \code{"\linkS4class{bea}"}.\cr
    \code{\link{summarize}} \tab Summarize fitted concentration distributions. Method for objects of class \code{"\linkS4class{ea}"} and \code{"\linkS4class{bea}"}.
    }

\emph{Dose-Response Modelling:}
  \tabular{ll}{
    \code{\link{drm}} \tab Maximum likelihood parametric modelling of dose-response data.\cr
    \code{\link{bdrm}} \tab Bayesian parametric modelling of dose-response data.\cr
    \code{\link{avg_drm}} \tab Model averaging for objects of class \code{"\linkS4class{drm}"}.\cr
    \code{\link{predict}} \tab Predict from fitted dose-response curve. Method for objects of class \code{"\linkS4class{drm}"}.\cr
    \code{\link{sim}} \tab Simulate from fitted dose-response curve. Method for objects of class \code{"\linkS4class{drm}"} and \code{"\linkS4class{bdrm}"}.\cr
    \code{\link{summarize}} \tab Summarize fitted dose-response curve. Method for objects of class \code{"\linkS4class{drm}"} and \code{"\linkS4class{bdrm}"}.
    }

\emph{Datasets:}
  \tabular{ll}{
    \code{\link{crypto}}	\tab \emph{Cryptosporidium} oocyst counts.\cr
    \code{\link{giardia}}	\tab \emph{Giardia} oocyst concentrations.\cr
    \code{\link{coliform}}	\tab Coliform presence/absence data.\cr
    \code{\link{campy}}		\tab \emph{Campylobacter} dose-response data.\cr
  }
}

\note{
The Bayesian functions in the \pkg{prevalence} package call on \pkg{JAGS (Just Another Gibbs Sampler)}, through the \pkg{\link[rjags:rjags]{rjags}} package. Therefore, JAGS has to be installed on the user's system.\cr\cr
JAGS can be downloaded from \url{http://sourceforge.net/projects/mcmc-jags/}
}

\keyword{package}

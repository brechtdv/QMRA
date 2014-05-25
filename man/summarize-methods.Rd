\name{summarize-methods}
\docType{methods}
\alias{summarize}
\alias{summarize-methods}
\alias{summarize,ea-method}
\alias{summarize,bea-method}
\alias{summarize,bdrm-method}

\title{Methods for Function \code{summarize} in Package \pkg{QMRA}}

\description{
Summarize objects of class \code{"\linkS4class{ea}"}, \code{"\linkS4class{bea}"}, and \code{"\linkS4class{bdrm}"}. Currently, this function returns the mean, standard error and 95\% confidence interval for the fitted parameters.}

\usage{
\S4method{summarize}{ea}(x, \dots)

\S4method{summarize}{bea}(x, \dots)

\S4method{summarize}{bdrm}(x, \dots)
}

\arguments{
  \item{x}{Object of class \code{"\linkS4class{ea}"}, \code{"\linkS4class{bea}"}, or \code{"\linkS4class{bdrm}"}.}
  \item{\dots}{Other arguments to be passed to \code{summarize}; currently ignored.}
}

\section{Methods}{
\describe{
  \item{\code{signature(x = "bdrm")}}{}
  \item{\code{signature(x = "bea")}}{}
  \item{\code{signature(x = "ea")}}{}
}}

\keyword{methods}

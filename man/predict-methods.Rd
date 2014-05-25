\name{predict-methods}
\docType{methods}
\alias{predict-methods}
\alias{predict,ANY-method}
\alias{predict,drm-method}
\title{Methods for Function \code{predict} in Package \pkg{QMRA}}
\description{
Predict dose-response from objects of class \code{"\linkS4class{drm}"}.
}

\usage{\S4method{predict}{drm}(object, dose, conf_level = 0.95, \dots)}

\arguments{
  \item{object}{Object of class \code{"\linkS4class{drm}"}.}
  \item{dose}{Vector of doses for which the response needs to be predicted.}
  \item{conf_level}{Confidence level used in the construction of confidence intervals; defaults to 0.95.}
  \item{\dots}{Other arguments to be passed to \code{predict}; currently ignored.}
}

\details{
The construction of confidence intervals for the predicted probability of response is based on the Delta method.
}

\section{Methods}{
\describe{
  \item{\code{signature(object = "drm")}}{}
}}

\seealso{
  \code{\linkS4class{drm}}
}

\examples{
## Fit Beta-Poisson dose-response model to Campylobacter dataset
fit <- drm(x = infected, n = total, dose = dose,
           data = campy, model = "bp")

## Predict probability of response at dose 1
predict(fit, dose = 1)
}

\keyword{methods}

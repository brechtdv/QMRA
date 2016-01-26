\name{drm-class}
\Rdversion{1.1}
\docType{class}
\alias{drm-class}
\alias{print,drm-method}
\alias{show,drm-method}

\title{Class \code{"drm"}}

\description{
The \code{"drm"} class represents output from dose-response models.
Objects of class \code{"drm"} are created by calls to \code{\link{drm}}.
}

\section{Slots}{
Objects of class \code{"drm"} contain the following slots: 
  \describe{
    \item{\code{call}:}{Object of class \code{"language"};
                        original call to exposure assessment function}
    \item{\code{data}:}{Object of class \code{"data.frame"};
                        original input data}
    \item{\code{family}:}{Object of class \code{"function"};
                          distribution family used to fit the model}
    \item{\code{gof}:}{Object of class \code{"gof"};
                       Chi-square goodness of fit}
    \item{\code{AIC}:}{Object of class \code{"numeric"};
                       resulting model \code{\link{AIC}}}
    \item{\code{mle_call}:}{Object of class \code{"language"};
                            the call to \code{\link{mle}}}
    \item{\code{coef}:}{Object of class \code{"numeric"};
                        estimated parameters}
    \item{\code{fullcoef}:}{Object of class \code{"numeric"};
                            fixed and estimated parameters}
    \item{\code{vcov}:}{Object of class \code{"matrix"};
                        approximate variance-covariance matrix}
    \item{\code{min}:}{Object of class \code{"numeric"};
                       minimum value of objective function}
    \item{\code{details}:}{Object of class \code{"list"},
                           as returned from \code{\link{optim}}}
    \item{\code{minuslogl}:}{Object of class \code{"function"};
                             the negative loglikelihood function}
    \item{\code{nobs}:}{Object of class \code{"integer"};
                        the number of observations (often \code{NA},
                        when not set in call explicitly).}
    \item{\code{method}:}{Object of class \code{"character"};
                          the optimization method used}
  }
}

\section{Extends}{
Class \code{"\linkS4class{mle}"}, directly.
}

\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "drm")}}
    \item{show}{\code{signature(object = "drm")}}
  }
}

\author{\email{brechtdv@gmail.com}}

\seealso{
\code{\linkS4class{mle}}
}

\keyword{classes}

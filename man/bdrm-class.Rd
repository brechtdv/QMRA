\name{bdrm-class}
\Rdversion{1.1}
\docType{class}
\alias{bdrm-class}
\alias{print,bdrm-method}
\alias{show,bdrm-method}

\title{Class \code{"bdrm"}}

\description{
The \code{"bdrm"} class represents output from Bayesian dose-response models. Objects of class \code{"bdrm"} are created by calls to \code{\link{bdrm}}.
}

\section{Slots}{
Objects of class \code{"bdrm"} contain the following slots: 
  \describe{
    \item{\code{call}:}{Object of class \code{"language"};
                        original call to exposure assessment function}
    \item{\code{data}:}{Object of class \code{"data.frame"};
                        original input data}
    \item{\code{family}:}{Object of class \code{"function"};
                          distribution family used to fit the model}
    \item{\code{par}:}{Object of class \code{"list"};
                       settings for Bayesian model fitting}
    \item{\code{model}:}{Object of class \code{"JAGS_model"};
                         fitted Bayesian model}
    \item{\code{mcmc}:}{Object of class \code{"list"};
                        resulting model estimates}
    \item{\code{diagnostics}:}{Object of class \code{"list"};
                               DIC and BGR diagnostics}
  }
}

\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "bdrm")}}
    \item{show}{\code{signature(object = "bdrm")}}
  }
}

\author{\email{Brecht.Devleesschauwer@UGent.be}}

\keyword{classes}

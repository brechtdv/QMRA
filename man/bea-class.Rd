\name{bea-class}
\Rdversion{1.1}
\docType{class}
\alias{bea-class}
\alias{print,bea-method}
\alias{show,bea-method}

\title{Class \code{"bea"}}

\description{
The \code{"bea"} class represents output from Bayesian exposure assessment
models. Objects of class \code{"bea"} are created by calls to
\code{\link{bea_count}}, \code{\link{bea_conc}} and \code{\link{bea_presence}}.
}

\section{Slots}{
Objects of class \code{"bea"} contain the following slots: 
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
    \item{print}{\code{signature(x = "bea")}: ... }
    \item{show}{\code{signature(object = "bea")}: ... }
  }
}

\author{\email{Brecht.Devleesschauwer@UGent.be}}

\seealso{
\code{\link{bea_count}}, \code{\link{bea_conc}}, \code{\link{bea_presence}}
}

\keyword{classes}

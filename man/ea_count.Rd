\name{ea_count}
\alias{ea_count}

\title{Exposure assessment from count data}

\description{
Parametric modelling of exposure concentration from count data.
}

\usage{
ea_count(x, q = 1, data,
         family = c("poisson", "negbin", "pln", "pig", "pgig"))
}

\arguments{
  \item{x}{The observed counts.}
  \item{q}{The quantities (e.g., volumes, masses) in which the counts were
           observed; defaults to 1.}
  \item{data}{An optional data frame, containing the variables in the model.
              If not found in \code{data}, the variables are taken from the
              environment from which \code{ea_count} is called.}
  \item{family}{A character string naming the family function;
                should be one of \code{"poisson"}, \code{"negbin"},
                \code{"pln"}, \code{"pig"}, \code{"pgig"}. See details below.}
}

\details{
  Distributions available:

\itemize{
  \item{\strong{Poisson:  }}{\code{family = "poisson"}}
  \item{\strong{Negative Binomial:  }}{\code{family = "negbin"}}
  \item{\strong{Poisson-LogNormal:  }}{\code{family = "pln"}}
  \item{\strong{Poisson-Inverse Gaussian:  }}{\code{family = "pig"}}
  \item{\strong{Poisson-Generalized Inverse Gaussian:  }}{\code{family = "pgig"}}
}}

\value{
An object of class \link[QMRA:ea-class]{ea}.
}

\author{
\email{Brecht.Devleesschauwer@UGent.be}
}

\references{
\itemize{
  \item{
    Haas CN, Rose JB, Gerba CP (1999)
    \emph{Quantitative Microbial Risk Assessment.}
	John Wiley & Sons, Inc.
	}
}}

\seealso{
  \code{\link{ea_presence}}, for modelling exposure from presence/absence data\cr
  \code{\link{ea_conc}}, for modelling exposure from concentration data
}

\examples{
ea_count(x = counts, q = volume, family = "negbin", data = crypto)
}

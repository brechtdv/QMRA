\name{ea_count}
\alias{ea_count}

\title{Exposure assessment from count data}

\description{
Parametric modelling of exposure concentration from count data.
}

\usage{
ea_count(x, q = 1, data,
         model = c("poisson", "p",
                   "negbin", "nb",
                   "poislognorm", "pln",
                   "poisinvgauss", "pig",
                   "poisgeninvgauss", "pgig"),
         \dots)
}

\arguments{
  \item{x}{The observed counts.}
  \item{q}{The quantities (e.g., volumes, masses) in which the counts were
           observed; defaults to 1.}
  \item{data}{An optional data frame, containing the variables in the model.
              If not found in \code{\link{data}}, the variables are taken from the
              environment from which \code{\link{ea_count}} is called.}
  \item{model}{A character string naming the count model. See details below.}
  \item{\dots}{Arguments to be passed to \code{\link{mle}}}
}

\details{
  Distributions available:

\itemize{
  \item{\strong{Poisson:  }}{\code{model = "poisson"} or \code{"p"}}
  \item{\strong{Negative Binomial:  }}{\code{model = "negbin"} or \code{"nb"}}
  \item{\strong{Poisson-LogNormal:  }}{\code{model = "poislognorm"} or \code{"pln"}}
  \item{\strong{Poisson-Inverse Gaussian:  }}{\code{model = "poisinvgauss"} or \code{"pig"}}
  \item{\strong{Poisson-Generalized Inverse Gaussian:  }}{\code{model = "poisgeninvgauss"} or \code{"pgig"}}
}}

\value{An object of class \code{"\linkS4class{ea}"}.}

\author{\email{brechtdv@gmail.com}}

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
ea_count(x = counts, q = volume, model = "negbin", data = crypto)
}

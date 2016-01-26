\name{ea_presence}
\alias{ea_presence}

\title{Exposure assessment from presence/absence data}

\description{
Parametric modelling of exposure concentration from presence/absence data.
}

\usage{
ea_presence(x, q = 1, replicates = rep(1, length(x)), data,
            model = c("poisson", "p"), \dots)
}

\arguments{
  \item{x}{A vector indicating the number of positive samples per examined
           quantity.}
  \item{q}{The quantities (e.g., volumes, masses) in which presence/absence
           was observed; defaults to 1.}
  \item{replicates}{The number of replicates of each quantity \code{q};
                    defaults to \code{rep(1, length(x))}.}
  \item{data}{An optional data frame, containing the variables in the model.
              If not found in \code{data}, the variables are taken from the
              environment from which \code{ea_presence} is called.}
  \item{model}{A character string naming the model to be fitted. See details below.}
  \item{\dots}{Arguments to be passed to \code{\link{mle}}}
}

\details{
  Available distributions:

\itemize{
  \item{\strong{Poisson:  }}{\code{model = "poisson"} or \code{"p"}}
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
  \code{\link{ea_count}}, for modelling exposure from count data\cr
  \code{\link{ea_conc}}, for modelling exposure from concentration data
}

\examples{
## exposure assessment from presence/absence data
ea_presence(x = positive, q = volume, rep = replicates, data = coliform)
}

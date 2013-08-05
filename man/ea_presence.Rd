\name{ea_presence}
\alias{ea_presence}

\title{Exposure assessment from presence/absence data}

\description{
Parametric modelling of exposure concentration from presence/absence data.
}

\usage{
ea_presence(x, q = 1, replicates = rep(1, length(x)), data,
            family = c("poisson"))
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
  \item{family}{A character string naming the family function;
                currently only \code{"poisson"} is supported.
                See details below.}
}

\details{
  Available distributions:

\itemize{
  \item{\strong{Poisson:  }}{\code{family = "poisson"}}
}}

\value{An object of class \code{"\linkS4class{ea}"}.}

\author{\email{Brecht.Devleesschauwer@UGent.be}}

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
ea_presence(x = positive, q = volume, rep = replicates, data = dilution)
}
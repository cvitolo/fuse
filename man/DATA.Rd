\name{DATA}
\alias{DATA}
\docType{data}
\title{
Example forcing time series
}
\description{
This dataset consists of a data.frame with 3 columns: P (precipitation time series), E(potential evapo-transpiration time series), Q(discharge time series).
}
\usage{data(DATA)}
\format{
  A data frame with 9480 observations on the following 3 variables.
  \describe{
    \item{\code{P}}{precipitation time series}
    \item{\code{E}}{potential evapo-transpiration time series}
    \item{\code{Q}}{discharge time series}
  }
}
\details{
This is an example dataset.
}
\examples{
data(DATA)
str(DATA)
}
\keyword{datasets}

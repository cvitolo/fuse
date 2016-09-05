\name{parameters}
\alias{parameters}
\docType{data}
\title{
Example FUSE parameter set
}
\description{
This dataset consists of a data.frame with 24 columns (1 each sampled parameter value)
}
\usage{data(parameters)}
\format{
  A data frame with 1 observations on the following 24 variables.
  \describe{
 \item{\code{rferr_add}}{additive rainfall error (mm)}
 \item{\code{rferr_mlt}}{multiplicative rainfall error (-)}
 \item{\code{maxwatr_1}}{depth of the upper soil layer (mm)}
 \item{\code{maxwatr_2}}{depth of the lower soil layer (mm)}
 \item{\code{fracten}}{fraction total storage in tension storage (-)}
 \item{\code{frchzne}}{fraction tension storage in recharge zone (-)}
 \item{\code{fprimqb}}{fraction storage in 1st baseflow reservoir (-)}
 \item{\code{rtfrac1}}{fraction of roots in the upper layer (-)}
 \item{\code{percrte}}{percolation rate (mm day-1)}
 \item{\code{percexp}}{percolation exponent (-)}
 \item{\code{sacpmlt}}{SACRAMENTO model percolation multiplier for dry soil layer (-)}
 \item{\code{sacpexp}}{SACRAMENTO model percolation exponential coefficient for dry soil layer (-)}
 \item{\code{percfrac}}{fraction of percolation to tension storage (-)}
 \item{\code{iflwrte}}{interflow rate (mm day-1)}
 \item{\code{baserte}}{baseflow rate (mm day-1)}
 \item{\code{qb_powr}}{baseflow exponent (-)}
 \item{\code{qb_prms}}{baseflow depletion rate (day-1)}
 \item{\code{qbrate_2a}}{baseflow depletion rate 1st reservoir (day-1)}
 \item{\code{qbrate_2b}}{baseflow depletion rate 2nd reservoir (day-1)}
 \item{\code{sareamax}}{maximum saturated area (-)}
 \item{\code{axv_bexp}}{ARNO/VIC model "b" exponent (-)}
 \item{\code{loglamb}}{mean value of the topographic index (m)}
 \item{\code{tishape}}{shape parameter for the topographic index Gamma distribution (-)}
 \item{\code{timedelay}}{time delay in runoff (days)}
  }
}
\details{
This is an example dataset.
}
\examples{
data(parameters)
str(parameters)
}
\keyword{datasets}

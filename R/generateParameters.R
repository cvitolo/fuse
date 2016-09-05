#' This function generates parameter sets for FUSE models (see ranges suggested by Clark et al., 2011).
#'
#' @param NumberOfRuns number of samples to generate, can be any integer
#' @param SamplingType sampling procedure to use, can be "URS" or "LHS"
#' @param rferr_add range of the additive rainfall error (mm), default is c(0,0)
#' @param rferr_mlt range of the multiplicative rainfall error (-), default is c(1,1)
#' @param frchzne range of the fraction tension storage in recharge zone (-), default is c(0.05,0.95)#'
#' @param fracten range of the fraction total storage in tension storage (-), default is c(0.05,0.95)
#' @param maxwatr_1 range of the depth of the upper soil layer (mm), default is c(25,500)
#' @param percfrac range of the fraction of percltn to tension storage (-), default is c(0.05,0.95)
#' @param fprimqb range of the fraction storage in 1st baseflow reservoir (-), default is c(0.05,0.95)
#' @param qbrate_2a range of the baseflow depletion rate 1st reservoir (day-1), default is c(0.001,0.25)
#' @param qbrate_2b range of the baseflow depletion rate 2nd reservoir (day-1), default is c(0.001,0.25)
#' @param qb_prms range of the baseflow depletion rate (day-1), default is c(0.001,0.25)
#' @param maxwatr_2 range of the depth of the lower soil layer (mm), default is c(50,5000)
#' @param baserte range of the baseflow rate (mm day-1), default is c(0.001,1000)
#' @param rtfrac1 range of the fraction of roots in the upper layer (-), default is c(0.05,0.95)
#' @param percrte range of the percolation rate (mm day-1), default is c(0.01,1000)
#' @param percexp range of the percolation exponent (-), default is c(1,20)
#' @param sacpmlt range of the SAC model percltn mult for dry soil layer (-), default is c(1,250)
#' @param sacpexp range of the SAC model percltn exp for dry soil layer (-), default is c(1,5)
#'
#' @param iflwrte range of the interflow rate (mm day-1), default is c(0.01,1000)
#' @param axv_bexp range of the ARNO/VIC b exponent (-), default is c(0.001,3)
#' @param sareamax range of the maximum saturated area (-), default is c(0.05,0.95)
#' @param loglamb range of the mean value of the topographic index (m), default is c(5,10)
#' @param tishape range of the shape param for the topo index Gamma dist (-), default is c(2,5)
#' @param qb_powr range of the baseflow exponent (-), default is c(1,10)
#' @param timedelay range of the time delay in runoff (days), default is c(0.01,5)
#' @param params2remove vector with names of parameters to remove from the latin hypercube.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # For reproducible results, use set.seed() before running this function.
#' set.seed(123)
#' parameters <- generateParameters(1)
#' }
#'

generateParameters <- function(NumberOfRuns, SamplingType = "LHS",
                               rferr_add = NULL,
                               rferr_mlt = NULL,
                               frchzne = NULL,
                               fracten = NULL,
                               maxwatr_1 = NULL,
                               percfrac = NULL,
                               fprimqb = NULL,
                               qbrate_2a = NULL,
                               qbrate_2b = NULL,
                               qb_prms = NULL,
                               maxwatr_2 = NULL,
                               baserte = NULL,
                               rtfrac1 = NULL,
                               percrte = NULL,
                               percexp = NULL,
                               sacpmlt = NULL,
                               sacpexp = NULL,
                               iflwrte = NULL,
                               axv_bexp = NULL,
                               sareamax = NULL,
                               loglamb = NULL,
                               tishape = NULL,
                               qb_powr = NULL,
                               timedelay = NULL,
                               params2remove = NULL) {

  # require(tgp)
  # rferr_add <- rferr_mlt <- maxwatr_1 <- maxwatr_2 <- fracten <- frchzne <-   fprimqb <- rtfrac1 <- percrte <- percexp <- sacpmlt <- sacpexp <- percfrac <-   iflwrte <- baserte <- qb_powr <- qb_prms <- qbrate_2a <- qbrate_2b <- sareamax <- axv_bexp <- loglamb <- tishape <- timedelay <- NULL

  if ( is.null(rferr_add) ) rferr_add = c(0,0) # additive rainfall error (mm)
  if ( is.null(rferr_mlt) ) rferr_mlt = c(1,1) # multiplicative rainfall error (-)
  if ( is.null(maxwatr_1) ) maxwatr_1 = c(25,500) # depth of the upper soil layer (mm)
  if ( is.null(maxwatr_2) ) maxwatr_2 = c(50,5000) # depth of the lower soil layer (mm)
  if ( is.null(fracten) ) fracten   = c(0.05,0.95) # fraction total storage in tension storage (-)
  if ( is.null(frchzne) ) frchzne   = c(0.05,0.95) # fraction tension storage in recharge zone (-)
  if ( is.null(fprimqb) ) fprimqb   = c(0.05,0.95) # fraction storage in 1st baseflow reservoir (-)
  if ( is.null(rtfrac1) ) rtfrac1   = c(0.05,0.95) # fraction of roots in the upper layer (-)
  if ( is.null(percrte) ) percrte   = c(0.01,1000) # percolation rate (mm day-1)
  if ( is.null(percexp) ) percexp   = c(1,20) # percolation exponent (-)
  if ( is.null(sacpmlt) ) sacpmlt   = c(1,250) # SAC model percltn mult for dry soil layer (-)
  if ( is.null(sacpexp) ) sacpexp   = c(1,5) # SAC model percltn exp for dry soil layer (-)
  if ( is.null(percfrac) ) percfrac  = c(0.05,0.95) # fraction of percltn to tension storage (-)
  if ( is.null(iflwrte) ) iflwrte   = c(0.01,1000) # interflow rate (mm day-1)
  if ( is.null(baserte) ) baserte   = c(0.001,1000) # baseflow rate (mm day-1)
  if ( is.null(qb_powr) ) qb_powr   = c(1,10) # baseflow exponent (-)
  if ( is.null(qb_prms) ) qb_prms   = c(0.001,0.25) # baseflow depletion rate (day-1)
  if ( is.null(qbrate_2a) ) qbrate_2a = c(0.001,0.25) # baseflow depletion rate 1st reservoir (day-1)
  if ( is.null(qbrate_2b) ) qbrate_2b = c(0.001,0.25) # baseflow depletion rate 2nd reservoir (day-1)
  if ( is.null(sareamax) ) sareamax  = c(0.05,0.95) # maximum saturated area (-)
  if ( is.null(axv_bexp) ) axv_bexp  = c(0.001,3) # ARNO/VIC b exponent (-)
  if ( is.null(loglamb) ) loglamb   = c(5,10) # mean value of the topographic index (m)
  if ( is.null(tishape) ) tishape   = c(2,5) # shape param for the topo index Gamma dist (-)
  if ( is.null(timedelay) ) timedelay = c(0.01,5) # time delay in runoff (days)

  # Define sample domain (min and max for each parameter)
  DefaultRanges <- data.frame(rbind(rferr_add,rferr_mlt,maxwatr_1,maxwatr_2,
                                    fracten,frchzne,fprimqb,rtfrac1,percrte,
                                    percexp,sacpmlt,sacpexp,percfrac,iflwrte,
                                    baserte,qb_powr,qb_prms,qbrate_2a,qbrate_2b,
                                    sareamax,axv_bexp,loglamb,tishape,timedelay)
                              )
  names(DefaultRanges) <- c("Min","Max")

  if (SamplingType == "URS") {
    # This option generates "N" parameter sets sampling the ranges using
    # a Uniform Random distribution.

    rferr_add <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["rferr_add","Min"],
                       max=DefaultRanges["rferr_add","Max"])
    rferr_mlt <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["rferr_mlt","Min"],
                       max=DefaultRanges["rferr_mlt","Max"])
    maxwatr_1 <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["maxwatr_1","Min"],
                       max=DefaultRanges["maxwatr_1","Max"])
    maxwatr_2 <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["maxwatr_2","Min"],
                       max=DefaultRanges["maxwatr_2","Max"])
    fracten   <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["fracten","Min"],
                       max=DefaultRanges["fracten","Max"])
    frchzne   <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["frchzne","Min"],
                       max=DefaultRanges["frchzne","Max"])
    fprimqb   <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["fprimqb","Min"],
                       max=DefaultRanges["fprimqb","Max"])
    rtfrac1   <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["rtfrac1","Min"],
                       max=DefaultRanges["rtfrac1","Max"])
    percrte   <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["percrte","Min"],
                       max=DefaultRanges["percrte","Max"])
    percexp   <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["percexp","Min"],
                       max=DefaultRanges["percexp","Max"])
    sacpmlt   <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["sacpmlt","Min"],
                       max=DefaultRanges["sacpmlt","Max"])
    sacpexp   <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["sacpexp","Min"],
                       max=DefaultRanges["sacpexp","Max"])
    percfrac  <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["percfrac","Min"],
                       max=DefaultRanges["percfrac","Max"])
    iflwrte   <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["iflwrte","Min"],
                       max=DefaultRanges["iflwrte","Max"])
    baserte   <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["baserte","Min"],
                       max=DefaultRanges["baserte","Max"])
    qb_powr   <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["qb_powr","Min"],
                       max=DefaultRanges["qb_powr","Max"])
    qb_prms   <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["qb_prms","Min"],
                       max=DefaultRanges["qb_prms","Max"])
    qbrate_2a <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["qbrate_2a","Min"],
                       max=DefaultRanges["qbrate_2a","Max"])
    qbrate_2b <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["qbrate_2b","Min"],
                       max=DefaultRanges["qbrate_2b","Max"])
    sareamax  <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["sareamax","Min"],
                       max=DefaultRanges["sareamax","Max"])
    axv_bexp  <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["axv_bexp","Min"],
                       max=DefaultRanges["axv_bexp","Max"])
    loglamb   <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["loglamb","Min"],
                       max=DefaultRanges["loglamb","Max"])
    tishape   <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["tishape","Min"],
                       max=DefaultRanges["tishape","Max"])
    timedelay <- stats::runif(NumberOfRuns,
                       min=DefaultRanges["timedelay","Min"],
                       max=DefaultRanges["timedelay","Max"])

    parameters <- data.frame("rferr_add"=rferr_add,
                             "rferr_mlt"=rferr_mlt,
                             "maxwatr_1"=maxwatr_1,
                             "maxwatr_2"=maxwatr_2,
                             "fracten"=fracten,
                             "frchzne"=frchzne,
                             "fprimqb"=fprimqb,
                             "rtfrac1"=rtfrac1,
                             "percrte"=percrte,
                             "percexp"=percexp,
                             "sacpmlt"=sacpmlt,
                             "sacpexp"=sacpexp,
                             "percfrac"=percfrac,
                             "iflwrte"=iflwrte,
                             "baserte"=baserte,
                             "qb_powr"=qb_powr,
                             "qb_prms"=qb_prms,
                             "qbrate_2a"=qbrate_2a,
                             "qbrate_2b"=qbrate_2b,
                             "sareamax"=sareamax,
                             "axv_bexp"=axv_bexp,
                             "loglamb"=loglamb,
                             "tishape"=tishape,
                             "timedelay"=timedelay)

    parameters[,params2remove] <- -999

  }

  if (SamplingType == "LHS") {
    # This option generates "N" parameter sets sampling the ranges using
    # a Latin Hypercube.

    rferr_add <- rferr_mlt <- maxwatr_1 <- maxwatr_2 <- fracten <- frchzne <- NA
    fprimqb <- rtfrac1 <- percrte <- percexp <- sacpmlt <- sacpexp <- NA
    percfrac <- iflwrte <- baserte <- qb_powr <- qb_prms <- qbrate_2a <- NA
    qbrate_2b <- sareamax <- axv_bexp <- loglamb <- tishape <- timedelay <- NA

    if (!is.null(params2remove)) {
      rows2remove <- which(row.names(DefaultRanges) %in% params2remove)
      newRanges <- DefaultRanges[-rows2remove,]
    }else{
      newRanges <- DefaultRanges
    }

    temp <- tgp::lhs( NumberOfRuns, as.matrix(newRanges) )
    for (cols in 1:length(params2remove)){
      temp <- cbind(temp,rep(NA,dim(temp)[1]))
    }
    temp <- data.frame(temp)
    names(temp) <- c(row.names(newRanges),params2remove)

    parameters <- temp[,c("rferr_add","rferr_mlt","maxwatr_1","maxwatr_2",
                       "fracten","frchzne","fprimqb","rtfrac1","percrte",
                       "percexp","sacpmlt","sacpexp","percfrac","iflwrte",
                       "baserte","qb_powr","qb_prms","qbrate_2a",
                       "qbrate_2b","sareamax","axv_bexp","loglamb",
                       "tishape","timedelay")]

    parameters[is.na(parameters)] <- -999

  }

return(parameters)

}

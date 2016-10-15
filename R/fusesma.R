#' FUSE Soil Moisture Accounting model
#'
#' @description Implementation of the Soil Moisture Accounting module of the framework for hydrological modelling FUSE. This loss module is derived from the Fortran version of FUSE by Martyn Clark (2011).
#'
#' @param DATA This is a data.frame containing the observed time series (zoo objects). It is structured into three columns containing: precipitation (P), potential evapo-transpiration (E) and streamflow discharge (Q).
#' @param mid This is the model identification number (see first column of \code{modlist}).
#' @param deltim This is the input time step. By default deltim = 1 for daily steps. Other options: 1/24 (hourly time step), 1/(24*4) (15 min time step).
#' @param rferr_add Additive rainfall error, default is 0 (mm day-1).
#' @param rferr_mlt Multiplicative rainfall error, default is 1 (-).
#' @param frchzne Fraction of tension storage in recharge zone (-).
#' @param fracten Fraction total storage as tension storage (-).
#' @param maxwatr_1 Maximum total storage in upper soil layer (mm).
#' @param percfrac Fraction of percolation to tension storage in the lower layer (-).
#' @param fprimqb Fraction of storage in the first baseflow reservoir (-).
#' @param qbrate_2a Baseflow depletion rate in the first reservoir (day-1).
#' @param qbrate_2b Baseflow depletion rate in the second reservoir (day-1).
#' @param qb_prms Baseflow depletion rate (day-1).
#' @param maxwatr_2 Maximum total storage in lower soil layer (mm).
#' @param baserte Baseflow rate (mm day-1).
#' @param rtfrac1 Fraction of roots in the upper layer (-).
#' @param percrte Percolation rate (mm day-1).
#' @param percexp Percolation exponent (-).
#' @param sacpmlt Sacramento model percolation multiplier for dry soil layer (-).
#' @param sacpexp Sacramento model percolation exponent for dry soil layer (-).
#' @param iflwrte Interflow rate (mm day-1).
#' @param axv_bexp ARNO/VIC 'b' exponent (-).
#' @param sareamax Maximum saturated area (-).
#' @param loglamb Mean value of the log-transformed topographic index (m).
#' @param tishape Shape parameter for the topo index gamma distribution (-).
#' @param qb_powr Baseflow exponent (-).
#' @param fracstate0 Initial saturation of water storages, default is 0.25.
#' @param absError Absolute solver error, default is 10^(-4).
#' @param relError Relative solver error, default is 10^(-4).
#' @param StatesFluxes By default StatesFluxes = FALSE which means the function only output is the instantaneous runoff. If StatesFluxes=TRUE, the output also contains the list of fluxes and state variables.
#'
#' @details fusesma.sim() is a function to generate an ensemble of SOIL MOISTURE ACCOUNTING models. It is compatible with the HYDROMAD framework (see hydromad package: http://hydromad.catchment.org/). fusesma.sim() can simulate several model structures. The default list is \code{modlist} contained in the data folder of this package. The parameter set varies depending on the selected model structure. Ranges of parameter values are in \code{fusesma.ranges}. For more information on suggested parameter ranges see Clark et al. 2011. Also see \code{GenerateFUSEParameters}. Flow can then be routed using the function \code{fuserouting.sim} (which is based on the Gamma function) or any other routing function.
#'
#'
#' @return The function returns an array of simulated "instantaneous" discharges. If necessary, \code{fuserouting.sim} can be run to obtain routed discharges using a two parameter Gamma distribution.
#'
#' @references{
#' Clark M. P., SlaterA. G., Rupp D. E., Woods R. A., Vrugt J. A., Gupta H. V., Wagener T. and Hay L. E. (2008), Framework for Understanding Structural Errors (FUSE): A modular framework to diagnose differences between hydrological models, Water Resour. Res. 44 p. 91-94.
#' Clark M. P., McMillan H. K., Collins D. B. G., Kavetski D. and Woods R. A. (2011), Hydrological field data from a modeller's perspective: Part 2: process-based evaluation of model hypotheses. Hydrological Processes, 25: 523-543. doi: 10.1002/hyp.7902
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(DATA)
#' set.seed(123)
#' parameters <- generateParameters(1)
#' U <- fusesma.sim(DATA, 60, 1, parameters$rferr_add, parameters$rferr_mlt,
#'                  parameters$frchzne, parameters$fracten,
#'                  parameters$maxwatr_1, parameters$percfrac,
#'                  parameters$fprimqb, parameters$qbrate_2a,
#'                  parameters$qbrate_2b, parameters$qb_prms,
#'                  parameters$maxwatr_2, parameters$baserte,
#'                  parameters$rtfrac1, parameters$percrte, parameters$percexp,
#'                  parameters$sacpmlt, parameters$sacpexp, parameters$iflwrte,
#'                  parameters$axv_bexp, parameters$sareamax,
#'                  parameters$loglamb, parameters$tishape, parameters$qb_powr)
#' }
#'

fusesma.sim <- function(DATA, mid, deltim,
                        rferr_add = 0,
                        rferr_mlt = 1,
                        frchzne, fracten, maxwatr_1, percfrac, fprimqb,
                        qbrate_2a, qbrate_2b, qb_prms, maxwatr_2, baserte,
                        rtfrac1, percrte, percexp, sacpmlt, sacpexp, iflwrte,
                        axv_bexp, sareamax, loglamb, tishape, qb_powr,
                        fracstate0 = 0.25,
                        absError = 10^(-4),   # absolute solver error (default)
                        relError = 10^(-4),   # relative solver error (default)
                        StatesFluxes = FALSE){

    ## Keep attributes, but work with raw matrix
    inAttr <- attributes(DATA[, 1])
    DATA <- zoo::coredata(DATA)

    stopifnot(c("P","E") %in% colnames(DATA))
    P <- DATA[,"P"]
    E <- DATA[,"E"]

    ## skip over missing values
    #bad <- is.na(P) | is.na(E)
    #P[bad] <- 0
    #E[bad] <- 0

    # list of availabe models is loaded by LazyData: true in DESCRIPTION
    modlist <- modlist                           # to remove NOTE in R CMD check

    # Make sure mid is an integer (needed for compatibility with hydromad)
    mid <- round(mid,0)

    # Read model structure [LIST]
    smodl <- list("rferr"=modlist[mid,2],
                  "arch1"=modlist[mid,3],
                  "arch2"=modlist[mid,4],
                  "qsurf"=modlist[mid,5],
                  "qperc"=modlist[mid,6],
                  "esoil"=modlist[mid,7],
                  "qintf"=modlist[mid,8],
                  "q_tdh"=modlist[mid,9])

    res <- .C("fusesmaSimR",
              as.double(P),
              as.integer(length(P)),
              as.double(E),
              as.integer(length(E)),
              as.integer(smodl),       # A (1d) array of ints (the model params)
              as.integer(length(smodl)),
              as.double(deltim),
              as.double(frchzne),
              as.double(fracten),
              as.double(maxwatr_1),
              as.double(percfrac),
              as.double(fprimqb),
              as.double(qbrate_2a),
              as.double(qbrate_2b),
              as.double(qb_prms),
              as.double(maxwatr_2),
              as.double(baserte),
              as.double(rtfrac1),
              as.double(percrte),
              as.double(percexp),
              as.double(sacpmlt),
              as.double(sacpexp),
              as.double(iflwrte),
              as.double(axv_bexp),
              as.double(sareamax),
              as.double(loglamb),
              as.double(tishape),
              as.double(qb_powr),
              stateResults = double( 10 * length(P) ),
              fluxResults = double( 20 * length(P) ),
              as.double(absError),
              as.double(relError),
              as.logical(FALSE),                        # correctNegVols = FALSE
              as.double(fracstate0),
              as.double(rferr_add),
              as.double(rferr_mlt),
              status = integer(1) , PACKAGE="fuse")

    # TODO: add module for correction of negative values
    # (now correctNegVols = FALSE)

    if (res$status != 0){
      stop("fusesmaSimR Failed!")
    }

    s <- res$stateResults
    dim(s) <- c(length(P),10)
    s <- data.frame( matrix(s,nrow=length(P),ncol=10,byrow=T) )

    f <- res$fluxResults
    dim(f) <- c(length(P),20)
    f <- data.frame( matrix(f,nrow=length(P),ncol=20,byrow=T) )

    results <- cbind(s,f)
    names(results) <- c("tens_1a", "tens_1b", "tens_1", "free_1", "watr_1",
                        "tens_2", "free_2a", "free_2b", "watr_2", "free_2",
                        "eff_ppt", "satarea", "qsurf", "evap_1a", "evap_1b",
                        "evap_1", "evap_2", "rchr2excs", "tens2free_1",
                        "tens2free_2", "qintf_1", "qperc_12", "qbase_2",
                        "qbase_2a", "qbase_2b", "oflow_1", "oflow_2",
                        "oflow_2a", "oflow_2b", "U")

    if (StatesFluxes == FALSE) {
      results <- results$U
      attributes(results) <- inAttr
    }

    return(results)
}

#' Function to define the parameter ranges for FUSE Soil Moisture Accounting module
#'
#' @description Parameter ranges for FUSE Soil Moisture Accounting module, as suggested in Clark et al. (2011).
#'
#' @return The function returns the list of parameter ranges for the FUSE Soil Moisture Accounting model.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fusesma.ranges()
#' }
#'

fusesma.ranges <- function() {
    list("rferr_add" = 0,                         # additive rainfall error (mm)
         "rferr_mlt" = 1,                    # multiplicative rainfall error (-)
         "maxwatr_1" = c(25, 500),          # depth of the upper soil layer (mm)
         "maxwatr_2" = c(50, 5000),         # depth of the lower soil layer (mm)
         "fracten"   = c(0.05, 0.95), # fraction total storage in tension storage (-)
         "frchzne"   = c(0.05, 0.95), # fraction tension storage in recharge zone (-)
         "fprimqb"   = c(0.05, 0.95), # fraction storage in 1st baseflow reservoir (-)
         "rtfrac1"   = c(0.05, 0.95), # fraction of roots in the upper layer (-)
         "percrte"   = c(0.01, 1000),              # percolation rate (mm day-1)
         "percexp"   = c(1, 20),                      # percolation exponent (-)
         "sacpmlt"   = c(1, 250), # SAC model percltn mult for dry soil layer (-)
         "sacpexp"   = c(1, 5),   # SAC model percltn exp for dry soil layer (-)
         "percfrac"  = c(0.05, 0.95), # fraction of percltn to tension storage (-)
         "iflwrte"   = c(0.01, 1000),                # interflow rate (mm day-1)
         "baserte"   = c(0.001, 1000),                # baseflow rate (mm day-1)
         "qb_powr"   = c(1, 10),                         # baseflow exponent (-)
         "qb_prms"   = c(0.001, 0.25),         # baseflow depletion rate (day-1)
         "qbrate_2a" = c(0.001, 0.25), # baseflow depletion rate 1st reservoir (day-1)
         "qbrate_2b" = c(0.001, 0.25), # baseflow depletion rate 2nd reservoir (day-1)
         "sareamax"  = c(0.05, 0.95),               # maximum saturated area (-)
         "axv_bexp"  = c(0.001, 3),                  # ARNO/VIC "b" exponent (-)
         "loglamb"   = c(5, 10),       # mean value of the topographic index (m)
         "tishape"   = c(2, 5))  # shape param for the topo index Gamma dist (-)
}


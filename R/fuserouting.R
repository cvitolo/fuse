#' FUSE Routing model
#'
#' @description Implementation of the Routing module of the framework for hydrological modelling FUSE, derived from the Fortran version of FUSE by Martyn Clark (2011).
#'
#' @param U This is the effective rainfall/instantaneous runoff, which is the sum of surface runoff, overflow, interflow, and baseflow. This vector can either be of class zoo or numeric.
#' @param mid This is the model identification number (see first column of \code{modlist}). This is a numeric value in the range [1, 1248].
#' @param deltim This is the input time step. By default deltim = 1 for daily steps. All options: deltim = 1 (daily time step), 1/24 (hourly time step), 1/24/4 (15 min time step). This is a numeric value.
#' @param timedelay Time delay in runoff (days). This is a numeric value.
#'
#' @details fuserouting.sim() is a routing module based on a two parameter Gamma distribution. It takes in input the instantaneous discharge and returns the routed discharge. It is compatible with the HYDROMAD framework (see hydromad package: http://hydromad.catchment.org/). For more information on suggested parameter ranges see Clark et al. 2011. Also see \code{GenerateFUSEParameters}.
#'
#' @return The function returns an array of simulated "routed" discharges. It can be used after calculating instantaneous discharges with \code{fusesma.sim}.
#'
#' @references{
#' Clark M. P., SlaterA. G., Rupp D. E., Woods R. A., Vrugt J. A., Gupta H. V., Wagener T. and Hay L. E. (2008), Framework for Understanding Structural Errors (FUSE): A modular framework to diagnose differences between hydrological models, Water Resour. Res. 44 p. 91-94.
#' Clark M. P., McMillan H. K., Collins D. B. G., Kavetski D. and Woods R. A. (2011), Hydrological field data from a modeller's perspective: Part 2: process-based evaluation of model hypotheses. Hydrological Processes, 25: 523-543. doi: 10.1002/hyp.7902.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(fuse_hydrological_timeseries)
#' set.seed(123)
#' parameters <- generateParameters(1)
#' U <- fusesma.sim(fuse_hydrological_timeseries,
#'                  60, 1, parameters$rferr_add, parameters$rferr_mlt,
#'                  parameters$frchzne, parameters$fracten,
#'                  parameters$maxwatr_1, parameters$percfrac,
#'                  parameters$fprimqb, parameters$qbrate_2a,
#'                  parameters$qbrate_2b, parameters$qb_prms,
#'                  parameters$maxwatr_2, parameters$baserte,
#'                  parameters$rtfrac1, parameters$percrte, parameters$percexp,
#'                  parameters$sacpmlt, parameters$sacpexp, parameters$iflwrte,
#'                  parameters$axv_bexp, parameters$sareamax,
#'                  parameters$loglamb, parameters$tishape, parameters$qb_powr)
#' Q <- fuserouting.sim(U, 60, 1, 0.5)
#' }
#'

fuserouting.sim <- function(U, mid, deltim, timedelay) {

   inAttr <- attributes(U)
   U <- zoo::coredata(U)

   # list of availabe models is loaded by LazyData: true in DESCRIPTION
   modlist <- modlist                            # to remove NOTE in R CMD check

   # Read model structure
   smodl <- c(modlist[mid,2],    # rferr
              modlist[mid,3],    # arch1
              modlist[mid,4],    # arch2
              modlist[mid,5],    # qsurf
              modlist[mid,6],    # qperc
              modlist[mid,7],    # esoil
              modlist[mid,8],    # qintf
              modlist[mid,9])    # q_tdh

   res <- .C(  "fuseRoutingSimR",
               as.double(U),
               as.integer(length(U)),
               as.integer(smodl),      # A (1d) array of ints (the model params)
               as.integer(length(smodl)),
               as.double(timedelay),
               as.double(deltim),
               X = double( length(U) ),
               as.integer(length(U)),
               status = integer(1), PACKAGE="fuse" )

   if (res$status != 0){
     stop("fuseRoutingSimR Failed!")
   }

   X <- res$X
   attributes(X) <- inAttr

   return(X)
}


#' Function to define the parameter ranges for GAMMA routing module
#'
#' @description Parameter ranges for GAMMA routing module, as suggested in Clark et al. (2011).
#'
#' @return The function returns the list of parameter ranges for the Routing model.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fuserouting.ranges()
#' }
#'

fuserouting.ranges <- function() {
   list("timedelay" = c(0.01, 5))        # time delay in runoff (days)
}


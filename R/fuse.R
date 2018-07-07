#' Combine FUSE soil moisture accounting module and routing module in one function
#'
#' @param DATA data.frame containing observations. It consists of 3 columns:
#' Rainfall (P), Potential Evapo-Transpiration (E) and Streamflow (Q)
#' @param mid model id number in Model List 2011(see below for details)
#' @param deltim observation time step (days)
#' @param ParameterSet list of parameters
#'
#' @details The list of parameters can be generated as follows: ParameterSet <- GeneratePsetsFUSE(1).
#'
#' @return Simulated streamflow discharge
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(fuse_hydrological_timeseries)
#' set.seed(123)
#' parameters <- generateParameters(1)
#' Q <- fuse(fuse_hydrological_timeseries, 60, 1, parameters)
#' }
#'

fuse <- function(DATA, mid, deltim, ParameterSet){
  
  # Check parameter names are correct
  pnames <- c("rferr_add", "rferr_mlt", "maxwatr_1", "maxwatr_2", "fracten",
              "frchzne", "fprimqb", "rtfrac1", "percrte", "percexp", "sacpmlt",
              "sacpexp", "percfrac", "iflwrte", "baserte", "qb_powr", "qb_prms",
              "qbrate_2a", "qbrate_2b", "sareamax", "axv_bexp", "loglamb",
              "tishape", "timedelay")
  if (all(names(ParameterSet) != pnames)){
    stop("The parameter set has incorrect names, please correct and re-try.")
  }

  U <- fusesma.sim(DATA,
                   mid,
                   deltim,
                   ParameterSet$rferr_add,
                   ParameterSet$rferr_mlt,
                   ParameterSet$frchzne,
                   ParameterSet$fracten,
                   ParameterSet$maxwatr_1,
                   ParameterSet$percfrac,
                   ParameterSet$fprimqb,
                   ParameterSet$qbrate_2a,
                   ParameterSet$qbrate_2b,
                   ParameterSet$qb_prms,
                   ParameterSet$maxwatr_2,
                   ParameterSet$baserte,
                   ParameterSet$rtfrac1,
                   ParameterSet$percrte,
                   ParameterSet$percexp,
                   ParameterSet$sacpmlt,
                   ParameterSet$sacpexp,
                   ParameterSet$iflwrte,
                   ParameterSet$axv_bexp,
                   ParameterSet$sareamax,
                   ParameterSet$loglamb,
                   ParameterSet$tishape,
                   ParameterSet$qb_powr)

  Q <- fuserouting.sim(U,
                       mid,
                       deltim,
                       timedelay=ParameterSet$timedelay)

  return(Q)

}

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
#' data(DATA)
#' set.seed(123)
#' parameters <- generateParameters(1)
#' Q <- fuse(DATA, 60, 1, parameters)
#' }
#'

fuse <- function(DATA, mid, deltim, ParameterSet){

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

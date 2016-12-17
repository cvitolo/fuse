#' Generate information on model structure components and parameter values
#'
#' @description fuseInfo is used to generate tables to summarise which model component and parameters are used given the model id. This is important for calculating the posterior distribution of parameter values and model structures.
#'
#' @param mid FUSE model structure ID number (integer from 1 to 1248).
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   fuseInfo(mid=60)
#' }
#'

fuseInfo <- function(mid) {

  # Model structure info #######################################################
  model <- data.frame(matrix(readmd2var(mid,number=TRUE),ncol=8,nrow=1))
  names(model) <- names(readmd2var(mid,number=TRUE))

  # Parameters info ############################################################
  parameters <- data.frame(matrix(FALSE,ncol=24,nrow=1))
  names(parameters) <- c("rferr_add","rferr_mlt","maxwatr_1","maxwatr_2",
                         "fracten","frchzne","fprimqb","rtfrac1","percrte",
                         "percexp","sacpmlt","sacpexp","percfrac","iflwrte",
                         "baserte","qb_powr","qb_prms","qbrate_2a",
                         "qbrate_2b","sareamax","axv_bexp","loglamb",
                         "tishape","timedelay")

  # (1) rainfall errors ########################################################
  if(model$rferr == 11) parameters$rferr_add <- TRUE                # additive_e
  if(model$rferr == 12) parameters$rferr_mlt <- TRUE                # multiplc_e

  # (2) upper-layer architecture ###############################################
  # Common parameters:
  parameters$maxwatr_1 <- TRUE            # maximum total storage in layer1 (mm)
  parameters$fracten   <- TRUE       # frac total storage as tension storage (-)
  # even with a single state, tension and free storage should be defined!

  # Model-specific parameters:

  # onestate_1, no other parameters needed
  # if(model$arch1 == 21) {}

  # tension1_1, no other parameters needed
  # if(model$arch1 == 22) {}

  # tension2_1 = tension storage sub-divided into recharge and excess
  if(model$arch1 == 23) {

    # PRMS: frac tension storage in recharge zone (-)
    parameters$frchzne <- TRUE

  }

  # (3) lower-layer architecture / baseflow ####################################
  # Common parameters:
  # maximum total storage in layer2 (mm)
  parameters$maxwatr_2 <- TRUE
  # The params below are always used to calculate qbsat,
  # this also means that loglamb and tishape are always calculated!
  # mean value of the log-transformed topographic index (m)
  parameters$loglamb   <- TRUE
  # shape parameter for the topo index Gamma distribution (-)
  parameters$tishape   <- TRUE
  # baseflow exponent (-)
  parameters$qb_powr   <- TRUE

  # Model-specific parameters:

  # fixedsiz_2
  # power-law relation, no parameters needed for the topo index distribution
  if(model$arch2 == 31) {

    # baseflow rate (mm day-1)
    parameters$baserte   <- TRUE

    # baseflow exponent (-)
    # parameters$qb_powr   <- TRUE

  }

  # tens2pll_2 = tension reservoir plus two parallel tanks
  if(model$arch2 == 32) {

    # fraction of percolation to tension storage (-)
    parameters$percfrac  <- TRUE

    # SAC: fraction of baseflow in primary resvr (-)
    parameters$fprimqb   <- TRUE

    # baseflow depletion rate for primary resvr (day-1)
    parameters$qbrate_2a <- TRUE

    # baseflow depletion rate for secondary resvr (day-1)
    parameters$qbrate_2b <- TRUE

  }

  # unlimfrc_2 = baseflow resvr of unlimited size (0-huge), frac rate
  if(model$arch2 == 33) {

    # baseflow depletion rate (day-1)
    parameters$qb_prms   <- TRUE

  }

  # unlimpow_2 = topmodel option, power-law transmissivity profile
  if(model$arch2 == 34) {

    # baseflow rate (mm day-1)
    parameters$baserte   <- TRUE

    # mean value of the log-transformed topographic index (m)
    # parameters$loglamb   <- TRUE

    # shape parameter for the topo index Gamma distribution (-)
    # parameters$tishape   <- TRUE

    # baseflow exponent (-)
    # parameters$qb_powr   <- TRUE

  }

  # (4) surface runoff #########################################################
  if(model$qsurf == 41) {

    # arno_x_vic = arno/xzang/vic parameterization (upper zone control)
    parameters$axv_bexp <- TRUE

  }

  if(model$qsurf == 42) {

    # prms_varnt = prms variant (fraction of upper tension storage)
    parameters$sareamax <- TRUE

  }

  # if(model$qsurf == 43) {
  #
  #   # tmdl_param = topmodel parameterization
  #   # need the topographic index if we don't have it for baseflow
  #   # if(model$arch2 == 32 || model$arch2 == 33 || model$arch2 == 31) {
  #   #   parameters$loglamb   <- TRUE
  #   #   parameters$tishape   <- TRUE
  #   #}
  #
  #   # if(model$arch2 == 32 || model$arch2 == 33 || model$arch2 == 35) {
  #   #   # need the topmodel power if we don't have it for baseflow
  #   #   # baseflow exponent (-), used to modify the topographic
  #   #   parameters$qb_powr   <- TRUE
  #   # }
  #
  # }

  # (5) percolation ############################################################

  # perc_f2sat, perc_w2sat = standard equation k(theta)**c
  if(model$qperc == 51||model$qperc == 53) {

    # percolation rate (mm day-1)
    parameters$percrte   <- TRUE

    # percolation exponent (-)
    parameters$percexp   <- TRUE

  }

  # perc_lower = perc defined by moisture content in lower layer (sac)
  if(model$qperc == 52) {

    # multiplier in the SAC model for dry lower layer (-)
    parameters$sacpmlt   <- TRUE

    # exponent in the SAC model for dry lower layer (-)
    parameters$sacpexp   <- TRUE
  }

  # (6) evaporation ############################################################

  # sequential, no additional parameters
  # if(model$esoil == 62) {}

  # rootweight
  if(model$esoil == 61) {

    # fraction of roots in the upper layer (-)
    parameters$rtfrac1 <- TRUE

  }

  # (7) interflow ##############################################################

  # intflwnone, no additional parameters
  # if(model$qintf == 71) {}

  # intflwsome
  if(model$qintf == 72) {

    # interflow rate (mm day-1)
    parameters$iflwrte <- TRUE

  }

  # (8) routing ################################################################

  # no_routing, no additional parameters
  # if(model$q_tdh == 81) {}

  # routgamma
  if(model$q_tdh == 82) {

    # timedelay (days)
    parameters$timedelay <- TRUE

  }

  return(cbind(mid, model, parameters))

}

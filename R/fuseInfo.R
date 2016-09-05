#' This function returns information on model structure components and parameters used, given FUSE model.
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

  # (1) rainfall errors
  if(model$rferr == 11) parameters$rferr_add <- TRUE                # additive_e
  if(model$rferr == 12) parameters$rferr_mlt <- TRUE                # multiplc_e

  # (2) upper-layer architecture
  # Common parameters:
  parameters$maxwatr_1 <- TRUE # maximum total storage in layer1 (mm)
  parameters$fracten   <- TRUE # frac total storage as tension storage (-)
  # even with a single state, tension and free storage should be defined!
  # Model-specific parameters:
  # if(model$arch1 == 21) {}            # onestate_1, no other parameters needed
  # if(model$arch1 == 22) {}            # tension1_1, no other parameters needed
  if(model$arch1 == 23) {
    # tension2_1 = tension storage sub-divided into recharge and excess
    parameters$frchzne <- TRUE # PRMS: frac tension storage in recharge zone (-)
  }

  # (3) lower-layer architecture / baseflow
  # Common parameters:
  parameters$maxwatr_2 <- TRUE # maximum total storage in layer2 (mm)
  # Model-specific parameters:
  if(model$arch2 == 31) {
    # fixedsiz_2
    # power-law relation, no parameters needed for the topo index distribution
    parameters$baserte   <- TRUE  # baseflow rate (mm day-1)
    parameters$qb_powr   <- TRUE  # baseflow exponent (-)
  }
  if(model$arch2 == 32) {
    # tens2pll_2 = tension reservoir plus two parallel tanks
    parameters$percfrac  <- TRUE  # fraction of percolation to tension storage (-)
    parameters$fprimqb   <- TRUE  # SAC: fraction of baseflow in primary resvr (-)
    parameters$qbrate_2a <- TRUE  # baseflow depletion rate for primary resvr (day-1)
    parameters$qbrate_2b <- TRUE  # baseflow depletion rate for secondary resvr (day-1)
  }
  if(model$arch2 == 33) {
    # unlimfrc_2 = baseflow resvr of unlimited size (0-huge), frac rate
    parameters$qb_prms   <- TRUE  # baseflow depletion rate (day-1)
  }
  if(model$arch2 == 34) {
    # unlimpow_2 = topmodel option, power-law transmissivity profile
    parameters$baserte   <- TRUE  # baseflow rate (mm day-1)
    parameters$loglamb   <- TRUE  # mean value of the log-transformed topographic index (m)
    parameters$tishape   <- TRUE  # shape parameter for the topo index Gamma distribution (-)
    parameters$qb_powr   <- TRUE  # baseflow exponent (-)
  }

  # (4) surface runoff
  if(model$qsurf == 41) {
    # arno_x_vic = arno/xzang/vic parameterization (upper zone control)
    parameters$axv_bexp <- TRUE
  }
  if(model$qsurf == 42) {
    # prms_varnt = prms variant (fraction of upper tension storage)
    parameters$sareamax <- TRUE
  }
  if(model$qsurf == 43) {
    # tmdl_param = topmodel parameterization
    # need the topographic index if we don't have it for baseflow
    if(model$arch2 == 32 || model$arch2 == 33 || model$arch2 == 31) {
      parameters$loglamb   <- TRUE
      parameters$tishape   <- TRUE
    }
    if(model$arch2 == 32 || model$arch2 == 33 || model$arch2 == 35) {
      # need the topmodel power if we don't have it for baseflow
      # baseflow exponent (-), used to modify the topographic
      parameters$qb_powr   <- TRUE
    }
  }

  # (5) percolation
  # perc_f2sat, perc_w2sat = standard equation k(theta)**c
  if(model$qperc == 51||model$qperc == 53) {
    parameters$percrte   <- TRUE             # percolation rate (mm day-1)
    parameters$percexp   <- TRUE             # percolation exponent (-)
  }

  # perc_lower = perc defined by moisture content in lower layer (sac)
  if(model$qperc == 52) {
    parameters$sacpmlt   <- TRUE             # multiplier in the SAC model for dry lower layer (-)
    parameters$sacpexp   <- TRUE             # exponent in the SAC model for dry lower layer (-)
  }

  # (6) evaporation
  # if(model$esoil == 62) {}              # sequential, no additional parameters
  if(model$esoil == 61) {
    # rootweight
    parameters$rtfrac1 <- TRUE   # fraction of roots in the upper layer (-)
  }

  # (7) interflow
  # if(model$qintf == 71) {}              # intflwnone, no additional parameters
  if(model$qintf == 72) {
    parameters$iflwrte <- TRUE   # interflow rate (mm day-1)
  }

  # (8) routing
  # if(model$q_tdh == 81) {}              # no_routing, no additional parameters
  if(model$q_tdh == 82) {
    parameters$timedelay <- TRUE   # timedelay (days)
  }

  return(cbind(model,parameters))

}

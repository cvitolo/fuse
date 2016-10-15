#' Prepare a table containing Model ID numbers, model options and parameter in use.
#'
#' @param reduced This is a logical value. If FALSE the entire model list is taken into account. If FALSE, a reduced version is selected (78 models, according to Clark et al. (2008)).
#'
#' @return table containing Model ID numbers, model options and parameter in use
#'
#' @export
#'
#' @examples
#' # model_param_table()
#'

model_param_table <- function(reduced = FALSE){

  modlist <- modlist

  if (reduced == TRUE) {

    # Rainfall error is not included in the inference
    # By convention, we keep only models with rferr = 12.
    modlist <- modlist[which(modlist[,"rferr"] == 12),]

    # Routing is always allowed
    # therefore we keep only models with q_tdh = 82.
    modlist <- modlist[which(modlist[,"q_tdh"] == 82),]

    # The effect of the evaporation scheme is negligible
    modlist <- modlist[-which(modlist[,"esoil"] == 61),]

    # The effect of the interflow is negligible
    modlist <- modlist[-which(modlist[,"qintf"] == 72),]

    # The number of models to take into consideration reduces to 78 (79?)
    # according to Clark et al. (2008).

    # Later, I found there are some model combinations which are
    # physically unrealistic, therefore could be removed:
    # modlist <- modlist[-which(modlist[,"arch2"] == 33 &
    #                           modlist[,"qperc"] == 52),]
    # modlist <- modlist[-which(modlist[,"arch2"] == 34 &
    #                           modlist[,"qperc"] == 52),]
    # modlist <- modlist[-which((modlist[,"arch2"] == 33 |
    #                              modlist[,"arch2"] == 34) &
    #                             modlist[,"qperc"] == 52),]

  }

  tempList <- lapply(X = modlist$mid, FUN = fuseInfo)

  extendedModlist <- dplyr::bind_rows(tempList)

  return(extendedModlist)

}

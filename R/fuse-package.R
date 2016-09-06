#' fuse: Framework for Understanding Structural Errors
#'
#' Implementation of the framework for hydrological modelling FUSE, based on the Fortran version described in Clark et al. (2008). The package consists of two modules: Soil Moisture Accounting module (\code{fuse.sim}) and Gamma routing module (\code{routgamma.sim}). It also contains default parameter ranges (see \code{fuse.ranges} and \code{routgamma.ranges}) and three data objects: DATA, parameters and modlist.
#'
#' @name fuse
#' @docType package
#'
#' @importFrom zoo coredata
#' @importFrom tgp lhs
#' @importFrom stats runif
#' @importFrom utils data
#'
#' @useDynLib fuse
#'
NULL

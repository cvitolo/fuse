/*! \mainpage Introduction
 *
 * FUSE hydrological framework is illustrated in Clark et al. (2008). It is an ensemble of four parent lumped models: PRMS, SACRAMENTO, TOPMODEL and ARNO/VIC. Each model is characterized by a different architecture of the upper and lower soil layers and for the parametrization of processes such as: evaporation, vertical percolation between soil layers, interflow, base flow and surface runoff. 
 *
 * FUSE can combine each element of one model with elements from other models to obtain several model structures. Model structures are selected through 8 options:
 * 
 * -# Rainfall error, which can be:
 *   -# additive
 *   -# multiplicative
 *   .
 * 
 * -# Architecture of the upper soil layer, which can be:
 *   -# upper layer defined by a single state variable
 *   -# upper layer broken up into tension and free storage
 *   -# upper layer defined by a tension storage sub-divided into recharge and excess
 *   .
 *
 * -# Architecture of the lower soil layer, which can be:
 *   -# lower layer defined by a single state variable, with a baseflow reservoir of fixed size
 *   -# lower layer defined by a single state variable, with a baseflow reservoir of unlimited size (frac rate)
 *   -# lower layer defined by a single state variable, with a baseflow reservoir of unlimited size (power recession)
 *   -# lower layer defined by a tension reservoir plus two parallel tanks (in this case the upper soil layer can only be like 2.a or 2.b)
 *   .
 *
 * -# Surface runoff scheme, which is based only on the saturation-excess mechanism. The saturated area is calsulated using:
 *   -# arno/xzang/vic parameterization (upper zone control)
 *   -# prms variant (fraction of upper tension storage)
 *   -# topmodel parameterization
 *   .
 *
 * 
 * -# Percolation scheme, in which water availability for percolation is:
 *   -# limited by the field capacity
 *   -# limited by the wilting point (in this case the upper soil layer can only be like 2.a)
 *   -# defined by moisture content in lower layer (sac)
 *   .
 * 
 * -# Evaporation scheme, which can be:
 *   -# sequential (the potential evaporative demand is first satisfied by evaporation from upper soil layer, then the residual from the lower soil layer)
 *   -# rootweight (evaporation is computed on the basis of the relative root fractions in each of the soil layers)
 *   . 
 * 
 * -# Interflow scheme, which can be:
 *   -# allowed
 *   -# denied (e.g. TOPMODEL and ARNO/VIC)
 *   .
 * 
 * -# Time delay in run-off, which can be:
 *   -# allowed (it uses a two-parameters gamma distribution to route runoff to the basin outlet)
 *   -# denied
 *   .
 * .
 * 
 * The original Fortran code was provided by Martyn Clark, then ported in R as part of the RHydro package (Reusser et al., 2012). The current library implements FUSE's core functions in C++ to speed up processing time.
 *
 * \section install_sec Reference
 *
 * Clark, M. P., A. G. Slater, D. E. Rupp, R. A. Woods, J. A. Vrugt, H. V. Gupta, T. Wagener, and L. E. Hay (2008),
    Framework for Understanding Structural Errors (FUSE): A modular framework to diagnose differences between hydrological models,
    Water Resour. Res., 44, W00B02, doi:10.1029/2007WR006735.
 * 
 * Reusser D. E., Buytaert W., Vitolo, C. (2012). RHydro – Hydrological models and tools to represent and analyze hydrological data in R. European Geosciences Union General Assembly 2012, Vienna, Austria.
 */

#ifndef FUSE_H
#define FUSE_H

/*! \file fuse.h
\brief The main fuse library.
fuse.h contains the top-level routing functions made available by the library.  
*/

#include "Fuse_global.h"
#include "smaStructures.h"
#include "modelStructure.h"

#include <iostream>

using namespace std;

typedef std::vector< double > state_type;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////                                   ///////////////////////////////////////
///////////////////////////////  Externally Accessible Functions  ///////////////////////////////////////
///////////////////////////////                                   ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/*! \brief Version.
    Output version information to stdout.
 */
FUSESHARED_EXPORT void version(){ cout << "1.0.4"; }

/*! \brief Version.
    Output version information to stdout.
 */
R_CALL void CALL_CONV versionR(){ version(); }

/*! \brief FUSE Routing model.
    \param instantaneousDischarge an array of effective rainfall (sum of surface runoff, overflow, interflow, and baseflow).
    \param instantaneousDischargeSize the size of instantaneousDischarge.
    \param modelId the index of the model configuration to use for routing.
    \param modelListFileName filename containing Model Structures.
    \param timeDelay the time delay in runoff (days).
    \param timeStep the data time step measured in days: timeStep = 1 (daily time step), 1/24 (hourly time step), 1/24/4 (15 min time step).
    \param routedDischarge a pointer to the array containing the resultant routed runoff.
    \param routedDischargeSize the length of routedDischarge.
    \return returns true on success, false on failure.
    \warning modelId is the index of a Zero-based array.
    \warning Note that if the '~' character (or similar) is used, it will need to be expanded before it
             is passed to this function.
    \warning routedDischarge needs to be allocated by the caller.

    The time delay in runoff is modeled using a two-parameter Gamma distribution with shape parameter = 2.5 and scale parameter = 1.

    Reference: Clark, M. P., A. G. Slater, D. E. Rupp, R. A. Woods, J. A. Vrugt, H. V. Gupta, T. Wagener, and L. E. Hay (2008),
    Framework for Understanding Structural Errors (FUSE): A modular framework to diagnose differences between hydrological models,
    Water Resour. Res., 44, W00B02, doi:10.1029/2007WR006735.

 */
FUSESHARED_EXPORT bool fuseRoutingSim(double* instantaneousDischarge,
                                      int instantaneousDischargeSize,        // Length of instantaneousDischarge
                                      int* modelParams,
                                      int modelParamsSize,
                                      double timeDelay,
                                      double timeStep,
                                      double* routedDischarge,                        // Result is stored here
                                      int routedDischargeSize);                        // length of routedDischarge




/*! \brief R Compatible FUSE Routing model.
    \param instantaneousDischarge an array of effective rainfall (sum of surface runoff, overflow, interflow, and baseflow).
    \param instantaneousDischargeSize the size of instantaneousDischarge.
    \param modelId the index of the model configuration to use for routing.
    \param modelListFileName filename containing Model Structures.
    \param timeDelay the time delay in runoff (days).
    \param timeStep the data time step measured in days: timeStep = 1 (daily time step), 1/24 (hourly time step), 1/24/4 (15 min time step).
    \param routedDischarge a pointer to the array containing the resultant routed runoff.
    \param routedDischargeSize the length of routedDischarge.
    \param status reports the status of the call.  Zero = Success, Non-Zero = Failure.
    \warning modelId is the index of a One-based array.
    \warning Note that if the '~' character (or similar) is used, it will need to be expanded before it
             is passed to this function.

    The time delay in runoff is modeled using a two-parameter Gamma distribution with shape parameter = 2.5 and scale parameter = 1. 

    Reference: Clark, M. P., A. G. Slater, D. E. Rupp, R. A. Woods, J. A. Vrugt, H. V. Gupta, T. Wagener, and L. E. Hay (2008),
    Framework for Understanding Structural Errors (FUSE): A modular framework to diagnose differences between hydrological models,
    Water Resour. Res., 44, W00B02, doi:10.1029/2007WR006735.

    A typical call to this function from R might look like:

    \code
    fuseRoutingSim <- function( instantaneousDischarge,
                                modelId,
                                modelListFileName,
                                timeDelay,
                                timeStep ){
  
        res <- .C(  "fuseRoutingSimR",
                    as.double(instantaneousDischarge),
                    as.integer(length(instantaneousDischarge)),
                    as.integer(modelId),
                    as.character(modelListFileName),
                    as.double(timeDelay),
                    as.double(timeStep),
                    routedDischarge = double( length(instantaneousDischarge) ),
                    as.integer(length(instantaneousDischarge)),
                    status = integer(1) )
  
        if (res$status != 0){
            stop("fuseRoutingSimR Failed!")
        }
  
        return(res$routedDischarge)
  
    }
    \endcode

 */
R_CALL void CALL_CONV fuseRoutingSimR(  double* instantaneousDischarge,
                                        int* instantaneousDischargeSize,     // Length of instantaneousDischarge
                                        int* modelParams,
                                        int* modelParamsSize, 
                                        double* timeDelay,
                                        double* timeStep,
                                        double* routedDischarge,                      // Result is returned here
                                        int* routedDischargeSize,                      // length of routedDischarge
                                        int* status);                  // Status (true=success, false=fail) returned here




/*! \brief FUSE Soil Moisture Accounting model.
    \param precipitation an array of precipitation
    \param precipitationLength the length of precipitation
    \param potentialEvapotrans an array of potential evapotranspiration
    \param potentialEvapotransLength the length of potentialEvapotrans
    \param modelId the model id
    \param modelListFileName filename containing model structure definitions.
    \param timeStep the data time step measured in days: timeStep = 1 (daily time step), 1/24 (hourly time step), 1/24/4 (15 min time step).
    \param frchzne frchzne  input model parameter.
    \param fracten fracten input model parameter.
    \param maxwatr_1 maxwatr_1 input model parameter.
    \param percfrac percfrac input model parameter.
    \param fprimqb fprimqb input model parameter.
    \param qbrate_2a qbrate_2a input model parameter.
    \param qbrate_2b qbrate_2b input model parameter.
    \param qb_prms qb_prms input model parameter.
    \param maxwatr_2 maxwatr_2 input model parameter.
    \param baserte baserte input model parameter.
    \param rtfrac1 rtfrac1 input model parameter.
    \param percrte percrte input model parameter.
    \param percexp percexp input model parameter.
    \param sacpmlt sacpmlt input model parameter.
    \param sacpexp sacpexp input model parameter.
    \param iflwrte iflwrte input model parameter.
    \param axv_bexp axv_bexp input model parameter.
    \param sareamax sareamax input model parameter.
    \param loglamb loglamb input model parameter.
    \param tishape tishape input model parameter.
    \param qb_powr qb_powr input model parameter.
    \param states pointer to an allocated 1-dimensional array of length (10 * precipitationLength) in which to store output states.
    \param fluxes pointer to an allocated 1-dimensional array of length (10 * precipitationLength) in which to store output fluxes.
    \param abs_error solver's absolute error (default value 1.0e-2).
    \param rel_error solver's relative error (default value 1.0e-2).
    \param correctNegativeVolumes whether or not distribute negative volume over non-zero volumes (on by default).
    \param fracstate0 fracstate0 input parameter, default value 0.25.
    \param rferr_add rferr_add input parameter, default value 0.0
    \param rferr_mlt rferr_mlt input parameter, default value 1.0
    \return true on success and false otherwise.

    Reference: Clark, M. P., A. G. Slater, D. E. Rupp, R. A. Woods, J. A. Vrugt, H. V. Gupta, T. Wagener, and L. E. Hay (2008), 
               Framework for Understanding Structural Errors (FUSE): A modular framework to diagnose differences between hydrological models,
               Water Resour. Res., 44, W00B02, doi:10.1029/2007WR006735.

 */
FUSESHARED_EXPORT bool fusesmaSim(  double* precipitation,
                                    int precipitationLength,
                                    double* potentialEvapotrans,
                                    int potentialEvapotransLength,
                                    int* modelParams,
                                    int modelParamsSize,
                                    double timeStep,
                                    double frchzne,
                                    double fracten,
                                    double maxwatr_1,
                                    double percfrac,
                                    double fprimqb,
                                    double qbrate_2a,
                                    double qbrate_2b,
                                    double qb_prms,
                                    double maxwatr_2,
                                    double baserte,
                                    double rtfrac1,
                                    double percrte,
                                    double percexp,
                                    double sacpmlt,
                                    double sacpexp,
                                    double iflwrte,
                                    double axv_bexp,
                                    double sareamax,
                                    double loglamb,
                                    double tishape,
                                    double qb_powr,
                                    double* states,
                                    double* fluxes,
                                    double abs_error=1.0e-2,
                                    double rel_error=1.0e-2,
                                    bool correctNegativeVolumes=true,
                                    double fracstate0=0.25,
                                    double rferr_add=0.0,
                                    double rferr_mlt=1.0,
                                    bool outputStates=true,
                                    bool outputFluxes=true,
                                    bool debug=false);




/*! \brief R Compatible FUSE Soil Moisture Accounting model.
    \param precipitation an array of precipitation
    \param precipitationLength the length of precipitation
    \param potentialEvapotrans an array of potential evapotranspiration
    \param potentialEvapotransLength the length of potentialEvapotrans
    \param modelId the model id
    \param modelListFileName filename containing model structure definitions.
    \param timeStep the data time step measured in days: timeStep = 1 (daily time step), 1/24 (hourly time step), 1/24/4 (15 min time step).
    \param frchzne frchzne  input model parameter.
    \param fracten fracten input model parameter.
    \param maxwatr_1 maxwatr_1 input model parameter.
    \param percfrac percfrac input model parameter.
    \param fprimqb fprimqb input model parameter.
    \param qbrate_2a qbrate_2a input model parameter.
    \param qbrate_2b qbrate_2b input model parameter.
    \param qb_prms qb_prms input model parameter.
    \param maxwatr_2 maxwatr_2 input model parameter.
    \param baserte baserte input model parameter.
    \param rtfrac1 rtfrac1 input model parameter.
    \param percrte percrte input model parameter.
    \param percexp percexp input model parameter.
    \param sacpmlt sacpmlt input model parameter.
    \param sacpexp sacpexp input model parameter.
    \param iflwrte iflwrte input model parameter.
    \param axv_bexp axv_bexp input model parameter.
    \param sareamax sareamax input model parameter.
    \param loglamb loglamb input model parameter.
    \param tishape tishape input model parameter.
    \param qb_powr qb_powr input model parameter.
    \param s pointer to an allocated 1-dimensional array of length (10 * precipitationLength) in which to store output states.
    \param f pointer to an allocated 1-dimensional array of length (20 * precipitationLength) in which to store output fluxes.
    \param abs_error solver's absolute error (default value 1.0e-2).
    \param rel_error solver's relative error (default value 1.0e-2).
    \param correctNegativeVolumes whether or not distribute negative volume over non-zero volumes (on by default).
    \param fracstate0 fracstate0 input parameter, default value 0.25.
    \param rferr_add rferr_add input parameter, default value 0.0
    \param rferr_mlt rferr_mlt input parameter, default value 1.0
    \param status reports the status of the call.  Zero = Success, Non-Zero = Failure.

    A typical call to this function from R might look like:

    \code
    fusesmaSimRLocal <- function(   precipitation,
                                    potentialEvapotrans,
                                    modelId,
                                    modelListFileName,
                                    timeStep,
                                    frchzne,
                                    fracten,
                                    maxwatr_1,
                                    percfrac,
                                    fprimqb,
                                    qbrate_2a,
                                    qbrate_2b,
                                    qb_prms,
                                    maxwatr_2,
                                    baserte,
                                    rtfrac1,
                                    percrte,
                                    percexp,
                                    sacpmlt,
                                    sacpexp,
                                    iflwrte,
                                    axv_bexp,
                                    sareamax,
                                    loglamb,
                                    tishape,
                                    qb_powr,
                                    absError,
                                    relError,
                                    correctNegVols,
                                    fracstate0,
                                    rferr_add,
                                    rferr_mlt ){

        # load the library
        # Notice the lack of .so file extension below
        dyn.load(paste("/full/path/to/libFuse", .Platform$dynlib.ext, sep = "")) # general statement

        res <- .C(  "fusesmaSimR",
                    as.double(precipitation),
                    as.integer(length(precipitation)),
                    as.double(potentialEvapotrans),
                    as.integer(length(potentialEvapotrans)),
                    as.integer(modelId),
                    as.character(modelListFileName),
                    as.double(timeStep),
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
                    stateResults = double( 10 * length(precipitation) ),
                    fluxResults = double( 20 * length(precipitation) ),
                    as.double(absError),
                    as.double(relError),
                    as.logical(correctNegVols),
                    as.double(fracstate0),
                    as.double(rferr_add),
                    as.double(rferr_mlt),
                    status = integer(1) )

        if (res$status != 0){
            stop("fusesmaSimR Failed!")
        }

        results <- list("states"=res$stateResults,"fluxes"=res$fluxResults)

        return(results)

    }
    \endcode

 */
R_CALL void CALL_CONV fusesmaSimR(    double* precipitation,
                            int* precipitationLength,
                            double* potentialEvapotrans,
                            int* potentialEvapotransLength,
                            int* modelParams,
                            int* modelParamsSize,
                            double* timeStep,
                            double* frchzne,
                            double* fracten,
                            double* maxwatr_1,
                            double* percfrac,
                            double* fprimqb,
                            double* qbrate_2a,
                            double* qbrate_2b,
                            double* qb_prms,
                            double* maxwatr_2,
                            double* baserte,
                            double* rtfrac1,
                            double* percrte,
                            double* percexp,
                            double* sacpmlt,
                            double* sacpexp,
                            double* iflwrte,
                            double* axv_bexp,
                            double* sareamax,
                            double* loglamb,
                            double* tishape,
                            double* qb_powr,
                            double* states,
                            double* fluxes,
                            double* abs_error,
                            double* rel_error,
                            bool* correctNegativeVolumes,
                            double* fracstate0,
                            double* rferr_add,
                            double* rferr_mlt,
                            int* status);                       // Status (true=success, false=fail) returned here




/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////                                   ///////////////////////////////////////
///////////////////////////////        Internal Functions         ///////////////////////////////////////
///////////////////////////////                                   ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////




/*! \brief Computes baseflow at saturation.
    \param lowerSoilLayerArch architecture of the lower soil layer, it can be made of: one single state (fixed or unlimited size) or a tension storage plus two parallel base flow reservoirs.
    \param inputParams input model parameters.
    \param maximumStorageInThePrimaryBaseflowReservoir maximum storage in the primary base flow reservoir (derived parameter).
    \param maximumStorageInTheSecondaryBaseflowReservoir maximum storage in the secondary base flow reservoir (derived parameter).
    \param meanOfThePowerTransformedTopographicIndex mean of the power-transformed topographic index (derived parameter).
    \return Baseflow at saturation.

     Used in the sacramento percolation model.

*/
double baseflowAtSaturation(double lowerSoilLayerArch,
                            SmaInputParameters inputParams,
                            double maximumStorageInThePrimaryBaseflowReservoir,
                            double maximumStorageInTheSecondaryBaseflowReservoir,
                            double meanOfThePowerTransformedTopographicIndex);




/*! \brief Computes derived model parameters (bucket sizes, etc.).
    \param smodl list of model components.
    \param inputParams input model parameters.
    \return List of derived parameters as an SmaDerivedParameters structure.

    See table 4 of Clark et al., 2008.

*/
SmaDerivedParameters computeDerivedModelParams(ModelStructure smodl,
                                               SmaInputParameters inputParams);




/*! \brief Compute Fluxes.
    \param smodl list of model components.
    \param precipitationAtT observed precipitation at time "t".
    \param potentialEvapotransAtT observed potential evapotranspiration at time "t".
    \param inputParams input model parameters.
    \param derivedParams derived model parameters.
    \param x model state variables at time "t".
    \return SmaFluxStructure structure of fluxes at time "t"

    See table 2 of Clark et al., 2008.

*/
SmaFluxStructure computeFluxes(ModelStructure* smodl,
                               double* precipitationAtT,
                               double* potentialEvapotransAtT,
                               SmaInputParameters* inputParams,
                               SmaDerivedParameters* derivedParams,
                               state_type* x);




/*! \brief Compute baseflow from the lower soil layer.
    \param lowerSoilLayerArch architecture of the lower soil layer, it can be made of: one single state (fixed or unlimited size) or a tension storage plus two parallel base flow reservoirs.
    \param inputParams input model parameters.
    \param baseflowAtSaturation derived model parameter.
    \param x model state variables at time "t".
    \param fractionOfBaseflowInThePrimaryBaseflowReservoir OUTPUT base flow from the primary reservoir.
    \param fractionOfBaseflowInTheSecondaryBaseflowReservoir OUTPUT base flow from the secondary reservoir.
    \param totalBaseflow OUTPUT total base flow.

    For details on the parametrization see section 4.2 of Clark et al., 2008.

*/
void fluxBaseflow(LowerSoilLayerArchType lowerSoilLayerArch,
                  SmaInputParameters* inputParams,
                  double baseflowAtSaturation,
                  state_type* x,
                  double* fractionOfBaseflowInThePrimaryBaseflowReservoir,
                  double* fractionOfBaseflowInTheSecondaryBaseflowReservoir,
                  double* totalBaseflow);




/*! \brief Compute evaporation from the upper and lower layers.
    \param upperArch architecture of the upper soil layer, it can be made of: one single state or tension (1 or 2) plus free storages.
    \param lowerSoilLayerArch architecture of the lower soil layer, it can be made of: one single state (fixed or unlimited size) or a tension storage plus two parallel base flow reservoirs.
    \param evaporationType evaporation parametrization, it can be based on: root weighting or sequential evaporation models.
    \param potentialEvapotransAtT observed potential evapotranspiration at time "t".
    \param inputParams input model parameters.
    \param derivedParams derived model parameters.
    \param x model state variables at time "t".
    \param evaporationFromTheUpperLayerRecharge OUTPUT evaporation from the primary tension store (recharge) in the upper soil layer.
    \param evaporationFromTheUpperLayerExcess OUTPUT evaporation from the secondary tension store (excess) in the upper soil layer.
    \param evaporationFromTheUpperLayer OUTPUT evaporation from the upper soil layer.
    \param evaporationFromTheLowerLayer OUTPUT evaporation from the lower soil layer.

    For details on the parametrization see section 4.3 of Clark et al., 2008.

*/
void fluxEvap(UpperSoilLayerArchType upperArch,
              LowerSoilLayerArchType lowerSoilLayerArch,
              EvaporationType evaporationType,
              double* potentialEvapotransAtT,
              SmaInputParameters* inputParams,
              SmaDerivedParameters* derivedParams,
              state_type* x,
              double* evaporationFromTheUpperLayerRecharge,
              double* evaporationFromTheUpperLayerExcess,
              double* evaporationFromTheUpperLayer,
              double* evaporationFromTheLowerLayer);




/*! \brief Compute percolation.
    \param percolationType percolation parametrization, it can be defined either by the moisture content in lower soil layer (sacramento model) or by water availability from field capacity/wilt point to saturation.
    \param inputParams input model parameters.
    \param derivedParams derived model parameters.
    \param x model state variables at time "t".
    \return Percolation.

    For details on the parametrization see section 4.4 of Clark et al., 2008.

*/
double fluxperc( PercolationType percolationType,
                 SmaInputParameters* inputParams,
                 SmaDerivedParameters* derivedParams,
                 state_type* x);




/*! \brief Computes miscellaneous fluxes.
    \param upperArch architecture of the upper soil layer, it can be made of: one single state or tension (1 or 2) plus free storages.
    \param lowerSoilLayerArch architecture of the lower soil layer, it can be made of: one single state (fixed or unlimited size) or a tension storage plus two parallel base flow reservoirs.
    \param x model state variables at time "t".
    \param inputParams input model parameters.
    \param derivedParams derived model parameters.
    \param effectiveRainfall effective precipitation.
    \param surfaceRunoff runoff.
    \param percolation percolation.
    \param flowFromRechargeToExcess pointer to object to receive flow from recharge to excess (mm day-1).
    \param flowFromTensionStorageToFreeStorageInTheUpperLayer pointer to object to receive flow from tension storage to free storage in the upper layer (mm day-1).
    \param flowFromTensionStorageToFreeStorageInTheLowerLayer pointer to object to receive flow from tension storage to free storage in the lower layer (mm day-1).
    \param overflowFromTheUpperSoilLayer pointer to object to receive overflow from the upper soil layer (mm day-1).
    \param overflowFromTheLowerSoilLayer pointer to object to receive overflow from the lower soil layer (mm day-1).
    \param overflowFromTheLowerSoilLayerPrimary pointer to object to receive overflow from the lower soil layer (mm day-1).
    \param overflowFromTheLowerSoilLayerSecondary pointer to object to receive overflow from the lower soil layer (mm day-1).
    \param optionX method to use (defaults to 1).

    fluxQMisc uses 2 methods:

    Method 1: OVERFLOW FLUXES AS A FRACTION OF INFLUXES
    Method 2: OVERFLOW FLUXES COMPUTED AS A RESIDUAL OF AVAILABLE STORAGE

    The desired method (1 or 2) is passed via the optionX parameter.

*/
void fluxQMisc(UpperSoilLayerArchType upperArch,
               LowerSoilLayerArchType lowerSoilLayerArch,
               state_type* x,
               SmaInputParameters* inputParams,
               SmaDerivedParameters* derivedParams,
               double effectiveRainfall,
               double surfaceRunoff,
               double percolation,
               double* flowFromRechargeToExcess,
               double* flowFromTensionStorageToFreeStorageInTheUpperLayer,
               double* flowFromTensionStorageToFreeStorageInTheLowerLayer,
               double* overflowFromTheUpperSoilLayer,
               double* overflowFromTheLowerSoilLayer,
               double* overflowFromTheLowerSoilLayerPrimary,
               double* overflowFromTheLowerSoilLayerSecondary,
               int optionX=1);




/*! \brief Computes the saturated area and surface runoff.
    \param runoffType runoff parametrization, it can be defined as in prms, arno or topmodel models.
    \param effectiveRainfall corrected observed rainfall, it can be: observed rainfall * multiplicative rainfall error or observed rainfall + additive rainfall error.
    \param inputParams input model parameters.
    \param derivedParams derived model parameters.
    \param x model state variables at time "t".
    \param saturatedArea pointer to allocated structure to receive Saturated Area.
    \param surfaceRunoff pointer to allocated structure to receive Surface Runoff.

    For details on the parametrization see section 4.7 of Clark et al., 2008.

*/
void fluxQSatExcess(RunoffType runoffType,
                    double effectiveRainfall,
                    SmaInputParameters* inputParams,
                    SmaDerivedParameters* derivedParams,
                    state_type* x,
                    double* saturatedArea,
                    double* surfaceRunoff);




/*! \brief Initialize model states (upper and lower soil layers).
    \param smodl list of model components.
    \param inputParams input model parameters.
    \param derivedParams derived model parameters.
    \param fracstate0 arbitrary value to set initial state variables (default = 0.25).
    \return Initial condition for state variables as an SmaStateStructure.

    Initial state variable conditions are set as 25% of the maximum storage (see Clark et al., 2008).

*/
SmaStateStructure initStates(ModelStructure smodl,
                             SmaInputParameters inputParams,
                             SmaDerivedParameters derivedParams,
                             double fracstate0);




/*! \brief Uses a logistic function to smooth the threshold at the top of a bucket.
    \param statex state variable.
    \param maxstatex maximum value of the state variable (from input/derived model parameters).
    \return Threshold at the top of a bucket.

    This function is used within fluxQMisc function.

*/
double logismooth(double stateX,
                  double maxStateX);




/*! \brief Model's Differential Equations.
    \param x model state variable at time "t".
    \param dxdt state variable variation at time "t".
    \param t time.

    Numerical method used to solve the ordinary differential equations: Runge-Kutta of the 4th order with relative/absolute tolerance of 10^(-2). 

*/
void mDifEqn( state_type &x,
              state_type &dxdt,
              double t );




/*! \brief Computes the mean of the power-transformed topographic index.
    \param shapeParameterForTheTopoIndexGammaDistribution shape parameter for the topographic index gamma distribution (input model parameter, unit = not applicable).
    \param meanValueOfTheLogTransformedTopographicIndex mean value of the log-transformed topographic index (input model parameter, unit = "m").
    \param baseflowExponent baseflow exponent (input model parameter, unit = not applicable).
    \param powerTransformedTopographicIndex power-transformed topographic index.
    \param meanOfThePowerTransformedTopographicIndex mean of the power-transformed topographic index (derived model parameter).

*/
void meanTipow(double shapeParameterForTheTopoIndexGammaDistribution,
               double meanValueOfTheLogTransformedTopographicIndex,
               double baseflowExponent,
               double* powerTransformedTopographicIndex,
               double* meanOfThePowerTransformedTopographicIndex);




/*! \brief Output internal fluxes.
    \param smodl list of model components.
    \param precipitation rain+snow melt time series.
    \param precipitationLength length of precipitation.
    \param potentialEvapotrans potential evapotranspiration time series.
    \param potentialEvapotransLength length of potentialEvapotrans.
    \param inputParams input model parameters.
    \param derivedParams derived model parameters.
    \param s input state variables.
    \param fluxes OUTPUT model fluxes.
    \param correctNegativeVolumes whether or not distribute negative volume over non-zero volumes.

    This function aims to correct negative volumes. Still under development.

*/
void outFluxes(ModelStructure* smodl,
               double* precipitation,
               int precipitationLength,
               double* potentialEvapotrans,
               int potentialEvapotransLength,
               SmaInputParameters* inputParams,
               SmaDerivedParameters* derivedParams,
               double *s,
               double *fluxes,
               bool correctNegativeVolumes);




/*! \brief Computes the fraction of runoff in future time steps.
    \param timeDelay time delay in runoff (input model parameter, unit = "days").
    \param timeStep time step in days (unit = "days").
    \param sizeFracFuture maximum number of future time steps (default = 500).
    \return returns an array of Fraction of runoff in future time steps (the size of the array is given in sizeFracFuture).

    This funtion is used within the Gamma Routing module.

*/
double* qTimeDelay(float timeDelay,
                   float timeStep,
                   int sizeFracFuture=500);



/*! \brief Update the state variables "while" solving the differential equations.
    \param x model state variables at time "t".

    This function updates the states according to the maximum size of each store.

*/
void upstates(state_type &x);




#endif // FUSE_H

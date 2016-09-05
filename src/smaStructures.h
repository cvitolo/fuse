#ifndef SMASTRUCTURES_H
#define SMASTRUCTURES_H

/*! \file smaStructures.h
\brief Data structures for SMA Input Parameters.
smaStructures.h contains various data structures for working with SMA Input Parameters.
*/




//! Model parameters.
/*! Input model parameters required by FUSE soil moisture accounting module (see table 3 of Clark et al., 2008). */
struct SmaInputParameters{
    double additiveRainfallError;                           /*!< Additive rainfall error (mm day-1, was rferr_add). */
    double multiplicativeRainfallError;                     /*!< Multiplicative rainfall error (-, was rferr_mlt). */
    double fracTensionStorageInRechargeZone;                /*!< Frac tension storage in recharge zone (used only in prms model) (-, was frchzne). */
    double fracTotalStorageAsTensionStorage;                /*!< Frac total storage as tension storage (-, was fracten). */
    double maximumTotalStorageInLayer1;                     /*!< Maximum total storage in layer1 (mm, was maxwatr_1). */
    double fractionOfPercolationToTensionStorage;           /*!< Fraction of percolation to tension storage (-, was percfrac). */
    double fractionOfBaseflowInPrimaryResvr;                /*!< Fraction of baseflow in primary resvr (used only in sac module) (-, was fprimqb). */
    double baseflowDepletionRateForPrimaryResvr;            /*!< Baseflow depletion rate for primary resvr (day-1, was qbrate_2a). */
    double baseflowDepletionRateForSecondaryResvr;          /*!< Baseflow depletion rate for secondary resvr (day-1, was qbrate_2b). */
    double baseflowDepletionRate;                           /*!< Baseflow depletion rate (day-1, was qb_prms). */
    double maximumTotalStorageInLayer2;                     /*!< Maximum total storage in layer2 (mm, was maxwatr_2). */
    double baseflowRate;                                    /*!< Baseflow rate (mm day-1, was baserte). */
    double fractionOfRootsInTheUpperLayer;                  /*!< Fraction of roots in the upper layer (-, was rtfrac1). */
    double percolationRate;                                 /*!< Percolation rate (mm day-1, was percrte). */
    double percolationExponent;                             /*!< PercolationExponent (-, was percexp). */
    double multiplierInTheSacModelForDryLowerLayer;         /*!< Multiplier in the sac model for dry lower layer (-, was sacpmlt). */
    double exponentInTheSacModelForDryLowerLayer;           /*!< Exponent in the sac model for dry lower layer (-, was sacpexp). */
    double interflowRate;                                   /*!< Interflow rate (mm day-1, was iflwrte). */
    double bExponent;                                       /*!< "b" exponent (used only in the arno/vic modules) (-, was axv_bexp). */
    double maximumSaturatedArea;                            /*!< Maximum saturated area (-, was sareamax). */
    double meanValueOfTheLogTransformedTopographicIndex;    /*!< Mean value of the log-transformed topographic index (m, was loglamb). */
    double shapeParameterForTheTopoIndexGammaDistribution;  /*!< Shape parameter for the topo index gamma distribution (-, was tishape). */
    double baseflowExponent;                                /*!< Baseflow exponent (-, was qb_powr). */
    double fractionOfSoilExcessToLowerZone;                 /*!< Fraction of soil excess to lower zone (-, was fraclowz). */
};




//! Derived model parameters.
/*! Derived model parameters, automatically calculated, in FUSE soil moisture accounting module (see table 4 of Clark et al., 2008). */
struct SmaDerivedParameters{
    double maximumTensionStorageInUpperLayer;               /*!< Maximum tension storage in upper layer (was maxtens_1). */
    double maximumTensionStorageInLowerLayer;               /*!< Maximum tension storage in lower layer (was maxtens_2). */
    double maximumFreeStorageInUpperLayer;                  /*!< Maximum free storage in upper layer (was maxfree_1). */
    double maximumFreeStorageInLowerLayer;                  /*!< Maximum free storage in lower layer (was maxfree_2). */
    double maximumStorageInThePrimaryTensionReservoir;      /*!< Maximum storage in the primary tension reservoir (recharge) (was maxtens_1a). */
    double maximumStorageInTheSecondaryTensionReservoir;    /*!< Maximum storage in the secondary tension reservoir (excess) (was maxtens_1b). */
    double maximumStorageInThePrimaryBaseflowReservoir;     /*!< Maximum storage in the primary baseflow reservoir (was maxfree_2a). */
    double maximumStorageInTheSecondaryBaseflowReservoir;   /*!< Maximum storage in the secondary baseflow reservoir (was maxfree_2b). */
    double powerTransformedTopographicIndex;                /*!< Power-transformed topographic index (was maxpow). */
    double meanOfThePowerTransformedTopographicIndex;       /*!< Mean of the power-transformed topographic index (was powlamb). */
    double baseflowAtSaturation;                            /*!< Baseflow at saturation (used in the sac percolation model) (was qbsat). */
    double fractionOfRootsInTheLowerSoilLayer;              /*!< Fraction of roots in the lower soil layer (was rtfrac2). */
};




//! Model state variables.
/*! Model state variables (units = "mm"), automatically calculated, in FUSE soil moisture accounting module (see table 1 of Clark et al., 2008). */
struct SmaStateStructure{
    double excessTensionStorageUpperLayer;          /*!< Excess Tension storage Upper Layer (Was tens_1a). */
    double rechargeTensionStorageUpperLayer;        /*!< Recharge Tension storage Upper Layer (Was tens_1b). */
    double tensionStorageUpperLayer;                /*!< Tension storage Upper Layer (Was tens_1). */
    double freeStorageUpperLayer;                   /*!< Free storage Upper Layer (Was free_1). */
    double totalUpperLayerStorage;                  /*!< Total Upper Layer Storage (Was watr_1). */
    double tensionStorageLowerLayer;                /*!< Tension storage Lower Layer (Was tens_2). */
    double freeStoragePrimaryBaseflowReservoir;     /*!< Free Storage Primary Baseflow Reservoir (Was free_2a). */
    double freeStorageSecondaryBaseflowReservoir;   /*!< Free Storage Secondary Baseflow Reservoir (Was free_2b). */
    double totalLowerLayerStorage;                  /*!< Total Lower Layer Storage (Was watr_2). */
    double freeStorageBaseflowReservoir;            /*!< Free Storage Baseflow Reservoir (Was free_2). */
};




//! Model fluxes.
/*! Model fluxes (units = "mm/d"), automatically calculated, in FUSE soil moisture accounting module (see table 2 of Clark et al., 2008). */
struct SmaFluxStructure{
    double effectiveRainfall;                                   /*!< effective rainfall (Was eff_ppt). */
    double saturatedArea;                                       /*!< saturated area (Was satarea). */
    double surfaceRunoff;                                       /*!< surface runoff (Was qrunoff). */
    double evaporationFromTheUpperLayerRecharge;                /*!< evaporation from the upper layer (recharge) (Was evap_1a). */
    double evaporationFromTheUpperLayerExcess;                  /*!< evaporation from the upper layer (excess) (Was evap_1b). */
    double evaporationFromTheUpperLayer;                        /*!< evaporation from the upper layer (Was evap_1). */
    double evaporationFromTheLowerLayer;                        /*!< evaporation from the lower layer (Was evap_2). */
    double flowFromRechargeToExcess;                            /*!< flow from recharge to excess (Was rchr2excs). */
    double flowFromTensionStorageToFreeStorageInTheUpperLayer;  /*!< flow from tension storage to free storage in the upper layer (Was tens2free_1). */
    double flowFromTensionStorageToFreeStorageInTheLowerLayer;  /*!< flow from tension storage to free storage in the lower layer (Was tens2free_2). */
    double interflowFromFreeWaterInTheUpperLayer;               /*!< interflow from free water in the upper layer (Was qintf_1). */
    double percolationFromTheUpperToLowerSoilLayers;            /*!< percolation from the upper to lower soil layers (Was qperc_12). */
    double totalBaseflow;                                       /*!< total baseflow (Was qbase_2). */
    double fractionOfBaseflowInThePrimaryBaseflowReservoir;     /*!< fraction of baseflow in the primary baseflow reservoir (Was qbase_2a). */
    double fractionOfBaseflowInTheSecondaryBaseflowReservoir;   /*!< fraction of baseflow in the secondary baseflow reservoir (Was qbase_2b). */
    double overflowFromTheUpperSoilLayer;                       /*!< overflow from the upper soil layer (Was oflow_1). */
    double overflowFromTheLowerSoilLayer;                       /*!< overflow from the lower soil layer (Was oflow_2). */
    double overflowFromTheLowerSoilLayerPrimary;                /*!< overflow from the lower soil layer (primary baseflow reservoir) (Was oflow_2a). */
    double overflowFromTheLowerSoilLayerSecondary;              /*!< overflow from the lower soil layer (secondary baseflow reservoir) (Was oflow_2b). */
};




#endif // SMASTRUCTURES_H

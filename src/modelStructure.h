#ifndef MODELSTRUCTURE_H
#define MODELSTRUCTURE_H

/*! \file modelStructure.h
\brief Tools and functionality for dealing with Model Structures.
modelStructure.h contains various functions and data types for working with Model Structures.
*/

#include "Fuse_global.h"

#include <iostream>
#include <vector>

using std::vector;
using std::cout;
using std::endl;




// Types

//! Rainfall error.
/*! This option allows to include rainfall error model parameters in the inference. In this way, the variance of rainfall errors can be estimated during hydrological model calibration rather than being fixed a priori. 

Reference: Renard, B., D. Kavetski, E. Leblois, M. Thyer, G. Kuczera, and S. W. Franks (2011), Toward a reliable decomposition of predictive uncertainty in hydrological modeling: Characterizing rainfall errors using conditional simulation, Water Resour. Res., 47, W11516, doi:10.1029/2011WR010643.

 */
enum RainfallErrorType{
    Additive        = 11,   /*!< additive. */
    Multiplicative  = 12    /*!< multiplicative. */
};




//! Architecture of the upper soil layer.
/*! The Upper Soil Layer is defined as the portion of subsurface which lays above the water table (unsaturated zone). */
enum UpperSoilLayerArchType{
    SingleState             = 21,   /*!< Single state. */
    SeparateTensionStorage  = 22,   /*!< Separate tension storage (free storage + tension storage). */
    CascadingBuckets        = 23    /*!< Cascading buckets (free storage + tension storage sub-divided into recharge and excess). */
};




//! Architecture of the lower soil layer.
/*! The Lower Soil Layer is defined as the portion of subsurface which lays below the water table (saturated zone). */
enum LowerSoilLayerArchType{
    FixedSizeBaseflowReservoir                      = 31,   /*!< Baseflow reservoir of fixed size. */
    TensionReservoirPlusTwoParallelTanks            = 32,   /*!< Tension reservoir plus two parallel tanks. */
    UnlimitedSizeBaseflowReservoirFracRate          = 33,   /*!< Baseflow reservoir of unlimited size (frac rate). */
    UnlimitedSizeBaseflowReservoirPowerRecession    = 34    /*!< Baseflow reservoir of unlimited size (power recession). */
};




//! Runoff.
/*! Parametrization to compute the saturated area and surface runoff. */
enum RunoffType{
    UnsaturatedZonePareto       = 41,   /*!< Unsaturated zone Pareto (arno/xzang/vic parameterization, upper zone control). */
    UnsaturatedZoneLinear       = 42,   /*!< Unsaturated zone linear (prms variant, fraction of upper tension storage). */
    SaturatedZoneTopographic    = 43    /*!< Saturated zone topographic (topmodel parameterization). */
};




//! Percolation.
/*! Schematization to compute percolation of water between Upper and Lower Soil Layers. */
enum PercolationType{
    DrainageAboveFieldCapacity  = 51,   /*!< Drainage above field capacity (water from field capacity to saturation available for percolation). */
    GravityDrainage             = 52,   /*!< Gravity drainage (defined by moisture content in lower layer). */
    SaturatedZoneControl        = 53    /*!< Saturated zone control (water from wilt point to saturation available for percolation). */
};




//! Evaporation.
/*! Schematization to compute evaporation from the Upper and Lower Soil Layers. */
enum EvaporationType{
    RootWeighting   = 61,   /*!< Sequential. */
    Sequential      = 62    /*!< Root weighting. */
};




//! Interflows.
/*! Lateral flux through the Upper Soil Layer to the stream. */
enum InterflowsType{
    InterflowDenied     = 71,   /*!< Interflow denied. */
    InterflowAllowed    = 72    /*!< Interflow allowed. */
};




//! Routing.
/*! Overland routing based on the Gamma Distribution with shape parameter 2.5 */
enum RoutingType{
    RoutingDenied                           = 81,   /*!< Routing denied. */
    RoutingAllowedUsingGammaDistribution    = 82    /*!< Routing allowed using Gamma distribution. */
};




//! Model Structure.
/*! The model structure contains information the model building options. */
struct ModelStructure{
    RainfallErrorType       rainfallError;      /*!< Rainfall error. */
    UpperSoilLayerArchType  upperSoilLayerArch; /*!< Architecture of the upper soil layer. */
    LowerSoilLayerArchType  lowerSoilLayerArch; /*!< Architecture of the lower soil layer. */
    RunoffType              runoff;             /*!< Runoff. */
    PercolationType         percolation;        /*!< Percolation. */
    EvaporationType         evaporation;        /*!< Evaporation. */
    InterflowsType          interflows;         /*!< Interflows. */
    RoutingType             routing;            /*!< Routing. */
};




// Functions

/*! \brief Loads Model Structures from a file.
    \param fileName name of the file containing the model configurations.
    \param modelStructures pointer to a vector created and passed by the caller which will be populated
       with the Model Structures found in fileName.
    \return returns true on success, false on failure.
    \warning Note that if the '~' character (or similar) is used, it will need to be expanded before it
         is passed to this function.
*/
bool getModelStructures(const char* fileName, vector<ModelStructure>* modelStructures);




#endif // MODELSTRUCTURE_H

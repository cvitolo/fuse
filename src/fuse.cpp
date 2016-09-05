// Use CRAN package BH for Boost headers

// [[fuse::depends(BH)]]
#include "fuse.h"
#include "modelStructure.h"
#include "smaStructures.h"

#include <boost/math/distributions/gamma.hpp>
#include <boost/numeric/odeint.hpp>

using boost::math::gamma_distribution;
using boost::math::cdf;
using boost::math::isnan;

// Global pointers for ODE solver
SmaInputParameters* globalInputParams = 0;
SmaDerivedParameters* globalDerivedParams = 0;
double* globalPrecipitation = 0;
int globalPrecipitationLength = 0;
double* globalPotentialEvapotrans = 0;
int globalPotentialEvapotransLength = 0;
ModelStructure* globalSmodl = 0;
int globalTimeCounter = 0;
double* globalSavedStates = 0;
bool globalWriteSavedStates = false;

struct streaming_observer
{
    std::ostream &m_out;
    streaming_observer( std::ostream &out ) : m_out( out ) {}

    void operator()( const state_type &x , double t ) const
    {
        // m_out << "In observer at t=" << t << endl;
        if(t > 0){
            upstates((state_type&)x);
            if(globalWriteSavedStates){
                for(int j=0; j<10; j++){
                    globalSavedStates[globalTimeCounter*10 + j] = x[j];
                }
            }
            globalTimeCounter += 1;
        }

        //m_out << t;
        //for( size_t i=0 ; i < x.size() ; ++i )
        //    m_out << "\t" << x[i];
        //m_out << "\n";
    }
};




bool fuseRoutingSim( double* instantaneousDischarge, int instantaneousDischargeSize, int* modelParams, int modelParamsSize, double timeDelay, double timeStep, double* routedDischarge, int routedDischargeSize ){

    if(modelParamsSize != 8){
        // Looks like the mode parameters are not sensible
        return false;
    }
    
    ModelStructure ms;
    
    ms.rainfallError        = (RainfallErrorType)modelParams[0]; // rferr
    ms.upperSoilLayerArch   = (UpperSoilLayerArchType)modelParams[1]; // arch1
    ms.lowerSoilLayerArch   = (LowerSoilLayerArchType)modelParams[2]; // arch2
    ms.runoff               = (RunoffType)modelParams[3]; // qsurf
    ms.percolation          = (PercolationType)modelParams[4]; // qperc
    ms.evaporation          = (EvaporationType)modelParams[5]; // esoil
    ms.interflows           = (InterflowsType)modelParams[6]; // qintf
    ms.routing              = (RoutingType)modelParams[7]; // q_tdh

    /*
    bool debug = true;
    if(debug){
        cout << endl << endl << "fuseRoutingSim called with the following parameters:" << endl;
        cout << "instantaneousDischarge at " << &instantaneousDischarge[0] << endl;
        cout << "instantaneousDischargeSize" << " : " << instantaneousDischargeSize << endl;
        cout << "modelId" << " : " << modelId << endl;
        cout << "modelListFileName at " << &(modelListFileName[0]) << endl;
        cout << "modelListFileName" << " : " << modelListFileName << endl;
        cout << "modelListFileName[0]" << " : " << modelListFileName[0] << endl;
        cout << "timeDelay" << " : " << timeDelay << endl;
        cout << "timeStep" << " : " << timeStep << endl;
        cout << "routedDischarge at " << &routedDischarge[0] << endl;
        cout << "routedDischargeSize" << " : " << routedDischargeSize << endl;
    }
    */

    // Assing standard parameters
    int sizeFracFuture = 500; // fraction of runoff in future time steps
    float alpha = 2.5;          // shape parameter

    RoutingType routingType = ms.routing;

    double* fracFuture = qTimeDelay(timeDelay,
                                    timeStep,
                                    sizeFracFuture);

    double* future  = new double[instantaneousDischargeSize];

    for(int i=0; i<instantaneousDischargeSize; i++){
        future[i] = 0.0;
        routedDischarge[i] = 0.0;
    }

    if( routingType == RoutingDenied ){
        for(int i=0; i<instantaneousDischargeSize; i++){
            routedDischarge[i] = instantaneousDischarge[i];
        }
    }
    else if( routingType == RoutingAllowedUsingGammaDistribution ){
        for(int i=0; i<instantaneousDischargeSize; i++){
            for(int j=0; j<sizeFracFuture; j++){
                future[j] = future[j] + instantaneousDischarge[i] * fracFuture[j];
            }
            routedDischarge[i] = future[0];

            // Shift everything back
            for(int j=0; j<(sizeFracFuture-1); j++){
                future[j] = future[j+1];
            }
            future[sizeFracFuture-1] = 0.0;
        }
    }

    // Clean up any temp arrays
    delete[] future;

    return true;

}




void fuseRoutingSimR( double* instantaneousDischarge, 
                      int* instantaneousDischargeSize, 
                      int* modelParams,
                      int* modelParamsSize, 
                      double* timeDelay, 
                      double* timeStep, 
                      double* routedDischarge, 
                      int* routedDischargeSize, 
                      int* status ){
    if(fuseRoutingSim(  instantaneousDischarge,
                        *instantaneousDischargeSize,
                        modelParams,
                        *modelParamsSize,
                        *timeDelay,
                        *timeStep,
                         routedDischarge,
                        *routedDischargeSize) ){
        *status = 0;
    }else{
        *status = 1;
    }

    for(int i=0; i<*routedDischargeSize; i++){
        if( ::isnan(routedDischarge[i]) ){
            routedDischarge[i] = -999.0;
        }
    }

}




bool fusesmaSim(    double* precipitation,
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
                    double *states,
                    double *fluxes,
                    double abs_error,
                    double rel_error,
                    bool correctNegativeVolumes,
                    double fracstate0,
                    double rferr_add,
                    double rferr_mlt,
                    bool outputStates,
                    bool outputFluxes,
                    bool debug){

    /*
      In the original code we had:
      stopifnot(c("P","E") %in% colnames(DATA))

      There's no concept of column names in C arrays so I have instead opted for passing
      P and E individually.
      */

    using namespace boost::numeric::odeint;

    /*if(debug){
        cout << endl << endl << "fusesmaSim called with the following parameters:" << endl;
        cout << "precipitation at " << &precipitation[0] << endl;
        cout << "precipitation[0]" << " : " << precipitation[0] << endl;
        cout << "potentialEvapotransLength" << " : " << potentialEvapotransLength << endl;
        cout << "modelId" << " : " << modelId << endl;
        cout << "modelListFileName at " << modelListFileName << endl;
        cout << "modelListFileName[0]" << " : " << modelListFileName[0] << endl;

        cout << "timeStep" << " : " << timeStep << endl;
        cout << "frchzne" << " : " << frchzne << endl;
        cout << "fracten" << " : " << fracten << endl;
        cout << "maxwatr_1" << " : " << maxwatr_1 << endl;
        cout << "percfrac" << " : " << percfrac << endl;
        cout << "fprimqb" << " : " << fprimqb << endl;
        cout << "qbrate_2a" << " : " << qbrate_2a << endl;
        cout << "qbrate_2b" << " : " << qbrate_2b << endl;
        cout << "qb_prms" << " : " << qb_prms << endl;
        cout << "maxwatr_2" << " : " << maxwatr_2 << endl;
        cout << "baserte" << " : " << baserte << endl;
        cout << "rtfrac1" << " : " << rtfrac1 << endl;
        cout << "percrte" << " : " << percrte << endl;
        cout << "percexp" << " : " << percexp << endl;
        cout << "sacpmlt" << " : " << sacpmlt << endl;
        cout << "sacpexp" << " : " << sacpexp << endl;
        cout << "iflwrte" << " : " << iflwrte << endl;
        cout << "axv_bexp" << " : " << axv_bexp << endl;
        cout << "sareamax" << " : " << sareamax << endl;
        cout << "loglamb" << " : " << loglamb << endl;
        cout << "tishape" << " : " << tishape << endl;
        cout << "qb_powr" << " : " << qb_powr << endl;

        cout << "states at " << &(states) << endl;
        cout << "states[0]" << " : " << states[0] << endl;
        cout << "fluxes at " << &(fluxes) << endl;
        cout << "fluxes[0]" << " : " << fluxes[0] << endl;
        cout << "outputStates" << " : " << outputStates << endl;
        cout << "outputStates" << " : " << outputStates << endl;
        cout << "fracstate0" << " : " << fracstate0 << endl;
        cout << "rferr_add" << " : " << rferr_add << endl;
        cout << "rferr_mlt" << " : " << rferr_mlt << endl;

    }*/

	if(modelParamsSize != 8){
        // Looks like the mode parameters are not sensible
        return false;
    }

    ModelStructure smodl;
    smodl.rainfallError        = (RainfallErrorType)modelParams[0]; // rferr
    smodl.upperSoilLayerArch   = (UpperSoilLayerArchType)modelParams[1]; // arch1
    smodl.lowerSoilLayerArch   = (LowerSoilLayerArchType)modelParams[2]; // arch2
    smodl.runoff               = (RunoffType)modelParams[3]; // qsurf
    smodl.percolation          = (PercolationType)modelParams[4]; // qperc
    smodl.evaporation          = (EvaporationType)modelParams[5]; // esoil
    smodl.interflows           = (InterflowsType)modelParams[6]; // qintf
    smodl.routing              = (RoutingType)modelParams[7]; // q_tdh

    SmaInputParameters inputParams;
    inputParams.additiveRainfallError = rferr_add;
    inputParams.multiplicativeRainfallError = rferr_mlt;
    inputParams.fracTensionStorageInRechargeZone = frchzne;
    inputParams.fracTotalStorageAsTensionStorage = fracten;
    inputParams.maximumTotalStorageInLayer1 = maxwatr_1;
    inputParams.fractionOfPercolationToTensionStorage = percfrac;
    inputParams.fractionOfBaseflowInPrimaryResvr = fprimqb;
    inputParams.baseflowDepletionRateForPrimaryResvr = qbrate_2a;
    inputParams.baseflowDepletionRateForSecondaryResvr = qbrate_2b;
    inputParams.baseflowDepletionRate = qb_prms;
    inputParams.maximumTotalStorageInLayer2 = maxwatr_2;
    inputParams.baseflowRate = baserte;
    inputParams.fractionOfRootsInTheUpperLayer = rtfrac1;
    inputParams.percolationRate = percrte;
    inputParams.percolationExponent = percexp;
    inputParams.multiplierInTheSacModelForDryLowerLayer = sacpmlt;
    inputParams.exponentInTheSacModelForDryLowerLayer = sacpexp;
    inputParams.interflowRate = iflwrte;
    inputParams.bExponent = axv_bexp;
    inputParams.maximumSaturatedArea = sareamax;
    inputParams.meanValueOfTheLogTransformedTopographicIndex = loglamb;
    inputParams.shapeParameterForTheTopoIndexGammaDistribution = tishape;
    inputParams.baseflowExponent = qb_powr;

    // compute derived parameters (bucket sizes, etc.)
    SmaDerivedParameters derivedParams = computeDerivedModelParams(smodl, inputParams);

    if(debug){
        cout << "Start of CSV model parameter debugging output..." << endl << endl;

        cout << "New Name,Old Name,Value" << endl;
        cout << "additiveRainfallError,rferr_add," << inputParams.additiveRainfallError << endl;
        cout << "multiplicativeRainfallError,rferr_mlt," << inputParams.multiplicativeRainfallError << endl;
        cout << "maximumTotalStorageInLayer1,maxwatr_1," << inputParams.maximumTotalStorageInLayer1 << endl;
        cout << "maximumTotalStorageInLayer2,maxwatr_2," << inputParams.maximumTotalStorageInLayer2 << endl;
        cout << "fracTotalStorageAsTensionStorage,fracten," << inputParams.fracTotalStorageAsTensionStorage << endl;
        cout << "fracTensionStorageInRechargeZone,frchzne," << inputParams.fracTensionStorageInRechargeZone << endl;
        cout << "fractionOfBaseflowInPrimaryResvr,fprimqb," << inputParams.fractionOfBaseflowInPrimaryResvr << endl;
        cout << "fractionOfRootsInTheUpperLayer,rtfrac1," << inputParams.fractionOfRootsInTheUpperLayer << endl;
        cout << "percolationRate,percrte," << inputParams.percolationRate << endl;
        cout << "percolationExponent,percexp," << inputParams.percolationExponent << endl;
        cout << "multiplierInTheSacModelForDryLowerLayer,sacpmlt," << inputParams.multiplierInTheSacModelForDryLowerLayer << endl;
        cout << "exponentInTheSacModelForDryLowerLayer,sacpexp," << inputParams.exponentInTheSacModelForDryLowerLayer << endl;
        cout << "fractionOfPercolationToTensionStorage,percfrac," << inputParams.fractionOfPercolationToTensionStorage << endl;
        cout << "interflowRate,iflwrte," << inputParams.interflowRate << endl;
        cout << "baseflowRate,baserte," << inputParams.baseflowRate << endl;
        cout << "baseflowExponent,qb_powr," << inputParams.baseflowExponent << endl;
        cout << "baseflowDepletionRate,qb_prms," << inputParams.baseflowDepletionRate << endl;
        cout << "baseflowDepletionRateForPrimaryResvr,qbrate_2a," << inputParams.baseflowDepletionRateForPrimaryResvr << endl;
        cout << "baseflowDepletionRateForSecondaryResvr,qbrate_2b," << inputParams.baseflowDepletionRateForSecondaryResvr << endl;
        cout << "maximumSaturatedArea,sareamax," << inputParams.maximumSaturatedArea << endl;
        cout << "bExponent,axv_bexp," << inputParams.bExponent << endl;
        cout << "meanValueOfTheLogTransformedTopographicIndex,loglamb," << inputParams.meanValueOfTheLogTransformedTopographicIndex << endl;
        cout << "shapeParameterForTheTopoIndexGammaDistribution,tishape," << inputParams.shapeParameterForTheTopoIndexGammaDistribution << endl;

        cout << "maximumTensionStorageInUpperLayer,maxtens_1," << derivedParams.maximumTensionStorageInUpperLayer << endl;
        cout << "maximumStorageInThePrimaryTensionReservoir,maxtens_1a," << derivedParams.maximumStorageInThePrimaryTensionReservoir << endl;
        cout << "maximumStorageInTheSecondaryTensionReservoir,maxtens_1b," << derivedParams.maximumStorageInTheSecondaryTensionReservoir << endl;
        cout << "maximumFreeStorageInUpperLayer,maxfree_1," << derivedParams.maximumFreeStorageInUpperLayer << endl;
        cout << "maximumTensionStorageInLowerLayer,maxtens_2," << derivedParams.maximumTensionStorageInLowerLayer << endl;
        cout << "maximumFreeStorageInLowerLayer,maxfree_2," << derivedParams.maximumFreeStorageInLowerLayer << endl;
        cout << "maximumStorageInThePrimaryBaseflowReservoir,maxfree_2a," << derivedParams.maximumStorageInThePrimaryBaseflowReservoir << endl;
        cout << "maximumStorageInTheSecondaryBaseflowReservoir,maxfree_2b," << derivedParams.maximumStorageInTheSecondaryBaseflowReservoir << endl;
        cout << "fractionOfRootsInTheLowerSoilLayer,rtfrac2," << derivedParams.fractionOfRootsInTheLowerSoilLayer << endl;
        cout << "baseflowAtSaturation,qbsat," << derivedParams.baseflowAtSaturation << endl;
        cout << "meanOfThePowerTransformedTopographicIndex,powlamb," << derivedParams.meanOfThePowerTransformedTopographicIndex << endl;
        cout << "powerTransformedTopographicIndex,maxpow," << derivedParams.powerTransformedTopographicIndex << endl;

        cout << "End of CSV model parameter debugging output..." << endl;
    }

    // initialize model states (upper layer and lower layer)
    SmaStateStructure state0 = initStates(smodl, inputParams, derivedParams, fracstate0);

    // Solve derivatives
    /*double* times = new double[precipitationLength];
    for(int i=1; i<=precipitationLength; i++){
        times[i] = double(i);
    }*/

    // Start Solver section

    // Set up pointers to parameters:
    globalInputParams = &inputParams;
    globalDerivedParams = &derivedParams;
    globalPrecipitation = precipitation;
    globalPrecipitationLength = precipitationLength;
    globalPotentialEvapotrans = potentialEvapotrans;
    globalPotentialEvapotransLength = potentialEvapotransLength;
    globalSmodl = &smodl;

    state_type x( 10 );

    // Transfer state0 into x
    x[0] = state0.excessTensionStorageUpperLayer;
    x[1] = state0.rechargeTensionStorageUpperLayer;
    x[2] = state0.tensionStorageUpperLayer;
    x[3] = state0.freeStorageUpperLayer;
    x[4] = state0.totalUpperLayerStorage;
    x[5] = state0.tensionStorageLowerLayer;
    x[6] = state0.freeStoragePrimaryBaseflowReservoir;
    x[7] = state0.freeStorageSecondaryBaseflowReservoir;
    x[8] = state0.totalLowerLayerStorage;
    x[9] = state0.freeStorageBaseflowReservoir;

    typedef controlled_runge_kutta< runge_kutta_cash_karp54< state_type > > stepper_type;

    //typedef controlled_runge_kutta< runge_kutta_dopri5< state_type > > dopri_stepper_type;
    //typedef dense_output_runge_kutta< dopri_stepper_type > dense_stepper_type;

    //const double dt = timeStep;
    /*integrate_const( dense_stepper_type(),
                        mDifEqn,
                        x,
                        0.0,
                        10.0, //(precipitationLength * timeStep),
                        1.0, // timeStep,
                        streaming_observer(std::cout) );*/

    globalTimeCounter = 0;
    globalWriteSavedStates = outputStates;
    globalSavedStates = states;
    // FIXME - our code assumes that the first time will be t=0 - check this
    integrate_const(stepper_type( default_error_checker< double >( abs_error , rel_error ) ),
                    mDifEqn,
                    x,
                    0.0,
                    (precipitationLength * timeStep), // e.g. 10 array elements, 1 second time step, 10 seconds run time, outputs at 0, 1, ..., 10 (gives
                    timeStep,
                    streaming_observer(std::cout) );

    /*
    return true;

    // Note:  timeStep is the time step
    runge_kutta4< state_type > rk4;

    double t = 0.0;
    for( size_t i=0; i<precipitationLength; ++i, t += timeStep ){
        globalTimeCounter = i;
        // at this point, we call x "state a"
        rk4.do_step( mDifEqn, x, t, timeStep );
        // at this point, x is "state b" ( state a + (dxdt * timeStep) )

        // Call upstates
        upstates(x);

        if(states){
            // Write the states to an array for later
            for(int j=0; j<10; j++){
                s[j][i] = x[j];
            }
        }
    }

    // End solver section

    */

    // Compute fluxes

    outFluxes(&smodl,
              precipitation,
              precipitationLength,
              potentialEvapotrans,
              potentialEvapotransLength,
              &inputParams,
              &derivedParams,
              states,

              fluxes,
              correctNegativeVolumes);

    return true;

}




void fusesmaSimR(   double* precipitation,
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
                    int* status){

    if(fusesmaSim(   precipitation,
                    *precipitationLength,
                     potentialEvapotrans,
                    *potentialEvapotransLength,
                    modelParams,
                    *modelParamsSize,
                    *timeStep,
                    *frchzne,
                    *fracten,
                    *maxwatr_1,
                    *percfrac,
                    *fprimqb,
                    *qbrate_2a,
                    *qbrate_2b,
                    *qb_prms,
                    *maxwatr_2,
                    *baserte,
                    *rtfrac1,
                    *percrte,
                    *percexp,
                    *sacpmlt,
                    *sacpexp,
                    *iflwrte,
                    *axv_bexp,
                    *sareamax,
                    *loglamb,
                    *tishape,
                    *qb_powr,
                     states,
                     fluxes,
                    *abs_error,
                    *rel_error,
                    *correctNegativeVolumes,
                    *fracstate0,
                    *rferr_add,
                    *rferr_mlt) ){
        *status = 0;
    }else{
        *status = 1;
    }

    /*
            for(int i=0; i<precipitationValues.size(); i++){
                for(int j=0; j<10; j++){
                    SOutFile << outS[i*10 + j] << ",";
                }
                SOutFile << endl;
                for(int j=0; j<20; j++){
                    FOutFile << outF[i*20 + j] << ",";
                }
                FOutFile << endl;
                UOutFile << outU[i] << endl;
            }
      */

    for(int i=0; i<*precipitationLength; i++){
        for(int j=0; j<20; j++){
            if( ::isnan(fluxes[i*20 + j]) )
                fluxes[i*20 + j] = -999.0;
        }
        for(int j=0; j<10; j++){
            if( ::isnan(states[i*10 + j]) )
                states[i*10 + j] = -999.0;
        }
    }

}




double baseflowAtSaturation(double lowerSoilLayerArch,
                            SmaInputParameters inputParams,
                            double maximumStorageInThePrimaryBaseflowReservoir,
                            double maximumStorageInTheSecondaryBaseflowReservoir,
                            double meanOfThePowerTransformedTopographicIndex){

    /*
        Computes baseflow at saturation (used in the sac percolation model)
        Author: Claudia Vitolo

        Args:
          arch_2:
          inputParams:                        list of parameters obtained from assign_par
          maxfree_2a:
          maxfree_2b:
          powlamb:

        Returns:
          Baseflow at saturation.
    */

    double qbsat = 0.0;

    if(lowerSoilLayerArch == TensionReservoirPlusTwoParallelTanks){
        // tension reservoir plus two parallel tanks
        qbsat = inputParams.baseflowDepletionRateForPrimaryResvr * maximumStorageInThePrimaryBaseflowReservoir
                + inputParams.baseflowDepletionRateForSecondaryResvr * maximumStorageInTheSecondaryBaseflowReservoir;
    }

    if(lowerSoilLayerArch == UnlimitedSizeBaseflowReservoirFracRate){
        // baseflow resvr of unlimited size
        qbsat = inputParams.baseflowDepletionRate * inputParams.maximumTotalStorageInLayer2;
    }

    if(lowerSoilLayerArch == UnlimitedSizeBaseflowReservoirPowerRecession){
        // topmodel power-law transmissivity profile
        // this is a bit tricky.  the capacity of the aquifer is m*n, where m is a scaling parameter.
        // we have the capacity, i.e., maxwatr_2/1000., and need the topmodel "m" parameter
        double topmdm = (inputParams.maximumTotalStorageInLayer2 / 1000.0) / inputParams.baseflowExponent;
        // compute baseflow
        qbsat = inputParams.baseflowRate * ( topmdm / pow(meanOfThePowerTransformedTopographicIndex, inputParams.baseflowExponent) );
    }

    /* DELETED (depricated)
    if(lowerSoilLayerArch == 35) {   # topmodel exponential transmissivity profile
        topmdm <- inputParams$maxwatr_2 / 1000                          # for simplicity we use the capacity as the topmodel scaling parameter
        qbsat  <- inputParams$baserte * topmdm * exp(-inputParams$loglamb)   # compute baseflow
    }
    */

    if(lowerSoilLayerArch == FixedSizeBaseflowReservoir){
        // baseflow reservoir of fixed size
        qbsat = inputParams.baseflowRate;
    }

    return qbsat;

}




SmaDerivedParameters computeDerivedModelParams(ModelStructure smodl, SmaInputParameters inputParams){
    /*
        Computes derived model parameters (bucket sizes, etc.)
        Author: Claudia Vitolo

        Args:
          smodl:                         list of model components
          inputParams:                   list of parameters obtained from assign_par

        Returns:
          List of derived parameters.
    */

    SmaDerivedParameters dp;

    // BUCKETSIZE

    // derive maximum tension water in upper layer
    dp.maximumTensionStorageInUpperLayer = inputParams.fracTotalStorageAsTensionStorage * inputParams.maximumTotalStorageInLayer1;
    // derive maximum free water in upper layer
    dp.maximumFreeStorageInUpperLayer = (1 - inputParams.fracTotalStorageAsTensionStorage) * inputParams.maximumTotalStorageInLayer1;

    // derive maximum tension water in lower layer
    dp.maximumTensionStorageInLowerLayer = inputParams.fracTotalStorageAsTensionStorage * inputParams.maximumTotalStorageInLayer2;
    // derive maximum free water in lower layer
    dp.maximumFreeStorageInLowerLayer = (1 - inputParams.fracTotalStorageAsTensionStorage) * inputParams.maximumTotalStorageInLayer2;

    // derive capacities of the recharge and lower zone (only used if upper tension is divided in two)
    dp.maximumStorageInThePrimaryTensionReservoir = inputParams.fracTensionStorageInRechargeZone * dp.maximumTensionStorageInUpperLayer;
    dp.maximumStorageInTheSecondaryTensionReservoir = (1 - inputParams.fracTensionStorageInRechargeZone) * dp.maximumTensionStorageInUpperLayer;

    // derive capacities of the primary and secondary parallel baseflow reservoirs
    dp.maximumStorageInThePrimaryBaseflowReservoir = inputParams.fractionOfBaseflowInPrimaryResvr * dp.maximumFreeStorageInLowerLayer;
    dp.maximumStorageInTheSecondaryBaseflowReservoir = (1 - inputParams.fractionOfBaseflowInPrimaryResvr) * dp.maximumFreeStorageInLowerLayer;

    // fraction of roots in the lower layer (-)
    if ( smodl.evaporation == RootWeighting ) {
        dp.fractionOfRootsInTheLowerSoilLayer = 1.0 - inputParams.fractionOfRootsInTheUpperLayer;
    } else {
        dp.fractionOfRootsInTheLowerSoilLayer = 0.0;
    }

    // mean of the power-transformed topographic index
    meanTipow(inputParams.shapeParameterForTheTopoIndexGammaDistribution,
              inputParams.meanValueOfTheLogTransformedTopographicIndex,
              inputParams.baseflowExponent,

              &dp.powerTransformedTopographicIndex, //in fortran code maxpow/10, why???
              &dp.meanOfThePowerTransformedTopographicIndex);

    // compute baseflow at saturation (used in the sac percolation model)
    dp.baseflowAtSaturation = baseflowAtSaturation(smodl.lowerSoilLayerArch,
                                                   inputParams,
                                                   dp.maximumStorageInThePrimaryBaseflowReservoir,
                                                   dp.maximumStorageInTheSecondaryBaseflowReservoir,
                                                   dp.meanOfThePowerTransformedTopographicIndex);

    return dp;

}




SmaFluxStructure computeFluxes(ModelStructure* smodl,
                               double* precipitationAtT,
                               double* potentialEvapotransAtT,
                               SmaInputParameters* inputParams,
                               SmaDerivedParameters* derivedParams,
                               state_type* x){

    /*
        Compute Fluxes
        Author: Claudia Vitolo

        Args:
          smodl:                         list of model components
          mppt:                          rain+snow melt at time "t"
          mpet:                          potential evapotranspiration at time "t"
          inputParams:                   model parameters
          dparam:                        derived model parameters
          state:                         model states at time "t"

        Returns:
          List of fluxes at time "t"
    */

    SmaFluxStructure fs;

    // set all fluxes to zero at the start of time step
    fs.effectiveRainfall = 0.0;
    fs.saturatedArea = 0.0;
    fs.surfaceRunoff = 0.0;
    fs.evaporationFromTheUpperLayerRecharge = 0.0;
    fs.evaporationFromTheUpperLayerExcess = 0.0;
    fs.evaporationFromTheUpperLayer = 0.0;
    fs.evaporationFromTheLowerLayer = 0.0;
    fs.flowFromRechargeToExcess = 0.0;
    fs.flowFromTensionStorageToFreeStorageInTheUpperLayer = 0.0;
    fs.flowFromTensionStorageToFreeStorageInTheLowerLayer = 0.0;
    fs.interflowFromFreeWaterInTheUpperLayer = 0.0;
    fs.percolationFromTheUpperToLowerSoilLayers = 0.0;
    fs.totalBaseflow = 0.0;
    fs.fractionOfBaseflowInThePrimaryBaseflowReservoir = 0.0;
    fs.fractionOfBaseflowInTheSecondaryBaseflowReservoir = 0.0;
    fs.overflowFromTheUpperSoilLayer = 0.0;
    fs.overflowFromTheLowerSoilLayer = 0.0;
    fs.overflowFromTheLowerSoilLayerPrimary = 0.0;
    fs.overflowFromTheLowerSoilLayerSecondary = 0.0;

    // compute effective rainfall
    if(smodl->rainfallError == Additive)
        fs.effectiveRainfall = std::max(0.0, *precipitationAtT + inputParams->additiveRainfallError);
    if(smodl->rainfallError == Multiplicative)
        fs.effectiveRainfall = *precipitationAtT * inputParams->multiplicativeRainfallError;

    // compute excess of saturation

    fluxQSatExcess(smodl->runoff,
                   fs.effectiveRainfall,
                   inputParams,
                   derivedParams,
                   x,
                   &fs.saturatedArea,
                   &fs.surfaceRunoff );

    // compute evaporation
    fluxEvap(smodl->upperSoilLayerArch,
             smodl->lowerSoilLayerArch,
             smodl->evaporation,
             potentialEvapotransAtT,
             inputParams,
             derivedParams,
             x,
             &fs.evaporationFromTheUpperLayerRecharge,
             &fs.evaporationFromTheUpperLayerExcess,
             &fs.evaporationFromTheUpperLayer,
             &fs.evaporationFromTheLowerLayer);

    // compute interflow from free water in the upper layer
    if(smodl->interflows == InterflowAllowed)
        fs.interflowFromFreeWaterInTheUpperLayer = inputParams->interflowRate * ( x->at(3) / derivedParams->maximumFreeStorageInUpperLayer);
    if(smodl->interflows == InterflowDenied)
        fs.interflowFromFreeWaterInTheUpperLayer = 0.0;

    // compute percolation from the upper to lower soil layers
    fs.percolationFromTheUpperToLowerSoilLayers = fluxperc(smodl->percolation,
                                                           inputParams,
                                                           derivedParams,
                                                           x);

    // compute baseflow
    fluxBaseflow(smodl->lowerSoilLayerArch,
                 inputParams,
                 derivedParams->baseflowAtSaturation,
                 x,
                 &fs.fractionOfBaseflowInThePrimaryBaseflowReservoir,
                 &fs.fractionOfBaseflowInTheSecondaryBaseflowReservoir,
                 &fs.totalBaseflow);

    // compute overflow (miscellaneous fluxes)
    fluxQMisc(smodl->upperSoilLayerArch,
              smodl->lowerSoilLayerArch,
              x,
              inputParams,
              derivedParams,
              fs.effectiveRainfall,
              fs.surfaceRunoff,
              fs.percolationFromTheUpperToLowerSoilLayers,
              &fs.flowFromRechargeToExcess,
              &fs.flowFromTensionStorageToFreeStorageInTheUpperLayer,
              &fs.flowFromTensionStorageToFreeStorageInTheLowerLayer,
              &fs.overflowFromTheUpperSoilLayer,
              &fs.overflowFromTheLowerSoilLayer,
              &fs.overflowFromTheLowerSoilLayerPrimary,
              &fs.overflowFromTheLowerSoilLayerSecondary);

    return fs;

}




void fluxBaseflow(LowerSoilLayerArchType lowerSoilLayerArch,
                  SmaInputParameters* inputParams,
                  double baseflowAtSaturation,
                  state_type* x,
                  double* fractionOfBaseflowInThePrimaryBaseflowReservoir,
                  double* fractionOfBaseflowInTheSecondaryBaseflowReservoir,
                  double* totalBaseflow){
    /*
        Compute baseflow from the lower soil layer
        Author: Claudia Vitolo

        Args:
          arch2:         smodl$arch2
          inputParams:   useful model parameters
          qbsat:         dparam$qbsat
          free_2a:       Free Storage Primary Baseflow Reservoir
          free_2b:       Free Storage Secondary Baseflow Reservoir
          watr_2:        Total Storage Lower Layer

        Returns:
          Baseflow
    */

    *fractionOfBaseflowInThePrimaryBaseflowReservoir = 0;
    *fractionOfBaseflowInTheSecondaryBaseflowReservoir = 0;
    *totalBaseflow = 0;

    // baseflow reservoir of fixed size
    if(lowerSoilLayerArch == FixedSizeBaseflowReservoir)
        *totalBaseflow = inputParams->baseflowRate
                         * pow( (x->at(8) / inputParams->maximumTotalStorageInLayer2), inputParams->baseflowExponent );

    // tension reservoir plus two parallel tanks
    if(lowerSoilLayerArch == TensionReservoirPlusTwoParallelTanks){
        // qbrate_2a is a fraction (t-1)
        *fractionOfBaseflowInThePrimaryBaseflowReservoir =
                inputParams->baseflowDepletionRateForPrimaryResvr * x->at(6);
        // qbrate_2b is a fraction (t-1)
        *fractionOfBaseflowInTheSecondaryBaseflowReservoir = inputParams->baseflowDepletionRateForSecondaryResvr * x->at(7);
        // total baseflow
        *totalBaseflow  = *fractionOfBaseflowInThePrimaryBaseflowReservoir
                          + *fractionOfBaseflowInTheSecondaryBaseflowReservoir;
    }

    // baseflow resvr of unlimited size (0-huge), frac rate # qb_prms is a fraction (t-1)
    if(lowerSoilLayerArch == UnlimitedSizeBaseflowReservoirFracRate)
        *totalBaseflow = inputParams->baseflowDepletionRate * x->at(8);

    // baseflow resvr of unlimited size (0-huge), power recession
    if(lowerSoilLayerArch == UnlimitedSizeBaseflowReservoirPowerRecession)
        *totalBaseflow = baseflowAtSaturation * pow( (x->at(8) / inputParams->maximumTotalStorageInLayer2), inputParams->baseflowExponent);

}




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
              double* evaporationFromTheLowerLayer){

    /*
        Compute evaporation from the upper and lower layers
        Author: Claudia Vitolo

        Args:
          arch1:                         smodl$arch1
          arch2:                         smodl$arch2
          esoil:                         smodl$esoil
          fpet:                          mpet
          inputParams:                   model parameters
          dparam:                        derived model parameters
          tens_1a:                       tens_1a
          tens_1b:                       tens_1b
          tens_1:                        tens_1
          tens_2:                        tens_2

        Returns:
          Evaporation
    */

    // compute evaporation from the upper layer

    if(upperArch == CascadingBuckets){
        // tension storage sub-divided into recharge and excess
        if(evaporationType == Sequential){
            *evaporationFromTheUpperLayerRecharge = *potentialEvapotransAtT * x->at(0) / derivedParams->maximumStorageInThePrimaryTensionReservoir;
            *evaporationFromTheUpperLayerExcess = (*potentialEvapotransAtT - *evaporationFromTheUpperLayerRecharge) * x->at(1) / derivedParams->maximumStorageInTheSecondaryTensionReservoir;
            *evaporationFromTheUpperLayer = *evaporationFromTheUpperLayerRecharge + *evaporationFromTheUpperLayerExcess;
        }
        if(evaporationType == RootWeighting){
            *evaporationFromTheUpperLayerRecharge = *potentialEvapotransAtT * inputParams->fractionOfRootsInTheUpperLayer * x->at(0) / derivedParams->maximumStorageInThePrimaryTensionReservoir;
            *evaporationFromTheUpperLayerExcess = *potentialEvapotransAtT * derivedParams->fractionOfRootsInTheLowerSoilLayer * x->at(1) / derivedParams->maximumStorageInTheSecondaryTensionReservoir;
            *evaporationFromTheUpperLayer = *evaporationFromTheUpperLayerRecharge + *evaporationFromTheUpperLayerExcess;
        }
    }

    if(upperArch == SeparateTensionStorage || upperArch == SingleState){
        // single tension store or single state
        if(evaporationType == Sequential){
            *evaporationFromTheUpperLayer = *potentialEvapotransAtT * x->at(2) / derivedParams->maximumTensionStorageInUpperLayer;
        }
        if(evaporationType == RootWeighting){
            *evaporationFromTheUpperLayer = *potentialEvapotransAtT * inputParams->fractionOfRootsInTheUpperLayer * x->at(2) / derivedParams->maximumTensionStorageInUpperLayer;
        }
    }

    // compute evaporation from the lower layer
    if(lowerSoilLayerArch == TensionReservoirPlusTwoParallelTanks || lowerSoilLayerArch == FixedSizeBaseflowReservoir){
        // lower layer architecture
        if(upperArch == SeparateTensionStorage || upperArch == SingleState){
            // lower-layer evap is valid
            // use different evaporation schemes for the lower layer
            if(evaporationType == Sequential)
                *evaporationFromTheLowerLayer = (*potentialEvapotransAtT - *evaporationFromTheUpperLayer)
                                                * ( x->at(5) / derivedParams->maximumTensionStorageInLowerLayer );
            if(evaporationType == RootWeighting)
                *evaporationFromTheLowerLayer = *potentialEvapotransAtT * derivedParams->fractionOfRootsInTheLowerSoilLayer
                                                * (x->at(5) / derivedParams->maximumTensionStorageInLowerLayer);
        }
        if(upperArch == CascadingBuckets)
            // lower-layer evap is zero
            *evaporationFromTheLowerLayer = 0.0;
    }

    if(lowerSoilLayerArch == UnlimitedSizeBaseflowReservoirFracRate || lowerSoilLayerArch == UnlimitedSizeBaseflowReservoirPowerRecession)
        *evaporationFromTheLowerLayer = 0.0;

}




double fluxperc( PercolationType percolationType,
                 SmaInputParameters* inputParams,
                 SmaDerivedParameters* derivedParams,
                 state_type* x){

    /*
        Compute percolation
        Author: Claudia Vitolo

        Args:
          qperc:                         smodl$qperc
          inputParams:                   model parameters
          dparam:                        derived model parameters
          free_1:
          watr_1:
          watr_2:

        Returns:
          Percolation
    */

    double percolation = 0.0;

    if(percolationType == DrainageAboveFieldCapacity)
        percolation = inputParams->percolationRate
                      * pow( (x->at(3) / derivedParams->maximumFreeStorageInUpperLayer), inputParams->percolationExponent);

    if(percolationType == GravityDrainage){

        double lz_pd = 0.0;

        if( (x->at(8) / inputParams->maximumTotalStorageInLayer2) > 1.0 ){
            lz_pd = 1.0;
        } else {
            lz_pd = 1.0 + inputParams->multiplierInTheSacModelForDryLowerLayer
                    * pow( (1.0 - x->at(8) / inputParams->maximumTotalStorageInLayer2), inputParams->exponentInTheSacModelForDryLowerLayer);
        }

        percolation = derivedParams->baseflowAtSaturation * lz_pd
                      * (x->at(3) / derivedParams->maximumFreeStorageInUpperLayer);
    }

    if(percolationType == SaturatedZoneControl)
        percolation = inputParams->percolationRate
                      * pow( (x->at(4) / inputParams->maximumTotalStorageInLayer1), inputParams->percolationExponent);

    return percolation;
}




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
               int optionX){
    /*
        Computes miscellaneous fluxes using 2 methods:
        (x=1) OVERFLOW FLUXES AS A FRACTION OF INFLUXES
        (x=2) OVERFLOW FLUXES COMPUTED AS A RESIDUAL OF AVAILABLE STORAGE
        Author: Claudia Vitolo

        Args:
          arch1:                         architecture of the upper soil layer
          arch2:                         architecture of the lower soil layer
          state:                         model states at time "t"
          inputParams:                   model parameters
          dparam:                        derived model parameters
          eff_ppt:                       effective precipitation
          qrunoff:                       runoff
          qperc_12:                      percolation

        Returns:
          Overflow fluxes:
          rchr2excs   <- flow from recharge to excess (mm day-1)
          tens2free_1 <- flow from tension storage to free storage in the upper layer (mm day-1)
          tens2free_2 <- flow from tension storage to free storage in the lower layer (mm day-1)
          oflow_1     <- overflow from the upper soil layer (mm day-1)
          oflow_2     <- overflow from the lower soil layer (mm day-1)
          oflow_2a    <- overflow from the lower soil layer (mm day-1)
          oflow_2b    <- overflow from the lower soil layer (mm day-1)
    */

    *flowFromRechargeToExcess = 0.0;
    *flowFromTensionStorageToFreeStorageInTheUpperLayer = 0.0;
    *flowFromTensionStorageToFreeStorageInTheLowerLayer = 0.0;
    *overflowFromTheUpperSoilLayer = 0.0;
    *overflowFromTheLowerSoilLayer = 0.0;
    *overflowFromTheLowerSoilLayerPrimary = 0.0;
    *overflowFromTheLowerSoilLayerSecondary = 0.0;

    if(optionX == 1){
        if(upperArch == SingleState){
            // upper layer defined by a single state variable
            *flowFromRechargeToExcess = 0.0;
            // no tension stores
            *flowFromTensionStorageToFreeStorageInTheUpperLayer = 0.0;
            double w_func = logismooth( x->at(4), inputParams->maximumTotalStorageInLayer1 );
            // compute over-flow of free water
            *overflowFromTheUpperSoilLayer = w_func
                    * (effectiveRainfall - surfaceRunoff);
        }

        if(upperArch == SeparateTensionStorage){
            // upper layer broken up into tension and free storage
            // no separate recharge zone (flux should never be used)
            *flowFromRechargeToExcess = 0.0;
            double w_func = logismooth( x->at(2), derivedParams->maximumTensionStorageInUpperLayer );
            // compute flow from tension storage to free storage (mm s-1)
            *flowFromTensionStorageToFreeStorageInTheUpperLayer = w_func * (effectiveRainfall - surfaceRunoff);
            w_func = logismooth( x->at(3), derivedParams->maximumFreeStorageInUpperLayer );
            // compute over-flow of free water
            *overflowFromTheUpperSoilLayer = w_func *
                                             *flowFromTensionStorageToFreeStorageInTheUpperLayer;
        }

        if(upperArch == CascadingBuckets){
            // tension storage sub-divided into recharge and excess
            double w_func = logismooth( x->at(0), derivedParams->maximumStorageInThePrimaryTensionReservoir );
            // compute flow from recharge to excess (mm s-1)
            *flowFromRechargeToExcess = w_func * (effectiveRainfall - surfaceRunoff);
            w_func = logismooth( x->at(1), derivedParams->maximumStorageInTheSecondaryTensionReservoir );
            // compute flow from tension storage to free storage (mm s-1)
            *flowFromTensionStorageToFreeStorageInTheUpperLayer = w_func * *flowFromRechargeToExcess;
            w_func = logismooth( x->at(3), derivedParams->maximumFreeStorageInUpperLayer );
            // compute over-flow of free water
            *overflowFromTheUpperSoilLayer = w_func * *flowFromTensionStorageToFreeStorageInTheUpperLayer;
        }

        if(lowerSoilLayerArch == FixedSizeBaseflowReservoir){
            // no tension store
            *flowFromTensionStorageToFreeStorageInTheLowerLayer = 0.0;
            *overflowFromTheLowerSoilLayerPrimary = 0.0;
            *overflowFromTheLowerSoilLayerSecondary = 0.0;
            double w_func = logismooth( x->at(8), inputParams->maximumTotalStorageInLayer2 );
            // compute over-flow of free water
            *overflowFromTheLowerSoilLayer = w_func * percolation;
        }

        if(lowerSoilLayerArch == TensionReservoirPlusTwoParallelTanks){
            // tension reservoir plus two parallel tanks
            double w_func = logismooth( x->at(5), derivedParams->maximumTensionStorageInLowerLayer );
            // compute flow from tension storage to free storage (mm s-1)
            *flowFromTensionStorageToFreeStorageInTheLowerLayer = w_func
                                                                  * percolation
                                                                  * ( 1.0 - inputParams->fractionOfPercolationToTensionStorage);
            w_func = logismooth( x->at(6), derivedParams->maximumStorageInThePrimaryBaseflowReservoir );
            // compute over-flow of free water in the primary reservoir
            *overflowFromTheLowerSoilLayerPrimary = w_func * (percolation * (inputParams->fractionOfPercolationToTensionStorage / 2.0)
                                                             + *flowFromTensionStorageToFreeStorageInTheLowerLayer / 2.0);
            w_func = logismooth( x->at(7), derivedParams->maximumStorageInTheSecondaryBaseflowReservoir);
            // compute over-flow of free water in the secondary reservoir
            *overflowFromTheLowerSoilLayerSecondary = w_func *
                                                    ( percolation *
                                                      (inputParams->fractionOfPercolationToTensionStorage / 2.0)
                                                      + *flowFromTensionStorageToFreeStorageInTheLowerLayer / 2.0
                                                    );
            // compute total overflow
            *overflowFromTheLowerSoilLayer = *overflowFromTheLowerSoilLayerPrimary
                                             + *overflowFromTheLowerSoilLayerSecondary;
        }

        if(lowerSoilLayerArch == UnlimitedSizeBaseflowReservoirFracRate || lowerSoilLayerArch == UnlimitedSizeBaseflowReservoirPowerRecession){
            // unlimited size
            *flowFromTensionStorageToFreeStorageInTheLowerLayer = 0.0;
            *overflowFromTheLowerSoilLayer = 0.0;
            *overflowFromTheLowerSoilLayerPrimary = 0.0;
            *overflowFromTheLowerSoilLayerSecondary = 0.0;
        }

    }else if(optionX == 2){
        if(upperArch == SingleState){
            // upper layer defined by a single state variable
            *flowFromRechargeToExcess = 0.0;
            // no tension stores
            *flowFromTensionStorageToFreeStorageInTheUpperLayer = 0.0;
            // compute over-flow of free water
            *overflowFromTheUpperSoilLayer = max(0.0, (effectiveRainfall - surfaceRunoff) - (inputParams->maximumTotalStorageInLayer1 - x->at(4) ));
        }

        if(upperArch == SeparateTensionStorage){
            // upper layer broken up into tension and free storage
            // no separate recharge zone (flux should never be used)
            *flowFromRechargeToExcess = 0.0;
            // compute flow from tension storage to free storage (mm s-1)
            *flowFromTensionStorageToFreeStorageInTheUpperLayer = max(0.0, (effectiveRainfall - surfaceRunoff) - (derivedParams->maximumTensionStorageInUpperLayer - x->at(2) ));
            // compute over-flow of free water
            *overflowFromTheUpperSoilLayer = max(0.0, *flowFromTensionStorageToFreeStorageInTheUpperLayer - (derivedParams->maximumFreeStorageInUpperLayer - x->at(3) ));
        }

        if(upperArch == CascadingBuckets){
            // tension storage sub-divided into recharge and excess
            // compute flow from recharge to excess (mm s-1)
            *flowFromRechargeToExcess = max(0.0, (effectiveRainfall - surfaceRunoff) - (derivedParams->maximumStorageInThePrimaryTensionReservoir - x->at(0) ));
            // compute flow from tension storage to free storage (mm s-1)
            *flowFromTensionStorageToFreeStorageInTheUpperLayer = max(0.0, *flowFromRechargeToExcess - (derivedParams->maximumStorageInTheSecondaryTensionReservoir - x->at(1) ));
            // compute over-flow of free water
            *overflowFromTheUpperSoilLayer = max(0.0, *flowFromTensionStorageToFreeStorageInTheUpperLayer - (derivedParams->maximumFreeStorageInUpperLayer - x->at(3) ));
        }

        if(lowerSoilLayerArch == FixedSizeBaseflowReservoir){
            // no tension store
            *flowFromTensionStorageToFreeStorageInTheLowerLayer = 0.0;
            *overflowFromTheLowerSoilLayerPrimary = 0.0;
            *overflowFromTheLowerSoilLayerSecondary = 0.0;
            // compute over-flow of free water
            *overflowFromTheLowerSoilLayer = max(0.0, percolation - (inputParams->maximumTotalStorageInLayer2 - x->at(8) ));
        }

        if(lowerSoilLayerArch == TensionReservoirPlusTwoParallelTanks){
            // tension reservoir plus two parallel tanks
            // compute flow from tension storage to free storage (mm s-1)
            *flowFromTensionStorageToFreeStorageInTheLowerLayer = max(0.0, percolation * (1 - inputParams->fractionOfPercolationToTensionStorage) - (derivedParams->maximumTensionStorageInLowerLayer  - x->at(5) ));
            // compute over-flow of free water in the primary reservoir
            *overflowFromTheLowerSoilLayerPrimary = max(0.0, (percolation * (inputParams->fractionOfPercolationToTensionStorage / 2.0) + *flowFromTensionStorageToFreeStorageInTheLowerLayer / 2.0) - (derivedParams->maximumStorageInThePrimaryBaseflowReservoir - x->at(6)  ));
            // compute over-flow of free water in the secondary reservoir
            *overflowFromTheLowerSoilLayerSecondary = max(0.0, (percolation * (inputParams->fractionOfPercolationToTensionStorage / 2.0) + *flowFromTensionStorageToFreeStorageInTheLowerLayer / 2.0) - (derivedParams->maximumStorageInTheSecondaryBaseflowReservoir - x->at(7) ));
            // compute total overflow
            *overflowFromTheLowerSoilLayer = *overflowFromTheLowerSoilLayerPrimary + *overflowFromTheLowerSoilLayerSecondary;
        }

        if(lowerSoilLayerArch == UnlimitedSizeBaseflowReservoirFracRate || lowerSoilLayerArch == UnlimitedSizeBaseflowReservoirPowerRecession){
            // unlimited size
            *flowFromTensionStorageToFreeStorageInTheLowerLayer = 0.0;
            *overflowFromTheLowerSoilLayer = 0.0;
            *overflowFromTheLowerSoilLayerPrimary = 0.0;
            *overflowFromTheLowerSoilLayerSecondary = 0.0;
        }
    }

}




void fluxQSatExcess(RunoffType runoffType,
                    double effectiveRainfall,
                    SmaInputParameters* inputParams,
                    SmaDerivedParameters* derivedParams,
                    state_type* x,
                    double* saturatedArea,
                    double* surfaceRunoff){

    /*
        Computes the saturated area and surface runoff
        Author: Claudia Vitolo

        Args:
          mqsurf:      smodl$qsurf
          eff_ppt:     eff_ppt
          inputParams: inputParams
          dparam:      dparam
          tens_1:      tens_1
          watr_1:      watr_1
          watr_2:      watr_2

        Returns:
          Saturated Area and Surface Runoff
    */

    *saturatedArea = 0.0;
    *surfaceRunoff = 0.0;
    // avoid divide by zero
    double no_zero = 0.00000001;

    // saturated area method
    if(runoffType == UnsaturatedZonePareto)
        // arno/xzang/vic parameterization (upper zone control)
        *saturatedArea = 1.0 - pow( ( 1.0 - std::min( (x->at(4) / inputParams->maximumTotalStorageInLayer1), 1.0) ), inputParams->bExponent);

    if(runoffType == UnsaturatedZoneLinear)
        // prms variant (fraction of upper tension storage)
        *saturatedArea = std::min( (x->at(2) / derivedParams->maximumTensionStorageInUpperLayer), 1.0) * inputParams->maximumSaturatedArea;

    if(runoffType == SaturatedZoneTopographic){
        // topmodel parameterization (only valid for topmodel qb)
        // compute the minimum value of the topographic index where the basin is saturated
        // (this is correct, as maxwatr_2 is m*n -- units are meters**(1/n)
        double ti_sat = derivedParams->meanOfThePowerTransformedTopographicIndex / (x->at(8) / inputParams->maximumTotalStorageInLayer2 + no_zero);
        // compute the saturated area
        if(ti_sat > derivedParams->powerTransformedTopographicIndex){
            *saturatedArea = 0;
        }else{
            // convert the topographic index to log space, compute the saturated area (note: critical value of the topographic index is in log space)
            double ti_log = log( pow(ti_sat, inputParams->baseflowExponent) );
            // offset in the gamma distribution (the "3rd" parameter)
            double ti_off = 3.0;
            // shape of the gamma distribution (the "2nd" parameter)
            double ti_shp = inputParams->shapeParameterForTheTopoIndexGammaDistribution;
            // chi -- loglamb is the first parameter (mean)
            double ti_chi = (inputParams->meanValueOfTheLogTransformedTopographicIndex - ti_off)
                            / inputParams->shapeParameterForTheTopoIndexGammaDistribution;
            // argument to the incomplete gamma function
            double ti_arg = std::max(0.0, (ti_log - ti_off)) / ti_chi;
            // pgamma is the incomplete gamma function # fortran version: pgamma(scale,shape)
            gamma_distribution<> gd(ti_shp);
            *saturatedArea = 1.0 - cdf(gd, ti_arg);
            // We originally had : *saturatedArea = 1.0 - pgamma(ti_arg, ti_shp);
        }
    }

    // compute surface runoff
    *surfaceRunoff = effectiveRainfall * *saturatedArea;

}




SmaStateStructure initStates(ModelStructure smodl,
                             SmaInputParameters inputParams,
                             SmaDerivedParameters derivedParams,
                             double fracstate0){

    /*
        Initialize model states (upper layer and lower layer)
        Author: Claudia Vitolo

        Args:
          smodl:        model structure
          inputParams:  initial model parameters
          dparam:       derived model parameters
          fracstate0:   arbitrary value to set initial state variables (default = 0.25)

        Returns:
          Initial condition for state variables
    */

    double xmin = 1e-06;

    SmaStateStructure ic;

    // initialize model states (upper layer and lower layer)
    ic.excessTensionStorageUpperLayer           = -999.0; // Excess Tension storage Upper Layer
    ic.rechargeTensionStorageUpperLayer         = -999.0; // Recharge Tension storage Upper Layer
    ic.tensionStorageUpperLayer                 = -999.0; // Tension storage Upper Layer
    ic.freeStorageUpperLayer                    = -999.0; // Free storage Upper Layer
    ic.totalUpperLayerStorage                   = -999.0; // Total Upper Layer Storage
    ic.tensionStorageLowerLayer                 = -999.0; // Tension storage Lower Layer
    ic.freeStoragePrimaryBaseflowReservoir      = -999.0; // Free Storage Primary Baseflow Reservoir
    ic.freeStorageSecondaryBaseflowReservoir    = -999.0; // Free Storage Secondary Baseflow Reservoir
    ic.totalLowerLayerStorage                   = -999.0; // Total Lower Layer Storage
    ic.freeStorageBaseflowReservoir             = -999.0; // Free Storage Baseflow Reservoir

    // UPPER LAYER *****************************************************************
    if(smodl.upperSoilLayerArch == SingleState){
        // single state
        ic.totalUpperLayerStorage = inputParams.maximumTotalStorageInLayer1 * fracstate0;
        ic.freeStorageUpperLayer = derivedParams.maximumFreeStorageInUpperLayer * fracstate0;
        ic.tensionStorageUpperLayer = derivedParams.maximumTensionStorageInUpperLayer * fracstate0;
    }else{
        // upper layer broken up into multpiple storages
        ic.freeStorageUpperLayer = derivedParams.maximumFreeStorageInUpperLayer * fracstate0;
        ic.excessTensionStorageUpperLayer = derivedParams.maximumStorageInThePrimaryTensionReservoir * fracstate0;
        ic.rechargeTensionStorageUpperLayer = derivedParams.maximumStorageInTheSecondaryTensionReservoir * fracstate0;
        ic.tensionStorageUpperLayer = derivedParams.maximumTensionStorageInUpperLayer * fracstate0;
        ic.totalUpperLayerStorage = inputParams.maximumTotalStorageInLayer1 * fracstate0;
    }

    // LOWER LAYER ******************************************************************
    if(smodl.lowerSoilLayerArch == TensionReservoirPlusTwoParallelTanks){
        // tension reservoir plus two parallel tanks
        ic.freeStoragePrimaryBaseflowReservoir = derivedParams.maximumStorageInThePrimaryBaseflowReservoir * fracstate0;
        ic.freeStorageSecondaryBaseflowReservoir = derivedParams.maximumStorageInTheSecondaryBaseflowReservoir * fracstate0;
        ic.freeStorageBaseflowReservoir = derivedParams.maximumFreeStorageInLowerLayer * fracstate0;
        ic.totalLowerLayerStorage = inputParams.maximumTotalStorageInLayer2 * fracstate0;
        ic.tensionStorageLowerLayer = derivedParams.maximumTensionStorageInLowerLayer * fracstate0;
    }else{
        // single state
        ic.totalLowerLayerStorage = inputParams.maximumTotalStorageInLayer2 * fracstate0;
        ic.freeStorageBaseflowReservoir = derivedParams.maximumFreeStorageInLowerLayer * fracstate0;
        ic.tensionStorageLowerLayer = derivedParams.maximumTensionStorageInLowerLayer * fracstate0;
    }

    return ic;

}




double logismooth( double stateX, double maxStateX ){

    /*
        Uses a logistic function to smooth the threshold at the top of a bucket
        Author: Claudia Vitolo

        Args:
          statex:
          maxstatex:

        Returns:
          Threshold at the top of a bucket
    */

    // smoothing parameter
    double psmooth = 0.01;
    // actual smoothing
    double asmooth = psmooth * maxStateX;
    return 1.0 / ( 1.0 + std::exp( - (stateX - (maxStateX - asmooth * 5.0) ) / asmooth) );
}




void mDifEqn( state_type &x,
              state_type &dxdt,
              double t ){

    /*
        Model's Differential Equations
        Author: Claudia Vitolo

        Args:
          t:                             time
          state:                         model states at time "t"
          parameters:                    parameters for differential equations: model parameters + derived model parameters
          ppt:                           rain+snow melt time series
          pet:                           potential evapotranspiration time series
          smodl:                         list of model components

        Returns:
          List of parameters.
    */

    /*
    globalInputParams = &inputParams;
    globalDerivedParams = &derivedParams;
    globalPrecipitation = precipitation;
    globalPrecipitationLength = precipitationLength;
    globalPotentialEvapotrans = potentialEvapotrans;
    globalPotentialEvapotransLength = potentialEvapotransLength;
    globalSmodl = &smodl;
    */

    upstates(x);

    // cout << "In mDifEqn at t=" << t << endl;

    // compute fluxes:
    // cout << "In mdifeqn" << endl;
    SmaFluxStructure m_flux = computeFluxes(globalSmodl,
                                            &globalPrecipitation[ globalTimeCounter ],
                                            &globalPotentialEvapotrans[ globalTimeCounter ],
                                            globalInputParams,
                                            globalDerivedParams,
                                            &x);

    // initialize derivatives (rate of change) of all states
    double dtens_1a = 0.0;
    double dtens_1b = 0.0;
    double dfree_1 = 0.0;
    double dtens_1 = 0.0;
    double dwatr_1 = 0.0;
    double dtens_2 = 0.0;
    double dfree_2a = 0.0;
    double dfree_2b = 0.0;
    double dwatr_2 = 0.0;
    double dfree_2 = 0.0;

    // compute derivatives for states in the upper layer*********************************************************************
    // upper layer defined by a single state variable
    if(globalSmodl->upperSoilLayerArch == SingleState){
        dwatr_1 = m_flux.effectiveRainfall - m_flux.surfaceRunoff - m_flux.evaporationFromTheUpperLayer
                  - m_flux.percolationFromTheUpperToLowerSoilLayers - m_flux.interflowFromFreeWaterInTheUpperLayer
                  - m_flux.overflowFromTheUpperSoilLayer;

    }
    // upper layer broken up into tension and free storage
    if(globalSmodl->upperSoilLayerArch == SeparateTensionStorage){
        dtens_1 = m_flux.effectiveRainfall - m_flux.surfaceRunoff - m_flux.evaporationFromTheUpperLayer
                  - m_flux.flowFromTensionStorageToFreeStorageInTheUpperLayer;
        dfree_1 = m_flux.flowFromTensionStorageToFreeStorageInTheUpperLayer - m_flux.percolationFromTheUpperToLowerSoilLayers
                  - m_flux.interflowFromFreeWaterInTheUpperLayer - m_flux.overflowFromTheUpperSoilLayer;
    }
    // tension storage sub-divided into recharge and excess
    if(globalSmodl->upperSoilLayerArch == CascadingBuckets){
        dtens_1a = m_flux.effectiveRainfall - m_flux.surfaceRunoff - m_flux.evaporationFromTheUpperLayerRecharge
                   - m_flux.flowFromRechargeToExcess;
        dtens_1b = m_flux.flowFromRechargeToExcess - m_flux.evaporationFromTheUpperLayerExcess
                   - m_flux.flowFromTensionStorageToFreeStorageInTheUpperLayer;
        dfree_1 = m_flux.flowFromTensionStorageToFreeStorageInTheUpperLayer - m_flux.percolationFromTheUpperToLowerSoilLayers
                  - m_flux.interflowFromFreeWaterInTheUpperLayer - m_flux.overflowFromTheUpperSoilLayer;
    }

    // compute derivatives for states in the lower layer************************************************************
    // single state
    if(globalSmodl->lowerSoilLayerArch == FixedSizeBaseflowReservoir ||
       globalSmodl->lowerSoilLayerArch == UnlimitedSizeBaseflowReservoirFracRate ||
       globalSmodl->lowerSoilLayerArch == UnlimitedSizeBaseflowReservoirPowerRecession ){

        dwatr_2 = m_flux.percolationFromTheUpperToLowerSoilLayers - m_flux.evaporationFromTheLowerLayer
                  - m_flux.totalBaseflow - m_flux.overflowFromTheLowerSoilLayer;

    }
    // tension reservoir plus two parallel tanks
    if(globalSmodl->lowerSoilLayerArch == TensionReservoirPlusTwoParallelTanks){
        dtens_2 = m_flux.percolationFromTheUpperToLowerSoilLayers * (1-globalInputParams->fractionOfPercolationToTensionStorage)
                  - m_flux.evaporationFromTheLowerLayer - m_flux.flowFromTensionStorageToFreeStorageInTheLowerLayer;
        dfree_2a = m_flux.percolationFromTheUpperToLowerSoilLayers * (globalInputParams->fractionOfPercolationToTensionStorage / 2.0)
                  + m_flux.flowFromTensionStorageToFreeStorageInTheLowerLayer / 2.0
                  - m_flux.fractionOfBaseflowInThePrimaryBaseflowReservoir - m_flux.overflowFromTheLowerSoilLayerPrimary;
        dfree_2b = m_flux.percolationFromTheUpperToLowerSoilLayers * (globalInputParams->fractionOfPercolationToTensionStorage / 2.0)
                  + m_flux.flowFromTensionStorageToFreeStorageInTheLowerLayer / 2.0
                  - m_flux.fractionOfBaseflowInTheSecondaryBaseflowReservoir - m_flux.overflowFromTheLowerSoilLayerSecondary;
    }

    dxdt[0] = dtens_1a;
    dxdt[1] = dtens_1b;
    dxdt[2] = dtens_1;
    dxdt[3] = dfree_1;
    dxdt[4] = dwatr_1;
    dxdt[5] = dtens_2;
    dxdt[6] = dfree_2a;
    dxdt[7] = dfree_2b;
    dxdt[8] = dwatr_2;
    dxdt[9] = dfree_2;


    /*
      Example code:
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
    */

}




void meanTipow(double shapeParameterForTheTopoIndexGammaDistribution,
                double meanValueOfTheLogTransformedTopographicIndex,
                double baseflowExponent,
                double* powerTransformedTopographicIndex,
                double* meanOfThePowerTransformedTopographicIndex){

    /*
        Computes the mean of the power-transformed topographic index
        Author: Claudia Vitolo

        Args:
          ti_shp:
          log_lamb:
          qb_pwr:

        Returns:
          Mean of the power-transformed topographic index.
    */

    // internal variables
    int nbins = 2000;   // number of bins in pdf of topo index
    double ti_max = 50; // maximum possible log-transformed index

    // preliminaries: get parameters of the gamma distribution [ti_shp = shape of the gamma distribution  (the "2nd" parameter)]
    double ti_off = 3;                                 // offset in the gamma distribution (the "3rd" parameter)
    double ti_chi = (meanValueOfTheLogTransformedTopographicIndex - ti_off) / shapeParameterForTheTopoIndexGammaDistribution;    // chi -- loglamb is the first parameter (mean)

    // loop through the frequency distribution
    double lowerv = 0.0;
    double lowerp = 0.0;
    double avepow = 0.0;
    double powval = 0.0;

    gamma_distribution<> gd(shapeParameterForTheTopoIndexGammaDistribution);

    double maxGmarge = 0.0;

    for(int ibin=1; ibin<=nbins; ibin++){
        // get probability for the current bin
        double upperv = (double(ibin)/double(nbins)) * ti_max;      // upper value in frequency bin
        double gmarg2 = std::max(0.0, upperv - ti_off) / ti_chi;      // 1st argument to the gamma function
        double upperp = cdf(gd, gmarg2);                            // gammp is the incomplete gamma function
        double probin = upperp - lowerp;                            // probability of the current bin
        // get the scaled topographic index value
        double logval = 0.5 * (lowerv + upperv);                    // log-transformed index for the current bin
        powval = pow( exp(logval), (1.0/baseflowExponent));    // power-transformed index for the current bin
        avepow = avepow + powval * probin;                   // average power-transformed index
        // save the lower value and probability
        lowerv = upperv;                                            // lower value for the next bin
        lowerp = upperp;                                            // cumulative probability for the next bin
    }

    *powerTransformedTopographicIndex = powval;
    *meanOfThePowerTransformedTopographicIndex = avepow;

    return;

}




void outFluxes(ModelStructure* smodl,
               double* precipitation,
               int precipitationLength,
               double* potentialEvapotrans,
               int potentialEvapotransLength,
               SmaInputParameters* inputParams,
               SmaDerivedParameters* derivedParams,
               double* s,

               double* fluxes,
               bool correctNegativeVolumes){
    /*
        Output Fluxes
        Author: Claudia Vitolo

        Args:
          smodl:                        list of model components
          P:                            rain+snow melt time series
          E:                            potential evapotranspiration time series
          inputParams:                  model parameters
          dparam:                       derived model parameters
          state:                        state variables

        Returns:
          All the fluxes + effective rainfall (or instantaneous discharge)
      */

    // set all fluxes to zero at the start of time step (including the effective rainfall)
    for(int i=0; i<20; i++){
        for(int j=0; j<precipitationLength; j++){
            fluxes[j*20 + i] = 0.0;
        }
    }

    for(int i=0; i<precipitationLength; i++){

        state_type st( 10 );

        st[0] = s[i*10 + 0];
        st[1] = s[i*10 + 1];
        st[2] = s[i*10 + 2];
        st[3] = s[i*10 + 3];
        st[4] = s[i*10 + 4];
        st[5] = s[i*10 + 5];
        st[6] = s[i*10 + 6];
        st[7] = s[i*10 + 7];
        st[8] = s[i*10 + 8];
        st[9] = s[i*10 + 9];

        SmaFluxStructure fs = computeFluxes(smodl,
                                            &precipitation[i],
                                            &potentialEvapotrans[i],
                                            inputParams,
                                            derivedParams,
                                            &st);

        // eff_ppt
        fluxes[i*20 + 0] = fs.effectiveRainfall;

        // satarea
        fluxes[i*20 + 1] = fs.saturatedArea;

        // qrunoff
        fluxes[i*20 + 2] = fs.surfaceRunoff;

        // evap_1a
        fluxes[i*20 + 3] = fs.evaporationFromTheUpperLayerRecharge;

        // evap_1b
        fluxes[i*20 + 4] = fs.evaporationFromTheUpperLayerExcess;

        // evap_1
        fluxes[i*20 + 5] = fs.evaporationFromTheUpperLayer;

        // evap_2
        fluxes[i*20 + 6] = fs.evaporationFromTheLowerLayer;

        // rchr2excs
        fluxes[i*20 + 7] = fs.flowFromRechargeToExcess;

        // tens2free_1
        fluxes[i*20 + 8] = fs.flowFromTensionStorageToFreeStorageInTheUpperLayer;

        // tens2free_2
        fluxes[i*20 + 9] = fs.flowFromTensionStorageToFreeStorageInTheLowerLayer;

        // qintf_1
        fluxes[i*20 + 10] = fs.interflowFromFreeWaterInTheUpperLayer;

        // qperc_12
        fluxes[i*20 + 11] = fs.percolationFromTheUpperToLowerSoilLayers;

        // qbase_2
        fluxes[i*20 + 12] = fs.totalBaseflow;

        // qbase_2a
        fluxes[i*20 + 13] = fs.fractionOfBaseflowInThePrimaryBaseflowReservoir;

        // qbase_2b
        fluxes[i*20 + 14] = fs.fractionOfBaseflowInTheSecondaryBaseflowReservoir;

        // oflow_1
        fluxes[i*20 + 15] = fs.overflowFromTheUpperSoilLayer;

        // oflow_2
        fluxes[i*20 + 16] = fs.overflowFromTheLowerSoilLayer;

        // oflow_2a
        fluxes[i*20 + 17] = fs.overflowFromTheLowerSoilLayerPrimary;

        // oflow_2b
        fluxes[i*20 + 18] = fs.overflowFromTheLowerSoilLayerSecondary;

        // U
        fluxes[i*20 +19] = fs.surfaceRunoff +
                        fs.overflowFromTheUpperSoilLayer +
                        fs.interflowFromFreeWaterInTheUpperLayer +
                        fs.overflowFromTheLowerSoilLayer +
                        fs.totalBaseflow;

        double U = fs.surfaceRunoff +
                   fs.overflowFromTheUpperSoilLayer +
                   fs.interflowFromFreeWaterInTheUpperLayer +
                   fs.overflowFromTheLowerSoilLayer +
                   fs.totalBaseflow;

        if(correctNegativeVolumes && U < 0.0){
            double totalAbsoluteVolume = fabs(fs.surfaceRunoff) +
                                         fabs(fs.overflowFromTheUpperSoilLayer) +
                                         fabs(fs.interflowFromFreeWaterInTheUpperLayer) +
                                         fabs(fs.overflowFromTheLowerSoilLayer) +
                                         fabs(fs.totalBaseflow);
            fluxes[i*20 +2] = fs.surfaceRunoff - (fs.surfaceRunoff / totalAbsoluteVolume) * U;
            fluxes[i*20 +15] = fs.overflowFromTheUpperSoilLayer - (fs.overflowFromTheUpperSoilLayer / totalAbsoluteVolume) * U;
            fluxes[i*20 +10] = fs.interflowFromFreeWaterInTheUpperLayer - (fs.interflowFromFreeWaterInTheUpperLayer / totalAbsoluteVolume) * U;
            fluxes[i*20 +16] = fs.overflowFromTheLowerSoilLayer - (fs.overflowFromTheLowerSoilLayer / totalAbsoluteVolume) * U;
            fluxes[i*20 +12] = fs.totalBaseflow - (fs.totalBaseflow / totalAbsoluteVolume) * U;
            fluxes[i*20 +19] = 0.0;
            cout << "Corrected a negative volume of " << U << endl;
        }else{
            fluxes[i*20 +19] = U;
        }

    }
}




double* qTimeDelay(float timeDelay, float timeStep, int sizeFracFuture){

    double* fracFuture = new double[sizeFracFuture];
    for(int i=0; i<sizeFracFuture; i++){
        fracFuture[i] = 0.0;
    }

    float alpha = 2.5;                  // shape parameter
    float alamb = alpha / timeDelay;    // scale parameter
    double psave = 0.0;

    float shape = 2.5;
    gamma_distribution<> gd(shape);

    // loop through time steps and compute the fraction of runoff in future time steps
    for(int i=0; i<sizeFracFuture; i++){
        float tfuture = float(i+1) * timeStep;            // future time (units of days)
        double cumprob = cdf(gd, (alamb*tfuture));       // cumulative probability at jtim
        fracFuture[i] = std::max(double(0.0), (cumprob-psave)); // probability between jtim-1 and jtim
        psave = cumprob;                                // cumulative probability at jtim-1
    }

    // Determine sum
    double fracFutureSum = 0.0;
    for(int i=0; i<sizeFracFuture; i++){
        fracFutureSum += fracFuture[i];
    }

    double oneOverFracFutureSum = 1.0 / fracFutureSum;

    for(int i=0; i<sizeFracFuture; i++){
        fracFuture[i] *= oneOverFracFutureSum;
    }

    return fracFuture;

}




void upstates(state_type &x){

    double xmin = 1e-06;

    double tens_1a = -999.0;
    double tens_1b = -999.0;
    double tens_1 = -999.0;
    double free_1 = -999.0;
    double watr_1 = -999.0;
    double tens_2 = -999.0;
    double free_2a = -999.0;
    double free_2b = -999.0;
    double watr_2 = -999.0;
    double free_2 = -999.0;

    // UPPER LAYER------------------------------------------------------------------------
    if (globalSmodl->upperSoilLayerArch == SingleState) {
        // case('onestate_1'): upper layer defined by a single state variable
        // (update state)
        watr_1 = std::max( xmin * globalInputParams->maximumTotalStorageInLayer1, x[4] );           // total storage
        // (derive state)
        free_1 = std::max( xmin, x[4] - globalDerivedParams->maximumTensionStorageInUpperLayer );   // free storage
        tens_1 = std::min( x[4], globalDerivedParams->maximumTensionStorageInUpperLayer );          // tension storage
    }
    if (globalSmodl->upperSoilLayerArch == SeparateTensionStorage) {
        // case('tension1_1'): upper layer broken up into tension and free storage
        // (update state)
        tens_1 = std::max( xmin * globalDerivedParams->maximumTensionStorageInUpperLayer, x[2] );   // tension storage
        free_1 = std::max( xmin * globalDerivedParams->maximumFreeStorageInUpperLayer, x[3] );      // free storage
        // (derive state)
        watr_1 = tens_1 + free_1;   // total storage
    }
    if (globalSmodl->upperSoilLayerArch == CascadingBuckets) {
        // case('tension2_1') // tension storage sub-divided into recharge and excess
        // (update state)
        // 1st tension store
        tens_1a = std::max( xmin * globalDerivedParams->maximumStorageInThePrimaryTensionReservoir, x[0] );
        // 2nd tension store
        tens_1b = std::max( xmin * globalDerivedParams->maximumStorageInTheSecondaryTensionReservoir, x[1] );
        // free storage
        free_1 = std::max( xmin * globalDerivedParams->maximumFreeStorageInUpperLayer, x[3] );
        // (derive state)
        tens_1 = tens_1a + tens_1b; // tension storage
        watr_1 = tens_1  + free_1;  // total storage
    }

    // LOWER LAYER------------------------------------------------------------------------
    if(     globalSmodl->lowerSoilLayerArch == FixedSizeBaseflowReservoir ||
            globalSmodl->lowerSoilLayerArch == UnlimitedSizeBaseflowReservoirFracRate ||
            globalSmodl->lowerSoilLayerArch == UnlimitedSizeBaseflowReservoirPowerRecession) {
        //case('unlimfrc_2','unlimpow_2','topmdexp_2','fixedsiz_2') single baseflow reservoir
        // (update state)
        watr_2  = std::max( xmin * globalInputParams->maximumTotalStorageInLayer2, x[8] );          // total storage
        // (derive state)
        free_2  = std::max( 0.0, x[8] - globalDerivedParams->maximumTensionStorageInLowerLayer );   // free storage
        tens_2  = std::min( x[8], globalDerivedParams->maximumTensionStorageInLowerLayer );         // tension storage
    }
    if (globalSmodl->lowerSoilLayerArch == TensionReservoirPlusTwoParallelTanks) {
        //case('tens2pll_2') // tension reservoir plus two parallel tanks
        // (update state)
        // primary reservoir
        free_2a = std::max( xmin * globalDerivedParams->maximumStorageInThePrimaryBaseflowReservoir, x[6] );
        // secondary reservoir
        free_2b = std::max( xmin * globalDerivedParams->maximumStorageInTheSecondaryBaseflowReservoir, x[7] );
        // tension storage
        tens_2  = std::max( xmin * globalDerivedParams->maximumTensionStorageInLowerLayer, x[5] );
        // (derive state)
        free_2  = free_2a + free_2b;                           // free storage
        watr_2  = tens_2 + free_2;                             // total storage
    }

    // Populate the state
    x[0] = tens_1a;
    x[1] = tens_1b;
    x[2] = tens_1;
    x[3] = free_1;
    x[4] = watr_1;
    x[5] = tens_2;
    x[6] = free_2a;
    x[7] = free_2b;
    x[8] = watr_2;
    x[9] = free_2;

    return;

}

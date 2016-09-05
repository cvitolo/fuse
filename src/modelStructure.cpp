// Use CRAN package BH for Boost headers

// [[fuse::depends(BH)]]
#include "modelStructure.h"

#include <fstream>
#include <vector>
#include <sstream>
#include <string.h>
#include <stdlib.h>

using std::vector;
using std::ifstream;
using std::string;
using std::stringstream;
using std::getline;

bool getModelStructures(const char* fileName, vector<ModelStructure>* modelStructures){

    bool errorsEncountered = false;
    string item;
    string line;
    vector<string> columns;

    char delim;
    strncpy(&delim, ",", 1);

    modelStructures->clear();   // Clear out any existing model structure data

    ifstream inputFile(fileName);
    if( ! inputFile.is_open() )
        return false;

    int anticipatedIndex = 1;
    getline(inputFile, line); // Eat the header
    while( inputFile.good() ){
        getline(inputFile, line);
        stringstream ss(line);
        columns.clear();
        while( getline(ss, item, delim) ){
            columns.push_back(item);
        }
        if( columns.size() == 0 )
            continue;
        if( atoi( columns.at(1).c_str() ) != anticipatedIndex ){
            errorsEncountered = true;
            break;
        }

        ModelStructure ms;
        ms.evaporation          = (EvaporationType)atoi(columns.at(7).c_str());
        ms.interflows           = (InterflowsType)atoi(columns.at(8).c_str());
        ms.lowerSoilLayerArch   = (LowerSoilLayerArchType)atoi(columns.at(4).c_str());
        ms.percolation          = (PercolationType)atoi(columns.at(6).c_str());
        ms.rainfallError        = (RainfallErrorType)atoi(columns.at(2).c_str());
        ms.routing              = (RoutingType)atoi(columns.at(9).c_str());
        ms.runoff               = (RunoffType)atoi(columns.at(5).c_str());
        ms.upperSoilLayerArch   = (UpperSoilLayerArchType)atoi(columns.at(3).c_str());

        modelStructures->push_back(ms);
        anticipatedIndex += 1;
    }

    inputFile.close();

    /*
    bool debug = true;
    if(debug){
        cout << "Model structure is:" << endl;
        for(int i=0; i<modelStructures->size(); i++){
            ModelStructure m = modelStructures->at(i);
            cout << i << ":" << endl;
            cout << "routing" << m.routing << endl;
        }
    }
    */

    if(errorsEncountered){
        modelStructures->clear();
        return false;
    }else{
        return true;
    }
}

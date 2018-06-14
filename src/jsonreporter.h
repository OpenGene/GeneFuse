#ifndef JSON_REPORTER_H
#define JSON_REPORTER_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "read.h"
#include "fusion.h"
#include "match.h"
#include <iostream>
#include <fstream>
#include "fusionmapper.h"

using namespace std;

class JsonReporter{
public:
    JsonReporter(string filename, FusionMapper* mapper);
    ~JsonReporter();
    void run();

private:
    string mFilename;
    FusionMapper* mFusionMapper;
    ofstream mFile;
    vector<FusionResult> mFusionResults;
};

#endif
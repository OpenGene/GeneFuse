#ifndef HTML_REPORTER_H
#define HTML_REPORTER_H

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

class HtmlReporter{
public:
    HtmlReporter(string filename, FusionMapper* mapper);
    ~HtmlReporter();
    void run();

private:
    void printHeader();
    void printCSS();
    void printJS();
    void printFooter();
    void printHelper();
    void printFusions();
    void printFusion(int id, FusionResult& fusion);
    void printScanTargets();

private:
    string mFilename;
    FusionMapper* mFusionMapper;
    ofstream mFile;
    vector<FusionResult> mFusionResults;
};


#endif
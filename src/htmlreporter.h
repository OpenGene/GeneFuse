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

using namespace std;

class HtmlReporter{
public:
    HtmlReporter(string filename, vector<Fusion>& fusionList, vector<Match*> *fusionMatches);
    ~HtmlReporter();
    void run();

private:
    void printHeader();
    void printCSS();
    void printJS();
    void printFooter();
    void printHelper();
    void printFusions();
    void printFusion(int id, Fusion& fusion, vector<Match*>& matches);
    void printScanTargets();

private:
    string mFilename;
    vector<Fusion> mFusionList;
    vector<Match*>* mFusionMatches;
    ofstream mFile;
};


#endif
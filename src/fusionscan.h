#ifndef FUSION_SCAN_H
#define FUSION_SCAN_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "read.h"
#include "fusion.h"
#include "indexer.h"

using namespace std;

class FusionScan{
public:
    FusionScan(string fusionFile, string refFile, string read1File, string read2File, string html, string json, int threadNum);
    bool scan();

private:
    string mFusionFile;
    string mRead1File;
    string mRead2File;
    string mHtmlFile;
    string mJsonFile;
    string mRefFile;
    int mThreadNum;
};


#endif
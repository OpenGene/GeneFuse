#include "fusionscan.h"
#include "fastqreader.h"
#include <iostream>
#include "htmlreporter.h"
#include "sescanner.h"
#include "pescanner.h"
#include "util.h"

FusionScan::FusionScan(string fusionFile, string refFile, string read1File, string read2File, string html, int threadNum){
    mRead1File = read1File;
    mRead2File = read2File;
    mFusionFile = fusionFile;
    mRefFile = refFile;
    mHtmlFile = html;
    mThreadNum = threadNum;
}

bool FusionScan::scan(){
    vector<Fusion> fusions = Fusion::parseCsv(mFusionFile);
}

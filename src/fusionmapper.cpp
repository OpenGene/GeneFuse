#include "fusionmapper.h"
#include "editdistance.h"
#include <iostream>
#include <fstream>
#include "util.h"
#include <string.h>
#include <sstream>
#include "globalsettings.h"

FusionMapper::FusionMapper(string refFile, vector<Fusion>& fusions){
	mRefFile = refFile;
    mFusions = fusions;
    init();
}

FusionMapper::~FusionMapper(){
    if(mIndexer != NULL){
        delete mIndexer;
        mIndexer = NULL;
    }
}


void FusionMapper::init(){
    mIndexer = new Indexer(mRefFile, mFusions);
    mIndexer->makeIndex();
}

Match* FusionMapper::mapRead(Read* r, int distanceReq, int qualReq) {
    return NULL;
}
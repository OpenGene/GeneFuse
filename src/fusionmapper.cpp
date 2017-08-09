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
    vector<SeqMatch> mapping = mIndexer->mapRead(r);
    
    //we only focus on the reads that can be mapped to two genome positions
    if(mapping.size() < 2)
        return NULL;

    //if the left part of mapping result is reverse, use its reverse complement alternative and skip this one
    if(!mIndexer->leftIsForward(mapping)) {
        return NULL;
    }

    cout<<r->mName<<endl;
    cout<<r->mSeq.mStr<<endl;
    cout << mapping.size() << " mappings " << endl;
    vector<SeqMatch>::iterator iter;
    for(iter = mapping.begin(); iter!=mapping.end(); iter++){
        iter->print();
        cout << endl;
    }
    cout << endl;

    return NULL;
}
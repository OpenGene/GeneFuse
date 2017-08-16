#include "fusionmapper.h"
#include "editdistance.h"
#include <iostream>
#include <fstream>
#include "util.h"
#include <string.h>
#include <sstream>
#include "globalsettings.h"
#include "matcher.h"

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

FastaReader* FusionMapper::getRef() {
    if(mIndexer == NULL)
        return NULL;
    else
        return mIndexer->getRef();
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

    // TODO: set readBreak, refBreak, distance
    Match* m = new Match(r, 0, 0, 0);
    return m;
}

void FusionMapper::removeAlignables(vector<Match*> *fusionMatches, int size) {
    FastaReader* ref = getRef();
    if(ref == NULL)
        return ;

    vector<Sequence> seqs;

    // first pass to gather all sequences
    for(int i=0; i<size; i++) {
        for(int m=0; m< fusionMatches[i].size(); m++) {
            seqs.push_back(fusionMatches[i][m]->getRead()->mSeq);
        }
    }

    Matcher matcher(ref, seqs);

    int removed = 0;
    // second pass to remove alignable sequences
    for(int i=0; i<size; i++) {
        for(int m=fusionMatches[i].size()-1 ;m>=0; m--) {
            MatchResult* mr = matcher.match(fusionMatches[i][m]->getRead()->mSeq);
            if(mr != NULL) {
                fusionMatches[i][m]->getRead()->mSeq.print();
                cout<<endl;
                mr->print();
                delete fusionMatches[i][m];
                fusionMatches[i].erase(fusionMatches[i].begin() + m);
                removed++;
            }
        }
    }

    cout << "sequence number before removeAlignables: " << seqs.size() << endl;
    cout << "removed: "<< removed << endl;
}
#include "fusionmapper.h"
#include "editdistance.h"
#include <iostream>
#include <fstream>
#include "util.h"
#include <string.h>
#include <sstream>
#include "globalsettings.h"
#include "matcher.h"
#include <stdlib.h>

FusionMapper::FusionMapper(string refFile, string fusionFile){
    mRefFile = refFile;
    fusionList = Fusion::parseCsv(fusionFile);
    init();
}

FusionMapper::~FusionMapper(){
    if(mIndexer != NULL){
        delete mIndexer;
        mIndexer = NULL;
    }
    if(fusionMatches!=NULL) {
        //delete fusionMatches;
        //fusionMatches = NULL;
    }
}


void FusionMapper::init(){
    mIndexer = new Indexer(mRefFile, fusionList);
    mIndexer->makeIndex();

    mFusionMatchSize = fusionList.size() * fusionList.size();

    fusionMatches = new vector<Match*>[mFusionMatchSize];
    for(int i=0;i<mFusionMatchSize;i++){
        fusionMatches[i] = vector<Match*>();
    }
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
    if(!mIndexer->inRequiredDirection(mapping)) {
        return NULL;
    }

    /*cout<<r->mName<<endl;
    cout<<r->mSeq.mStr<<endl;
    cout << mapping.size() << " mappings " << endl;
    vector<SeqMatch>::iterator iter;
    for(iter = mapping.begin(); iter!=mapping.end(); iter++){
        iter->print();
        cout << endl;
    }
    cout << endl;*/

    // TODO: set int readBreak, int leftContig, int leftPos, int rightContig, int rightPos
    Match* m = makeMatch(r, mapping);
    return m;
}

Match* FusionMapper::makeMatch(Read* r, vector<SeqMatch>& mapping) {
    if(mapping.size()!=2)
        return NULL;
    SeqMatch left = mapping[0];
    SeqMatch right = mapping[1];
    if(left.seqStart > right.seqStart) {
        left = mapping[1];
        right = mapping[0];
    }
    int readBreak = (left.seqEnd + right.seqStart)/2;
    GenePos leftGP = left.startGP;
    GenePos rightGP = right.startGP;
    leftGP.position += readBreak;
    rightGP.position += readBreak+1;
    int gap = right.seqStart - left.seqEnd - 1;
    return new Match(r, readBreak, leftGP, rightGP, gap);
}
    
void FusionMapper::addMatch(Match* m) {
    int leftContig = m->mLeftGP.contig;
    int rightContig = m->mRightGP.contig;
    int index = fusionList.size() * rightContig + leftContig;
    fusionMatches[index].push_back(m);
}

void FusionMapper::removeAlignables() {
    FastaReader* ref = getRef();
    if(ref == NULL)
        return ;

    vector<Sequence> seqs;

    // first pass to gather all sequences
    for(int i=0; i<mFusionMatchSize; i++) {
        for(int m=0; m< fusionMatches[i].size(); m++) {
            seqs.push_back(fusionMatches[i][m]->getRead()->mSeq);
        }
    }

    loginfo( string("sequence number before removeAlignables: ") + string(int2str(seqs.size())));

    Matcher matcher(ref, seqs);

    int removed = 0;
    // second pass to remove alignable sequences
    for(int i=0; i<mFusionMatchSize; i++) {
        for(int m=fusionMatches[i].size()-1 ;m>=0; m--) {
            MatchResult* mr = matcher.match(fusionMatches[i][m]->getRead()->mSeq);
            if(mr != NULL) {
                //fusionMatches[i][m]->getRead()->mSeq.print();
                //cout<<endl;
                //mr->print();
                delete fusionMatches[i][m];
                fusionMatches[i].erase(fusionMatches[i].begin() + m);
                removed++;
            }
        }
    }
    loginfo( string("removed: ") + string( int2str(removed )));
}

void FusionMapper::sortMatches() {
    // sort the matches to make the pileup more clear
    for(int i=0;i<mFusionMatchSize;i++){
        sort(fusionMatches[i].begin(), fusionMatches[i].end(), Match::greater); 
    }
}

void FusionMapper::freeMatches() {
    // free it
    for(int i=0;i<mFusionMatchSize;i++){
        for(int m=fusionMatches[i].size()-1 ;m>=0; m--)
            delete fusionMatches[i][m];
        fusionMatches[i].clear();
    }
}

void FusionMapper::clusterMatches() {
    for(int i=0;i<mFusionMatchSize;i++){
        vector<FusionResult> frs;
        for(int m=0 ;m<fusionMatches[i].size(); m++){
            bool found = false;
            Match* match = fusionMatches[i][m];
            for(int f=0; f<frs.size(); f++) {
                if(frs[f].support(match)){
                    frs[f].addMatch(match);
                    found = true;
                    break;
                }
            }
            if(!found) {
                FusionResult fr;
                fr.addMatch(match);
                frs.push_back(fr);
            }
        }
        for(int f=0; f<frs.size(); f++) {
            frs[f].calcFusionPoint();
            frs[f].print(fusionList);
            mFusionResults.push_back(frs[f]);
        }
    }
    cerr<<"mFusionResults size: " << mFusionResults.size() << endl;
}

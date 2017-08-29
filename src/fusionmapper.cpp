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

Match* FusionMapper::mapRead(Read* r, bool& mapable, int distanceReq, int qualReq) {
    vector<SeqMatch> mapping = mIndexer->mapRead(r);
    
    //we only focus on the reads that can be mapped to two genome positions
    if(mapping.size() < 2){
        mapable = false;
        return NULL;
    }

    mapable = true;

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
    Match* match = new Match(r, readBreak, leftGP, rightGP, gap);

    calcDistance(match);

    return match;
}

void FusionMapper::calcDistance(Match* match) {
    string seq = match->mRead->mSeq.mStr;

    int readBreak = match->mReadBreak;
    int leftLen = readBreak+1;
    int rightLen = seq.length() - (readBreak+1);

    string leftSeq = seq.substr(0, leftLen);
    string rightSeq = seq.substr(readBreak+1, rightLen);

    //Gene& leftGene = fusionList[match->mLeftGP.contig].mGene;
    //Gene& rightGene = fusionList[match->mRightGP.contig].mGene;

    match->mLeftDistance = calcED(leftSeq, match->mLeftGP.contig, match->mLeftGP.position - leftLen + 1, match->mLeftGP.position);
    match->mRightDistance = calcED(rightSeq, match->mRightGP.contig, match->mRightGP.position, match->mRightGP.position + rightLen - 1);
}


int FusionMapper::calcED(string seq, int contig, int start, int end) {
    // check start and end are in same strand
    if( (start>=0 && end<=0) || (start<=0 && end>=0) ) {
        return -1;
    }

    string& fusionSeq = mIndexer->mFusionSeq[contig];

    // check the overflow
    if(abs(start)>=fusionSeq.length() || abs(end)>=fusionSeq.length())
        return -2;

    string str = seq;
    if(start < 0) {
        Sequence s(seq);
        Sequence rc = ~s;
        str = rc.mStr;

        int tmp = start;
        start = -end;
        end = -tmp;
    }

    string refstr = fusionSeq.substr(start, end-start+1);

    return edit_distance(str.c_str(), str.length(), refstr.c_str(), refstr.length());
}
    
void FusionMapper::addMatch(Match* m) {
    int leftContig = m->mLeftGP.contig;
    int rightContig = m->mRightGP.contig;
    int index = fusionList.size() * rightContig + leftContig;
    fusionMatches[index].push_back(m);
}

void FusionMapper::filterMatches() {
    // calc the sequence number before any filtering
    int total = 0;
    for(int i=0; i<mFusionMatchSize; i++)
        total += fusionMatches[i].size();

    loginfo( string("sequence number before filtering: ") + string( int2str(total )));

    removeByComplexity();
    removeByDistance();
    removeIndels();
    removeAlignables();
}


void FusionMapper::removeByComplexity() {
    int removed = 0;
    for(int i=0; i<mFusionMatchSize; i++) {
        for(int m=fusionMatches[i].size()-1 ;m>=0; m--) {
            string seq = fusionMatches[i][m]->mRead->mSeq.mStr;
            int readBreak = fusionMatches[i][m]->mReadBreak;
            if( isLowComplexity(seq.substr(0, readBreak+1)) 
                || isLowComplexity(seq.substr(readBreak+1, seq.length() - (readBreak+1) )) ) {
                delete fusionMatches[i][m];
                fusionMatches[i].erase(fusionMatches[i].begin() + m);
                removed++;
            }
        }
    }
    loginfo( string("removeByComplexity: ") + string( int2str(removed )));
}

bool FusionMapper::isLowComplexity(string str) {
    if(str.length() < 20)
        return true;

    int diffCount = 0;
    for(int i=0;i<str.length()-1;i++) {
        if( str[i] != str[i+1])
            diffCount++;
    }

    if(diffCount < 7)
        return true;

    return false;
}

void FusionMapper::removeByDistance() {
    // diff should be less than DIFF_THRESHOLD
    const int DIFF_THRESHOLD = 5;
    int removed = 0;
    for(int i=0; i<mFusionMatchSize; i++) {
        for(int m=fusionMatches[i].size()-1 ;m>=0; m--) {
            if(fusionMatches[i][m]->mLeftDistance + fusionMatches[i][m]->mRightDistance >= DIFF_THRESHOLD) {
                delete fusionMatches[i][m];
                fusionMatches[i].erase(fusionMatches[i].begin() + m);
                removed++;
            }
        }
    }
    loginfo( string("removeByDistance: ") + string( int2str(removed )));
}

void FusionMapper::removeIndels() {
    // diff should be greather than INDEL_THRESHOLD
    const int INDEL_THRESHOLD = 50;
    int removed = 0;
    for(int i=0; i<mFusionMatchSize; i++) {
        for(int m=fusionMatches[i].size()-1 ;m>=0; m--) {
            if(fusionMatches[i][m]->mLeftGP.contig == fusionMatches[i][m]->mRightGP.contig 
                && abs(fusionMatches[i][m]->mLeftGP.position - fusionMatches[i][m]->mRightGP.position) <= INDEL_THRESHOLD) {
                delete fusionMatches[i][m];
                fusionMatches[i].erase(fusionMatches[i].begin() + m);
                removed++;
            }
        }
    }
    loginfo( string("removeIndels: ") + string( int2str(removed )));
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
    loginfo( string("removeAlignables: ") + string( int2str(removed )));
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
            frs[f].calcUnique();
            if(frs[f].mUnique >= GlobalSettings::uniqueRequirement) {
                frs[f].print(fusionList);
                mFusionResults.push_back(frs[f]);
            }
        }
    }
    loginfo("found " + int2str(mFusionResults.size()) + " fusions");
}

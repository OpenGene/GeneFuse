#include "match.h"
#include <vector>

Match::Match(Read* r, int readBreak, GenePos leftGP, GenePos rightGP, int gap, int distance, bool reversed){
    mRead = new Read(*r);
    mDistance = distance;
    mGap = gap;
    mReadBreak = readBreak;
    mLeftGP = leftGP;
    mRightGP = rightGP;
}

Match::~Match(){
    delete mRead;
    mRead = NULL;
    for(int i=0;i<mOriginalReads.size();i++){
        delete mOriginalReads[i];
        mOriginalReads[i] = NULL;
    }
}

void Match::addOriginalRead(Read* r){
    mOriginalReads.push_back(new Read(*r));
}

void Match::addOriginalPair(ReadPair* pair){
    mOriginalReads.push_back(new Read(*pair->mLeft));
    mOriginalReads.push_back(new Read(*pair->mRight));
}

void Match::print(){
    cout<<"pos: "<<mReadBreak;
    if(mReversed)
        cout<<", reverse";
    else
        cout<<", forward";
    cout<<endl;
    vector<int> breaks;
    // we use mReadBreak+1 since it means the length of (0...mReadBreak)
    breaks.push_back( mReadBreak+1 );
    mRead->printWithBreaks(breaks);
}

void Match::printHtmlTD(ofstream& file, int leftlen, int centerlen, int rightlen){
    file<<"<a title='"<<mRead->mName<<"'>";
    file<<"d:" << mDistance;
    if(mReversed)
        file<<", <--";
    else
        file<<", -->";

    file<<"</a></span>";

    vector<int> breaks;
    breaks.push_back(max(mReadBreak-leftlen, 0));
    breaks.push_back( mReadBreak );
    breaks.push_back( mReadBreak+centerlen );
    breaks.push_back( min(mReadBreak+centerlen+rightlen, mRead->length()));
    mRead->printHtmlTDWithBreaks(file, breaks);
}

void Match::printReadsToFile(ofstream& file){
    for(int i=0;i<mOriginalReads.size();i++){
        mOriginalReads[i]->printFile(file);
    }
}

void Match::setReversed(bool flag){
    mReversed = flag;
}

int Match::countUnique(vector<Match*>& matches) {
    if(matches.size()==0)
        return 0;
    int count = 1;
    Match* cur = matches[0];
    for(int i=1;i<matches.size();i++){
        Match* m = matches[i];
        if( *m > *cur || *m < *cur) {
            cur = m;
            count++;
        }
    }
    return count;
}

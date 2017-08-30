#include "match.h"
#include <vector>
#include "util.h"

Match::Match(Read* r, int readBreak, GenePos leftGP, GenePos rightGP, int gap, bool reversed){
    mRead = new Read(*r);
    mGap = gap;
    mReadBreak = readBreak;
    mLeftGP = leftGP;
    mRightGP = rightGP;
    mLeftDistance = 0;
    mRightDistance = 0;
    mOverallDistance = 0;
    mReversed = reversed;
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
    cout<<"break:"<<mReadBreak+1;
    cout<<", diff:("<<mLeftDistance<<" "<<mRightDistance<<")";
    if(mReversed)
        cout<<", read direction: reversed complement";
    else
        cout<<", read direction: original direction";
    cout << ", name: " << mRead->mName.substr(1, mRead->mName.length()-1);
    cout << endl;
    cout << mRead->mSeq.mStr.substr(0, mReadBreak+1);
    cout << " ";
    cout << mRead->mSeq.mStr.substr(mReadBreak+1, mRead->length() - (mReadBreak+1));
    cout << endl;
}

void Match::printHtmlTD(ofstream& file){
    //file<<"d:" << mDistance;
    if(mReversed)
        file<<" <--";
    else
        file<<" -->";

    file<<"</a></span>";

    file<<"</td><td>(" << int2str(mLeftDistance) << ", " << int2str(mRightDistance) << ")</td>";

    vector<int> breaks;
    breaks.push_back( mReadBreak+1 );
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

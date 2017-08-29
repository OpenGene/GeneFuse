#include "fusionresult.h"

FusionResult::FusionResult() {

}

FusionResult::~FusionResult() {
    
}

void FusionResult::addMatch(Match* m) {
    mMatches.push_back(m);
}

bool FusionResult::support(Match* m) {
    for(int i=0; i<mMatches.size(); i++) {
        if(supportSame(m, mMatches[i]))
            return true;
    }
    return false;
}

bool FusionResult::supportSame(Match* m1, Match* m2) {
    const int T=3;

    if( abs(m1->mLeftGP.position - m2->mLeftGP.position) > T)
        return false;
    if( abs(m1->mRightGP.position - m2->mRightGP.position) > T)
        return false;
    if( m1->mLeftGP.contig != m2->mLeftGP.contig)
        return false;
    if( m1->mRightGP.contig != m2->mRightGP.contig)
        return false;

    return true;
}

void FusionResult::calcFusionPoint() {
    if(mMatches.size() == 0)
        return ;

    // if we can find an exact match with 0 gap, then use it
    // else we use the mean one
    long leftTotal = 0;
    long rightTotal = 0;
    for(int i=0; i<mMatches.size(); i++) {
        Match* match = mMatches[i];
        if(match->mGap == 0){
            mLeftGP = match->mLeftGP;
            mRightGP = match->mRightGP;

            adjustFusionBreak();
            return ;
        }
        leftTotal += match->mLeftGP.position;
        rightTotal += match->mRightGP.position;
    }

    mLeftGP.contig = mMatches[0]->mLeftGP.contig;
    mLeftGP.position = leftTotal/(long)mMatches.size();
    mRightGP.contig = mMatches[0]->mRightGP.contig;
    mRightGP.position = rightTotal/(long)mMatches.size();

    adjustFusionBreak();

}
    
void FusionResult::calcUnique() {
    mUnique = 1;
    // since it is sorted, so just check every match with previous one
    for(int i=1; i<mMatches.size(); i++) {
        if(mMatches[i]->mReadBreak != mMatches[i-1]->mReadBreak ||
            mMatches[i]->mRead->length() != mMatches[i-1]->mRead->length())
            mUnique ++;
    }
}

void FusionResult::adjustFusionBreak() {
    for(int i=0; i<mMatches.size(); i++) {
        int shift = mLeftGP.position - mMatches[i]->mLeftGP.position;
        mMatches[i]->mReadBreak += shift;
        mMatches[i]->mLeftGP.position += shift;
        mMatches[i]->mRightGP.position += shift;
        /*if(shift != 0) {
            cout << "after shift:" << mMatches[i]->mRightGP.position << ", mRightGP.position:" << mRightGP.position << endl;
        }*/
    }
}

void FusionResult::print(vector<Fusion>& fusions) {
    cout << endl << "#Fusion: ";
    cout << fusions[mLeftGP.contig].pos2str(mLeftGP.position) << "_";
    cout << fusions[mRightGP.contig].pos2str(mRightGP.position) ;
    cout << " (total: " << mMatches.size() << ", unique:" << mUnique <<")";
    cout << endl;
    for(int i=0; i<mMatches.size(); i++) {
        cout << ">" << i+1 << ", ";
        mMatches[i]->print();
    }
}
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
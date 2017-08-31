#ifndef FUSIONRESULT_H
#define FUSIONRESULT_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "match.h"
#include "fusion.h"


using namespace std;

class FusionResult{
public:
    FusionResult();
    ~FusionResult();

    void print(vector<Fusion>& fusions);
    void calcFusionPoint();
    void calcUnique();
    void makeTitle(vector<Fusion>& fusions);
    void makeReference(string& refL, string& refR);
    void adjustFusionBreak();
    void addMatch(Match* m);
    bool support(Match* m);
    bool isDeletion();
    bool canBeMapped();
    bool canBeMatched(string& s1, string& s2);
    static bool supportSame(Match* m1, Match* m2);

private:
    string getRefSeq(string& ref, int start, int end);
    int calcED(Match* m, int shift, int& leftED, int& rightED);

public:
    GenePos mLeftGP;
    GenePos mRightGP;
    vector<Match*> mMatches;
    int mUnique;
    string mTitle;
    string mLeftRef;
    string mRightRef;
    string mLeftRefExt;
    string mRightRefExt;
    string mLeftPos;
    string mRightPos;
};


#endif

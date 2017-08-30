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
    void adjustFusionBreak();
    void addMatch(Match* m);
    bool support(Match* m);
    static bool supportSame(Match* m1, Match* m2);

public:
    GenePos mLeftGP;
    GenePos mRightGP;
    vector<Match*> mMatches;
    int mUnique;
    string mTitle;
};


#endif

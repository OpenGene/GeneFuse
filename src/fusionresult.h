#ifndef FUSIONRESULT_H
#define FUSIONRESULT_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "match.h"


using namespace std;

class FusionResult{
public:
    FusionResult();
    ~FusionResult();

    void print();
    void calcFusionPoint();
    void addMatch(Match* m);
    bool support(Match* m);
    static bool supportSame(Match* m1, Match* m2);

public:
    GenePos mLeftGP;
    GenePos mRightGP;
    vector<Match*> mMatches;
};


#endif

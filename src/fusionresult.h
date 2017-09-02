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
    void updateInfo(vector<Fusion>& fusions);
    void makeReference(string& refL, string& refR);
    void adjustFusionBreak();
    void addMatch(Match* m);
    bool support(Match* m);
    bool isDeletion();
    bool canBeMapped();
    bool canBeMatched(string& s1, string& s2);
    bool isQualified();
    static bool supportSame(Match* m1, Match* m2);
    void printFusionProteinHTML(ofstream& file);

private:
    string getRefSeq(string& ref, int start, int end);
    int calcED(Match* m, int shift, int& leftED, int& rightED);
    bool isLeftProteinForward();
    bool isRightProteinForward();
    void calcLeftExonIntronNumber();
    void calcRightExonIntronNumber();
    void printLeftProteinHTML(ofstream& file);
    void printRightProteinHTML(ofstream& file);
    void printExonIntronTD(ofstream& file, bool isExon, bool forward, int number, float percent, string style);

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
    Gene mLeftGene;
    Gene mRightGene;
    bool mLeftIsExon;
    bool mRightIsExon;
    int mLeftExonOrIntronID;
    int mRightExonOrIntronID;
    float mLeftExonNum;
    float mLeftIntronNum;
    float mRightExonNum;
    float mRightIntronNum;
};


#endif

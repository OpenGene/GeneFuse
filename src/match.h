#ifndef MATCH_H
#define MATCH_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "read.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <functional>
#include "common.h"

using namespace std;

/*
*                          
*  ref:                 leftGP       rightGP
* read:           ------------       ------------
*                    readBreak       readBreak+1
*/

class Match{
public:
    Match(Read* r, int readBreak, GenePos leftGP, GenePos rightGP, int gap, int distance = 0, bool reversed = false);
    ~Match();
    void print();
    void printHtmlTD(ofstream& file, int leftlen, int centerlen, int rightlen);
    void printReadsToFile(ofstream& file);
    void setReversed(bool flag);
    void addOriginalRead(Read* r);
    void addOriginalPair(ReadPair* pair);
    Read* getRead() {return mRead;}

    inline bool operator <(const Match& other) const 
    {
        return mReadBreak < other.mReadBreak || (mReadBreak == other.mReadBreak && mRead->length() > other.mRead->length());
    }
    inline bool operator >(const Match& other) const
    {
        return mReadBreak > other.mReadBreak || (mReadBreak == other.mReadBreak && mRead->length() < other.mRead->length());
    }
    inline static bool less(const Match* m1, const Match* m2)
    {
        return *m1 < *m2;
    }
    inline static bool greater(const Match* m1, const Match* m2)
    {
        return *m1 > *m2;
    }

    static int countUnique(vector<Match*>& matches);

public:
    Read* mRead;
    vector<Read*> mOriginalReads;
    int mDistance;
    // the gap between left and right segment after segmentation
    int mGap;
    bool mReversed;
    int mReadBreak;
    GenePos mLeftGP;
    GenePos mRightGP;
};


#endif
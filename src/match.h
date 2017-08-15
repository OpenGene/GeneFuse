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

using namespace std;

/*
*                          refBreak
*  ref: ---------------------- ++++++++++++++++++++++++++++++
* read:           ------------ ------------
*                         readBreak
*/

class Match{
public:
    Match(Read* r, int readBreak, int refBreak, int distance, bool reversed = false);
    ~Match();
    void print(int leftlen, int centerlen, int rightlen);
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

private:
    Read* mRead;
    vector<Read*> mOriginalReads;
    int mDistance;
    bool mReversed;
    int mReadBreak;
    int mRefBreak;
};


#endif
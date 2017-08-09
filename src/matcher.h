#ifndef MATCHER_H
#define MATCHER_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "read.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "read.h"
#include "fastareader.h"
#include <map>
#include <vector>
#include <cctype>
#include "common.h"

using namespace std;

struct MatchResult{
    GenePos startGP;
    Sequence sequence;
    bool reversed;
    vector<int> mismatches;
    inline void print() {
        cout<<startGP.contig<<":"<<startGP.position<<(reversed?" reversed":" forward")<<endl;
        for(int m=0;m<mismatches.size();m++)
            cout << mismatches[m] << " ";
        cout << endl;
    }
};

class Matcher{
public:
    Matcher(FastaReader* ref, vector<Sequence>& seqs);
    ~Matcher();
    void makeIndex();
    void indexContig(int ctg, string seq, int start);
    void printStat();

    MatchResult* mapToIndex(Sequence& seq);
    MatchResult* match(Sequence& seq);

    static inline long gp2long(const GenePos& gp);
    static inline GenePos long2gp(const long val);
    static inline unsigned int makeKmer(string & seq, int pos, bool& valid);
    static inline int base2num(char c);
    static inline GenePos shift(const GenePos& gp, int i);

    static bool test();


private:
    void makeMask(unsigned char* mask, unsigned char flag, int seqlen, int start, int kmerSize);
    void initBloomFilter(vector<Sequence>& seqs);
    void initBloomFilterWithSeq(Sequence& seq);

public:
    map<unsigned int, vector<GenePos>> mKmerPositions;
    vector<string> mContigNames;


private:
    FastaReader* mReference;
    int mUniquePos;
    int mDupePos;
    unsigned char* mBloomFilterArray;
};


#endif

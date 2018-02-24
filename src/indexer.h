#ifndef INDEXER_H
#define INDEXER_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "read.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "fusion.h"
#include <map>
#include <vector>
#include <cctype>
#include "common.h"


using namespace std;
const unsigned char MATCH_TOP = 3;
const unsigned char MATCH_SECOND = 2;
const unsigned char MATCH_NONE = 1;
const unsigned char MATCH_UNKNOWN = 0; 

#pragma pack() 

struct SeqMatch{
    int seqStart;
    int seqEnd;
    GenePos startGP;
    inline void print() {
        cout << seqStart << "-" << seqEnd << "|";
        cout << startGP.contig << ":" << startGP.position;
    }
};

class Indexer{
public:
    Indexer(string refFile, vector<Fusion>& fusions);
    ~Indexer();
    void makeIndex();
    void indexContig(int ctg, string seq, int start);
    void printStat();
    FastaReader* getRef();

    static inline long makeKmer(string & seq, int pos, long lastKmer, int step = 1);
    static inline long gp2long(const GenePos& gp);
    static inline GenePos shift(const GenePos& gp, int i);
    static inline GenePos long2gp(const long val);

    // map the read onto the reference
    // return a map, with key is the first-base-mapped GenePos encoded in long, value is the count of this GenePos
    // GenePos encoded with 0 means not mapped
    vector<SeqMatch> mapRead(Read* r);
    // check a fusion of two mapping mates is: forward at left
    bool inRequiredDirection(vector<SeqMatch>& mapping);

    static bool test();


private:
    void makeMask(unsigned char* mask, unsigned char flag, int seqlen, int start, int kmerSize);
    vector<SeqMatch> segmentMask(unsigned char* mask, int seqlen, GenePos gp1, GenePos gp2);
    void fillBloomFilter();

public:
    map<long, GenePos> mKmerPos;
    unsigned char* mBloomFilter;
    vector<vector<GenePos>> mDupeList;
    vector<string> mFusionSeq;



private:
    string mRefFile;
    FastaReader* mReference;
    vector<Fusion> mFusions;
    int mUniquePos;
    int mDupePos;
};


#endif

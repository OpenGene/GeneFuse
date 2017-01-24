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


using namespace std;

const int KMER = 16;

#pragma pack(2) 
// if contig is -1, means this is a dupe entry, and position will be the position in the dupList
struct GenePos{
    short contig;
    int position;
};
#pragma pack() 

class Indexer{
public:
    Indexer(string refFile, vector<Fusion>& fusions);
    ~Indexer();
    void makeIndex();
    void indexContig(int ctg, string seq);
    void printStat();

    static inline long makeKmer(string & seq, int pos);
    static inline long gp2long(const GenePos& gp);
    static inline GenePos long2gp(const long val);

    vector<GenePos> mapRead(Read* r);

    static bool test();

public:
    map<long, GenePos> mKmerPos;
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

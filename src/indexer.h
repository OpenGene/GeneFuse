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
    Indexer(string refFile);
    ~Indexer();
    void makeIndex(vector<Fusion>& fusions);
    void indexContig(int ctg, string seq);
    void printStat();

public:
    map<long, GenePos> mKmerPos;
    vector<vector<GenePos>> mDupeList;

private:
    inline long makeKmer(string & seq, int pos);

private:
    string mRefFile;
    FastaReader* mReference;
    int mUniquePos;
    int mDupePos;
};


#endif

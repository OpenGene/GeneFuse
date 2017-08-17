#ifndef GENE_H
#define GENE_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "read.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "fastareader.h"

using namespace std;

struct Exon{
    int id;
    int start;
    int end;
};

class Gene{
public:
    Gene(string name, string chr, int start, int end);
    Gene(const Gene& other);
    // WAR to make a default constructor
    Gene();
    bool valid();
    bool isReversed() { return mReversed; }
    void addExon(Exon exon);
    void addExon(int id, int start, int end);
    void print();
    static Gene parse(string linestr);

public:
    string mName;
    string mChr;
    int mStart;
    int mEnd;
    vector<Exon> mExons;
    bool mReversed;
};

#endif

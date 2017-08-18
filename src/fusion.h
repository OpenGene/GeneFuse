#ifndef FUSION_H
#define FUSION_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "read.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "fastareader.h"
#include "gene.h"


using namespace std;

class Fusion{
public:
    Fusion(Gene gene);

    static vector<Fusion> parseCsv(string filename);

    void print();
    void printHtml(ofstream& file);
    static bool test();
    bool isReversed() {return mGene.isReversed();}
    string pos2str(int pos);

public:
    Gene mGene;
};


#endif

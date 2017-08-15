#ifndef FUSIONMAPPER_H
#define FUSIONMAPPER_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "read.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "indexer.h"
#include "match.h"


using namespace std;

class FusionMapper{
public:
    FusionMapper(string refFile, vector<Fusion>& fusions);
    ~FusionMapper();

    Match* mapRead(Read* r, int distanceReq = 2, int qualReq=20);
    FastaReader* getRef();

    void removeAlignables(vector<Match*> *fusionMatches, int size);

private:
    void init();

public:
    string mRefFile;
    Indexer* mIndexer;
    vector<Fusion> mFusions;
};


#endif

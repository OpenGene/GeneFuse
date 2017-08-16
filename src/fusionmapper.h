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
    FusionMapper(string refFile, string fusionFile);
    ~FusionMapper();

    Match* mapRead(Read* r, int distanceReq = 2, int qualReq=20);
    FastaReader* getRef();

    void removeAlignables();
    void sortMatches();
    void freeMatches();
    void clusterMatches();
    void addMatch(Match* m);

private:
    void init();
    Match* makeMatch(Read* r, vector<SeqMatch>& mapping);

public:
    string mRefFile;
    int mFusionMatchSize;
    Indexer* mIndexer;
    vector<Fusion> fusionList;
    vector<Match*> *fusionMatches;
};


#endif

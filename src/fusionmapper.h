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
#include "fusionresult.h"


using namespace std;

class FusionMapper{
public:
    FusionMapper(string refFile, string fusionFile);
    ~FusionMapper();

    Match* mapRead(Read* r, bool& mapable, int distanceReq = 2, int qualReq=20);
    FastaReader* getRef();
    void filterMatches();
    void sortMatches();
    void freeMatches();
    void clusterMatches();
    void addMatch(Match* m);

private:
    void init();
    Match* makeMatch(Read* r, vector<SeqMatch>& mapping);
    void calcDistance(Match* match);
    int calcED(string seq, int contig, int start, int end);

    void removeAlignables();
    void removeByDistance();
    void removeIndels();
    void removeByComplexity();
    bool isLowComplexity(string str);

public:
    string mRefFile;
    int mFusionMatchSize;
    Indexer* mIndexer;
    vector<Fusion> fusionList;
    vector<Match*> *fusionMatches;
    vector<FusionResult> mFusionResults;
};


#endif

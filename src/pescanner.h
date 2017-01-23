#ifndef PE_SCANNNER_H
#define PE_SCANNNER_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "read.h"
#include "fusion.h"
#include "match.h"
#include <cstdlib>
#include <condition_variable>
#include <mutex>
#include <thread>


using namespace std;

struct ReadPairPack {
    ReadPair** data;
    int count;
};

typedef struct ReadPairPack ReadPairPack;

struct ReadPairRepository {
    ReadPairPack** packBuffer;
    size_t readPos;
    size_t writePos;
    size_t readCounter;
    std::mutex mtx;
    std::mutex readCounterMtx;
    std::condition_variable repoNotFull;
    std::condition_variable repoNotEmpty;
};

typedef struct ReadPairRepository ReadPairRepository;

class PairEndScanner{
public:
    PairEndScanner(string fusionFile, string refFile, string read1File, string read2File, string html="", int threadnum=1);
    bool scan();
    void textReport(vector<Fusion>& fusionList, vector<Match*> *fusionMatches);
    void htmlReport(vector<Fusion>& fusionList, vector<Match*> *fusionMatches);

private:
    bool scanPairEnd(ReadPairPack* pack);
    void initPackRepository();
    void destroyPackRepository();
    void producePack(ReadPairPack* pack);
    void consumePack();
    void producerTask();
    void consumerTask();
    void pushMatch(int i, Match* m);

private:
    string mFusionFile;
    string mRefFile;
    string mRead1File;
    string mRead2File;
    string mHtmlFile;
    ReadPairRepository mRepo;
    bool mProduceFinished;
    vector<Fusion> fusionList;
    vector<Match*> *fusionMatches;
    std::mutex mFusionMtx;
    int mThreadNum;
};


#endif
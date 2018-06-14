#ifndef SE_SCANNNER_H
#define SE_SCANNNER_H

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
#include "fusionmapper.h"


using namespace std;

struct ReadPack {
    Read** data;
    int count;
};

typedef struct ReadPack ReadPack;

struct ReadRepository {
    ReadPack** packBuffer;
    size_t readPos;
    size_t writePos;
    size_t readCounter;
    std::mutex mtx;
    std::mutex readCounterMtx;
    std::condition_variable repoNotFull;
    std::condition_variable repoNotEmpty;
};

typedef struct ReadRepository ReadRepository;

class SingleEndScanner{
public:
    SingleEndScanner(string fusionFile, string refFile, string read1File, string html, string json, int threadnum);
    ~SingleEndScanner();
    bool scan();
    void textReport();
    void htmlReport();
    void jsonReport();

private:
    bool scanSingleEnd(ReadPack* pack);
    void initPackRepository();
    void destroyPackRepository();
    void producePack(ReadPack* pack);
    void consumePack();
    void producerTask();
    void consumerTask();
    void pushMatch(Match* m);

private:
    string mFusionFile;
    string mRefFile;
    string mRead1File;
    string mRead2File;
    string mHtmlFile;
    string mJsonFile;
    ReadRepository mRepo;
    bool mProduceFinished;
    std::mutex mFusionMtx;
    int mThreadNum;
    FusionMapper* mFusionMapper;
};


#endif
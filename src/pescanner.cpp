#include "pescanner.h"
#include "fastqreader.h"
#include <iostream>
#include "htmlreporter.h"
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include "util.h"

PairEndScanner::PairEndScanner(string fusionFile, string refFile, string read1File, string read2File, string html, int threadNum){
    mRead1File = read1File;
    mRead2File = read2File;
    mFusionFile = fusionFile;
    mRefFile = refFile;
    mHtmlFile = html;
    mProduceFinished = false;
    mThreadNum = threadNum;
    mFusionMapper = NULL;
}

PairEndScanner::~PairEndScanner() {
    if(mFusionMapper != NULL) {
        delete mFusionMapper;
        mFusionMapper = NULL;
    }
}

bool PairEndScanner::scan(){

    mFusionMapper = new FusionMapper(mRefFile, mFusionFile);

    initPackRepository();
    std::thread producer(std::bind(&PairEndScanner::producerTask, this));

    std::thread** threads = new thread*[mThreadNum];
    for(int t=0; t<mThreadNum; t++){
        threads[t] = new std::thread(std::bind(&PairEndScanner::consumerTask, this));
    }

    producer.join();
    for(int t=0; t<mThreadNum; t++){
        threads[t]->join();
    }

    for(int t=0; t<mThreadNum; t++){
        delete threads[t];
        threads[t] = NULL;
    }

    mFusionMapper->removeAlignables();
    mFusionMapper->sortMatches();
    mFusionMapper->clusterMatches();

    textReport(mFusionMapper->fusionList, mFusionMapper->fusionMatches);
    htmlReport(mFusionMapper->fusionList, mFusionMapper->fusionMatches);

    mFusionMapper->freeMatches();
    return true;
}

void PairEndScanner::pushMatch(Match* m){
    std::unique_lock<std::mutex> lock(mFusionMtx);
    mFusionMapper->addMatch(m);
    lock.unlock();
}

bool PairEndScanner::scanPairEnd(ReadPairPack* pack){
    for(int p=0;p<pack->count;p++){
        ReadPair* pair = pack->data[p];
        Read* r1 = pair->mLeft;
        Read* r2 = pair->mRight;
        Read* rcr1 = NULL;
        Read* rcr2 = NULL;
        Read* merged = pair->fastMerge();
        Read* mergedRC = NULL;
        bool mapable = false;
        // if merged successfully, we only search the merged
        if(merged != NULL) {
            Match* matchMerged = mFusionMapper->mapRead(merged, mapable);
            if(matchMerged){
                matchMerged->addOriginalPair(pair);
                pushMatch(matchMerged);
            } else if(mapable){
                mergedRC = merged->reverseComplement();
                Match* matchMergedRC = mFusionMapper->mapRead(mergedRC, mapable);
                if(matchMergedRC){
                    matchMergedRC->addOriginalPair(pair);
                    pushMatch(matchMergedRC);
                }
                delete mergedRC;
            }

            delete pair;
            delete merged;
            continue;
        }
        // else still search R1 and R2 separatedly
        mapable = false;
        Match* matchR1 = mFusionMapper->mapRead(r1, mapable);
        if(matchR1){
            matchR1->addOriginalPair(pair);
            pushMatch(matchR1);
        } else if(mapable){
            rcr1 = r1->reverseComplement();
            Match* matchRcr1 = mFusionMapper->mapRead(rcr1, mapable);
            if(matchRcr1){
                matchRcr1->addOriginalPair(pair);
                matchRcr1->setReversed(true);
                pushMatch(matchRcr1);
            }
            delete rcr1;
        }
        mapable = false;
        Match* matchR2 = mFusionMapper->mapRead(r2, mapable);
        if(matchR2){
            matchR2->addOriginalPair(pair);
            pushMatch(matchR2);
        } else if(mapable) {
            rcr2 = r2->reverseComplement();
            Match* matchRcr2 = mFusionMapper->mapRead(rcr2, mapable);
            if(matchRcr2){
                matchRcr2->addOriginalPair(pair);
                matchRcr2->setReversed(true);
                pushMatch(matchRcr2);
            }
            delete rcr2;
        }
        delete pair;
    }

    delete pack->data;
    delete pack;

    return true;
}

void PairEndScanner::initPackRepository() {
    mRepo.packBuffer = new ReadPairPack*[PACK_NUM_LIMIT];
    memset(mRepo.packBuffer, 0, sizeof(ReadPairPack*)*PACK_NUM_LIMIT);
    mRepo.writePos = 0;
    mRepo.readPos = 0;
    mRepo.readCounter = 0;
    
}

void PairEndScanner::destroyPackRepository() {
    delete mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void PairEndScanner::producePack(ReadPairPack* pack){
    std::unique_lock<std::mutex> lock(mRepo.mtx);
    while(((mRepo.writePos + 1) % PACK_NUM_LIMIT)
        == mRepo.readPos) {
        mRepo.repoNotFull.wait(lock);
    }

    mRepo.packBuffer[mRepo.writePos] = pack;
    mRepo.writePos++;

    if (mRepo.writePos == PACK_NUM_LIMIT)
        mRepo.writePos = 0;

    mRepo.repoNotEmpty.notify_all();
    lock.unlock();
}

void PairEndScanner::consumePack(){
    ReadPairPack* data;
    std::unique_lock<std::mutex> lock(mRepo.mtx);
    // read buffer is empty, just wait here.
    while(mRepo.writePos == mRepo.readPos) {
        if(mProduceFinished){
            lock.unlock();
            return;
        }
        mRepo.repoNotEmpty.wait(lock);
    }

    data = mRepo.packBuffer[mRepo.readPos];
    (mRepo.readPos)++;
    lock.unlock();

    scanPairEnd(data);


    if (mRepo.readPos >= PACK_NUM_LIMIT)
        mRepo.readPos = 0;

    mRepo.repoNotFull.notify_all();
}

void PairEndScanner::producerTask()
{
    int slept = 0;
    ReadPair** data = new ReadPair*[PACK_SIZE];
    memset(data, 0, sizeof(ReadPair*)*PACK_SIZE);
    FastqReaderPair reader(mRead1File, mRead2File);
    int count=0;
    while(true){
        ReadPair* read = reader.read();
        if(!read){
            // the last pack
            if(count>0){
                ReadPairPack* pack = new ReadPairPack;
                pack->data = data;
                pack->count = count;
                producePack(pack);
            }
            data = NULL;
            break;
        }
        data[count] = read;
        count++;
        // a full pack
        if(count == PACK_SIZE){
            ReadPairPack* pack = new ReadPairPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            //re-initialize data for next pack
            data = new ReadPair*[PACK_SIZE];
            memset(data, 0, sizeof(ReadPair*)*PACK_SIZE);
            // reset count to 0
            count = 0;
            // if the consumer is far behind this producer, sleep and wait to limit memory usage
            while(mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT){
                //cout<<"sleep"<<endl;
                slept++;
                usleep(1000);
            }
        }
    }

    std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
    mProduceFinished = true;
    lock.unlock();

    // if the last data initialized is not used, free it
    if(data != NULL)
        delete data;
}

void PairEndScanner::consumerTask()
{
    while(true) {
        std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
        if(mProduceFinished && mRepo.writePos == mRepo.readPos){
            lock.unlock();
            break;
        }
        if(mProduceFinished){
            consumePack();
            lock.unlock();
        } else {
            lock.unlock();
            consumePack();
        }
    }
}

void PairEndScanner::textReport(vector<Fusion>& fusionList, vector<Match*> *fusionMatches) {
}

void PairEndScanner::htmlReport(vector<Fusion>& fusionList, vector<Match*> *fusionMatches) {
    if(mHtmlFile == "")
        return;

    HtmlReporter reporter(mHtmlFile, fusionList, fusionMatches);
    reporter.run();
}

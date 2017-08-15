#include "sescanner.h"
#include "fastqreader.h"
#include <iostream>
#include "htmlreporter.h"
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include "util.h"

SingleEndScanner::SingleEndScanner(string fusionFile, string refFile, string read1File, string html, int threadNum){
    mRead1File = read1File;
    mFusionFile = fusionFile;
    mRefFile = refFile;
    mHtmlFile = html;
    mProduceFinished = false;
    mThreadNum = threadNum;
    mFusionMapper = NULL;
}

SingleEndScanner::~SingleEndScanner() {
    if(mFusionMapper != NULL) {
        delete mFusionMapper;
        mFusionMapper = NULL;
    }
}

bool SingleEndScanner::scan(){

    fusionList = Fusion::parseCsv(mFusionFile);
    mFusionMapper = new FusionMapper(mRefFile, fusionList);

    fusionMatches = new vector<Match*>[fusionList.size()];
    for(int i=0;i<fusionList.size();i++){
        fusionMatches[i] = vector<Match*>();
    }

    initPackRepository();
    std::thread producer(std::bind(&SingleEndScanner::producerTask, this));

    std::thread** threads = new thread*[mThreadNum];
    for(int t=0; t<mThreadNum; t++){
        threads[t] = new std::thread(std::bind(&SingleEndScanner::consumerTask, this));
    }

    producer.join();
    for(int t=0; t<mThreadNum; t++){
        threads[t]->join();
    }

    for(int t=0; t<mThreadNum; t++){
        delete threads[t];
        threads[t] = NULL;
    }

    mFusionMapper->removeAlignables(fusionMatches, fusionList.size());

    // sort the matches to make the pileup more clear
    for(int i=0;i<fusionList.size();i++){
        sort(fusionMatches[i].begin(), fusionMatches[i].end(), Match::greater); 
    }

    textReport(fusionList, fusionMatches);
    htmlReport(fusionList, fusionMatches);

    // free memory
    for(int i=0;i<fusionList.size();i++){
        fusionMatches[i].clear();
    }
    return true;
}

void SingleEndScanner::pushMatch(int i, Match* m){
    std::unique_lock<std::mutex> lock(mFusionMtx);
    fusionMatches[i].push_back(m);
    lock.unlock();
}

bool SingleEndScanner::scanSingleEnd(ReadPack* pack){
    for(int p=0;p<pack->count;p++){
        Read* r1 = pack->data[p];
        Read* rcr1 = r1->reverseComplement();
        Match* matchR1 = mFusionMapper->mapRead(r1);
        if(matchR1){
            matchR1->addOriginalRead(r1);
            pushMatch(0, matchR1);
        }
        Match* matchRcr1 = mFusionMapper->mapRead(rcr1);
        if(matchRcr1){
            matchRcr1->addOriginalRead(r1);
            matchRcr1->setReversed(true);
            pushMatch(0, matchRcr1);
        }
        delete r1;
        delete rcr1;
    }

    delete pack->data;
    delete pack;

    return true;
}

void SingleEndScanner::initPackRepository() {
    mRepo.packBuffer = new ReadPack*[PACK_NUM_LIMIT];
    memset(mRepo.packBuffer, 0, sizeof(ReadPack*)*PACK_NUM_LIMIT);
    mRepo.writePos = 0;
    mRepo.readPos = 0;
    mRepo.readCounter = 0;
    
}

void SingleEndScanner::destroyPackRepository() {
    delete mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void SingleEndScanner::producePack(ReadPack* pack){
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

void SingleEndScanner::consumePack(){
    ReadPack* data;
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

    scanSingleEnd(data);


    if (mRepo.readPos >= PACK_NUM_LIMIT)
        mRepo.readPos = 0;

    mRepo.repoNotFull.notify_all();
}

void SingleEndScanner::producerTask()
{
    int slept = 0;
    Read** data = new Read*[PACK_SIZE];
    memset(data, 0, sizeof(Read*)*PACK_SIZE);
    FastqReader reader1(mRead1File);
    int count=0;
    while(true){
        Read* read = reader1.read();
        if(!read){
            // the last pack
            if(count>0){
                ReadPack* pack = new ReadPack;
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
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            //re-initialize data for next pack
            data = new Read*[PACK_SIZE];
            memset(data, 0, sizeof(Read*)*PACK_SIZE);
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

void SingleEndScanner::consumerTask()
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

void SingleEndScanner::textReport(vector<Fusion>& fusionList, vector<Match*> *fusionMatches) {
    /*
    //output result
    for(int i=0;i<fusionList.size();i++){
        vector<Match*> matches = fusionMatches[i];
        if(matches.size()>0){
            cout<<endl<<"---------------"<<endl;
            fusionList[i].print();
            for(int m=0; m<matches.size(); m++){
                cout<<m+1<<", ";
                matches[m]->print(fusionList[i].mLeft.length(), fusionList[i].mCenter.length(), fusionList[i].mRight.length());
            }
        }
    }
    */
}

void SingleEndScanner::htmlReport(vector<Fusion>& fusionList, vector<Match*> *fusionMatches) {
    if(mHtmlFile == "")
        return;

    HtmlReporter reporter(mHtmlFile, fusionList, fusionMatches);
    reporter.run();
}

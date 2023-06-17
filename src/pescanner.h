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
#include "fusionmapper.h"


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
    PairEndScanner(string fusionFile, string refFile, string read1File, string read2File, string html, string json, int threadnum);
    ~PairEndScanner();
    bool scan();
    void textReport();
    void htmlReport(FusionMapper* mFusionMapper);
    void jsonReport(FusionMapper* mFusionMapper);

private:
    bool scanPairEnd(ReadPairPack* pack);
    void initPackRepository();
    void destroyPackRepository();
    void producePack(ReadPairPack* pack);
    void consumePack();
    void producerTask();
    void consumerTask();
    void pushMatch(Match* m, FusionMapper* mFusionMapper);

private:
    string mFusionFile;
    string mRefFile;
    string mRead1File;
    string mRead2File;
    string mHtmlFile;
    string mJsonFile;
    ReadPairRepository mRepo;
    bool mProduceFinished;
    std::mutex mFusionMtx;
    int mThreadNum;
    //FusionMapper* mFusionMapper;
    vector<FusionMapper *> mFusionMapperV;
};


#endif

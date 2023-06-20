#include "pescanner.h"
#include "fastqreader.h"
#include <iostream>
#include "htmlreporter.h"
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include "util.h"
#include "jsonreporter.h"
#define PACK_SIZE 1024

PairEndScanner::PairEndScanner(string fusionFile, string refFile, string read1File, string read2File, string html, string json, int threadNum){
    mRead1File = read1File;
    mRead2File = read2File;
    mFusionFile = fusionFile;
    mRefFile = refFile;
    mHtmlFile = html;
    mJsonFile = json;
    mThreadNum = threadNum;
    process_number = 0;
}

PairEndScanner::~PairEndScanner() {
}

void PairEndScanner::scanPairEndWrapper(vector<ReadPair* >pack, FusionMapper* mFusionMapper){
    if(process_number >= mThreadNum){
        unique_lock<mutex> lock(mtx);
        while (process_number >= mThreadNum)
            condition.wait(lock);
    }
    process_number += 1;

    scanPairEnd(pack, mFusionMapper);

    process_number -= 1;
    condition.notify_one();
}

bool PairEndScanner::scan(){

    FastqReaderPair reader(mRead1File, mRead2File);
    int count = 0;
    while(true){
        ReadPair* read = reader.read();
        if(!read)
            break;
        ReadPairV.push_back(read);
        count++;
    }
    cerr << "!!! load " << count << " readpairs" << endl;

    int csv = 0;
    string suf;
    string bakHtml = mHtmlFile;
    int pos = bakHtml.find(".html");
    if(pos >= 0)
        bakHtml.replace(pos, 5, "");
    string bakJson = mJsonFile;

    ifstream file;
    file.open(mFusionFile.c_str(), ifstream::in);
    const int maxLine = 4096;
    char line[maxLine];
    FusionMapper* mFusionMapper;
    vector<ReadPair *> pack; 
    vector<vector<ReadPair *>> packV;
    vector<thread> th_set;

    while(file.getline(line, maxLine)){
        csv++;
        suf = "_" + to_string(csv);
        cerr << ">>> " << suf << endl; 
        if(bakHtml != "")
            mHtmlFile = bakHtml + suf + ".html";
        if(bakJson != "")
            mJsonFile = bakJson + suf;

        mFusionMapper = new FusionMapper(mRefFile, line);
        int j = 0;
        for(int i = 0; i < count; i++) {
            pack.push_back(ReadPairV[i]);

            if ((i + 1) % PACK_SIZE == 0 || (i + 1) == count){
                packV.push_back(pack);
                pack.clear();

                th_set.push_back(thread(&PairEndScanner::scanPairEndWrapper, this, packV[j], mFusionMapper));
                j++;
            }
            if ((i + 1) % (PACK_SIZE * mThreadNum * 4) == 0 || (i + 1) == count){
                for(auto &th : th_set)
                    th.join();
                th_set.clear();
                process_number = 0;

                packV.clear();
                j = 0;
            }
            if ((i + 1) % 1000000 == 0)
                cerr << "^^^ parse 1M pairs" <<endl;
        }

        mFusionMapper->filterMatches();
        mFusionMapper->sortMatches();
        mFusionMapper->clusterMatches();

        htmlReport(mFusionMapper);
        jsonReport(mFusionMapper);

        mFusionMapper->freeMatches();
    }
    return true;
}

void PairEndScanner::pushMatch(Match* m, FusionMapper* mFusionMapper){
    std::unique_lock<std::mutex> lock(mFusionMtx);
    mFusionMapper->addMatch(m);
    lock.unlock();
}

bool PairEndScanner::scanPairEnd(vector<ReadPair *> pack, FusionMapper* mFusionMapper){
    for(int p = 0; p < pack.size(); p++){
        ReadPair* pair = pack[p];
        Read* r1 = pair->mLeft;
        Read* r2 = pair->mRight;
        Read* rcr1 = NULL;
        Read* rcr2 = NULL;
        Read* merged = pair->fastMerge();
        Read* mergedRC = NULL;
        bool mapable = false;

        if(merged != NULL) {
            Match* matchMerged = mFusionMapper->mapRead(merged, mapable);
            if(matchMerged){
                matchMerged->addOriginalPair(pair);
                pushMatch(matchMerged, mFusionMapper);
            } else if(mapable){
                mergedRC = merged->reverseComplement();
                Match* matchMergedRC = mFusionMapper->mapRead(mergedRC, mapable);
                if(matchMergedRC){
                    matchMergedRC->addOriginalPair(pair);
                    pushMatch(matchMergedRC, mFusionMapper);
                }
                delete mergedRC;
            }

            delete merged;
            continue;
        }

        mapable = false;
        Match* matchR1 = mFusionMapper->mapRead(r1, mapable);
        if(matchR1){
            matchR1->addOriginalPair(pair);
            pushMatch(matchR1, mFusionMapper);
        } else if(mapable){
            rcr1 = r1->reverseComplement();
            Match* matchRcr1 = mFusionMapper->mapRead(rcr1, mapable);
            if(matchRcr1){
                matchRcr1->addOriginalPair(pair);
                matchRcr1->setReversed(true);
                pushMatch(matchRcr1, mFusionMapper);
            }
            delete rcr1;
        }
        mapable = false;
        Match* matchR2 = mFusionMapper->mapRead(r2, mapable);
        if(matchR2){
            matchR2->addOriginalPair(pair);
            pushMatch(matchR2, mFusionMapper);
        } else if(mapable) {
            rcr2 = r2->reverseComplement();
            Match* matchRcr2 = mFusionMapper->mapRead(rcr2, mapable);
            if(matchRcr2){
                matchRcr2->addOriginalPair(pair);
                matchRcr2->setReversed(true);
                pushMatch(matchRcr2, mFusionMapper);
            }
            delete rcr2;
        }
    }

    return true;
}

void PairEndScanner::htmlReport(FusionMapper* mFusionMapper) {
    if(mHtmlFile == "")
        return;

    HtmlReporter reporter(mHtmlFile, mFusionMapper);
    reporter.run();
}

void PairEndScanner::jsonReport(FusionMapper* mFusionMapper) {
    if(mJsonFile == "")
        return;

    JsonReporter reporter(mJsonFile, mFusionMapper);
    reporter.run();
}

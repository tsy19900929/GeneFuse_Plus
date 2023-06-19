#include <stdio.h>
#include "fastqreader.h"
#include "unittest.h"
#include "fusionscan.h"
#include <time.h>
#include "cmdline.h"
#include <sstream>
#include "util.h"
#include "globalsettings.h"

string command;

int main(int argc, char* argv[]){
    if (argc == 2 && strcmp(argv[1], "test")==0){
        UnitTest tester;
        tester.run();
        return 0;
    }
    cmdline::parser cmd;
    cmd.add<string>("read1", '1', "read1 file name", true, "");
    cmd.add<string>("read2", '2', "read2 file name", true, "");
    cmd.add<string>("fusion", 'f', "list of csv files", true, "");
    cmd.add<string>("ref", 'r', "reference fasta file name", true, "");
    cmd.add<int>("unique", 'u', "specify the least supporting read number is required to report a fusion, default is 2", false, 2);
    cmd.add<string>("html", 'h', "file name to store HTML report; 1th csv file -> _1.html, 2nd csv file -> _2.html ... ", false, "genefuse.html");
    cmd.add<string>("json", 'j', "file name to store JSON report, default is genefuse.json", false, "genefuse.json");
    cmd.add<int>("thread", 't', "worker thread number, default is 8", false, 8);
    cmd.add<int>("deletion", 'd', "specify the least deletion length of a intra-gene deletion to report, default is 50", false, 50);
    cmd.add("output_deletions", 'D', "long deletions are not output by default, enable this option to output them");
    cmd.add("output_untranslated_fusions", 'U', "the fusions that cannot be transcribed or translated are not output by default, enable this option to output them");
    cmd.parse_check(argc, argv);
    string r1file = cmd.get<string>("read1");
    string r2file = cmd.get<string>("read2");
    string fusionFile = cmd.get<string>("fusion");
    string html = cmd.get<string>("html");
    string json = cmd.get<string>("json");
    string refFile = cmd.get<string>("ref");
    int threadNum = cmd.get<int>("thread");
    int unique = cmd.get<int>("unique");
    int deletion = cmd.get<int>("deletion");
    bool outputDeletion = cmd.exist("output_deletions");
    bool outputUntranslated = cmd.exist("output_untranslated_fusions");

    GlobalSettings::setUniqueRequirement(unique);
    GlobalSettings::setDeletionThreshold(deletion);
    GlobalSettings::setOutputDeletions(outputDeletion);
    GlobalSettings::setOutputUntranslated(outputUntranslated);


    if(ends_with(refFile, ".gz") || ends_with(refFile, ".gz")) {
        cout << "reference fasta file should not be compressed.\nplease unzip "<<refFile<<" and try again."<<endl;
        exit(-1);
    }

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    check_file_valid(refFile);
    check_file_valid(r1file);
    if(r2file != "")
        check_file_valid(r2file);
    if(fusionFile != "")
        check_file_valid(fusionFile);

    loginfo("start with " + string( int2str(threadNum )) + " threads");

    time_t t1 = time(NULL);

    FusionScan fs(fusionFile, refFile, r1file, r2file, html, json, threadNum);
    fs.scan();

    time_t t2 = time(NULL);
    printf("\n# %s\n", command.c_str());
    printf("# genefuse v%s, time used: %ld seconds\n", FUSIONSCAN_VER, (t2-t1));

    loginfo("done");
}

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
    cmd.add<string>("read2", '2', "read2 file name", false, "");
    cmd.add<string>("fusion", 'f', "fusion file name, in CSV format", false, "");
    cmd.add<string>("ref", 'r', "reference fasta file name", true, "");
    cmd.add<int>("unique", 'u', "how many unique supporting read is required to report a fusion, default is 2", false, 2);
    cmd.add<string>("html", 'h', "filename of html report, no html report if not specified", false, "");
    cmd.add<int>("thread", 't', "worker thread number, default is 4", false, 4);
    cmd.parse_check(argc, argv);
    string r1file = cmd.get<string>("read1");
    string r2file = cmd.get<string>("read2");
    string fusionFile = cmd.get<string>("fusion");
    string html = cmd.get<string>("html");
    string refFile = cmd.get<string>("ref");
    int threadNum = cmd.get<int>("thread");
    int unique = cmd.get<int>("unique");

    GlobalSettings::setUniqueRequirement(unique);


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

    clock_t t1 = clock();

    FusionScan fs(fusionFile, refFile, r1file, r2file, html, threadNum);
    fs.scan();

    clock_t t2 = clock();
    printf("\n# %s\n", command.c_str());
    printf("# fusionscan v%s, time used: %f ms\n", FUSIONSCAN_VER, (t2-t1)/1000.0);

    loginfo("done");
}
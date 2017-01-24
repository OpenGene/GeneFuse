#include "unittest.h"
#include "editdistance.h"
#include "sequence.h"
#include "fastqreader.h"
#include "fastareader.h"
#include "overlap.h"
#include "read.h"
#include "fusion.h"
#include "indexer.h"
#include <time.h>

UnitTest::UnitTest(){

}

void UnitTest::run(){
    bool passed = true;
    passed &= editdistance_test();
    passed &= Sequence::test();
    passed &= Overlap::test();
    passed &= Fusion::test();
    passed &= Indexer::test();
    printf("\n==========================\n");
    printf("%s\n\n", passed?"PASSED":"FAILED");
}
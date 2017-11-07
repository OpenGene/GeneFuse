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
    passed &= report(editdistance_test(), "editdistance_test");
    passed &= report(Sequence::test(), "Sequence::test");
    passed &= report(Overlap::test(), "Overlap::test");
    passed &= report(Fusion::test(), "Fusion::test");
    passed &= report(Indexer::test(), "Indexer::test");
    passed &= report(FastaReader::test(), "FastaReader::test");
    printf("\n==========================\n");
    printf("%s\n\n", passed?"ALL PASSED":"FAILED");
}

bool UnitTest::report(bool result, string message) {
    printf("%s:%s\n\n", message.c_str(), result?" PASSED":" FAILED");
    return result;
}
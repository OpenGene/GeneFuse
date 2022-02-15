#ifndef GLOBALSETTINGS_H
#define GLOBALSETTINGS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

class GlobalSettings{
public:
    GlobalSettings();

public:
    inline static void setMarkedOnlyForVCF(bool flag){
        markedOnlyForVCF = flag;
    }
    inline static void setUniqueRequirement(int val){
        uniqueRequirement = val;
    }
    inline static void setDeletionThreshold(int val){
        deletionThreshold = val;
    }
    inline static void setOutputDeletions(bool flag){
        outputDeletions = flag;
    }
    inline static void setOutputUntranslated(bool flag){
        outputUntranslated = flag;
    }

public:
    static bool markedOnlyForVCF;
    static int uniqueRequirement;
    static int deletionThreshold;
    static bool outputDeletions;
    static bool outputUntranslated;
    static int skipKeyDupThreshold;
    static int majorGeneKeyRequirement;
    static int minorGeneKeyRequirement;
    static int mismatchThreshold;
};


#endif
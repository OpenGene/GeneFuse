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

public:
    static bool markedOnlyForVCF;
    static int uniqueRequirement;
    static int deletionThreshold;
};


#endif
#include "globalsettings.h"

bool GlobalSettings::markedOnlyForVCF = false;
int GlobalSettings::uniqueRequirement = 2;
int GlobalSettings::deletionThreshold = 50;
bool GlobalSettings::outputDeletions = false;
bool GlobalSettings::outputUntranslated = false;
int GlobalSettings::skipKeyDupThreshold = 5;
int GlobalSettings::majorGeneKeyRequirement = 40;
int GlobalSettings::minorGeneKeyRequirement = 20;
int GlobalSettings::mismatchThreshold = 10;
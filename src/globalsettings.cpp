#include "globalsettings.h"

bool GlobalSettings::markedOnlyForVCF = false;
int GlobalSettings::uniqueRequirement = 2;
int GlobalSettings::deletionThreshold = 50;
bool GlobalSettings::outputDeletions = false;
bool GlobalSettings::outputUntranslated = false;
int GlobalSettings::skipKeyDupThreshold = 10;
int GlobalSettings::majorGeneKeyRequirement = 30;
int GlobalSettings::minorGeneKeyRequirement = 10;
int GlobalSettings::mismatchThreshold = 10;
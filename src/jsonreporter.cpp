#include "jsonreporter.h"
#include "common.h"
#include <chrono>
#include "globalsettings.h"

JsonReporter::JsonReporter(string filename, FusionMapper* mapper){
    mFusionMapper = mapper;
    mFusionResults = mapper->mFusionResults;
    mFilename = filename;
    mFile.open(mFilename.c_str(), ifstream::out);
}

JsonReporter::~JsonReporter(){
    mFile.close();
}

extern string getCurrentSystemTime();
extern string command;

void JsonReporter::run() {
    mFile << "{" << endl;
    mFile << "\t\"command\":\"" << command << "\"," << endl;
    mFile << "\t\"version\":\"" << FUSIONSCAN_VER << "\"," << endl;
    mFile << "\t\"time\":\"" << getCurrentSystemTime() << "\"," << endl;
    mFile << "\t\"fusions\":{";

    bool isFirstMut = true;
    for(int i=0;i<mFusionResults.size();i++){
        FusionResult fusion = mFusionResults[i];
        vector<Match*> matches = fusion.mMatches;
        if(!GlobalSettings::outputDeletions && fusion.isDeletion())
            continue;
        if(fusion.isLeftProteinForward() != fusion.isRightProteinForward()) {
            if(!GlobalSettings::outputUntranslated)
                continue;
        }
        
        if(isFirstMut) {
            mFile << endl;
            isFirstMut = false;
        }
        else
            mFile << "," << endl;

        mFile << "\t\t\"" <<  fusion.mTitle << "\":{" << endl;
            mFile << "\t\t\t\"" <<  "left" << "\":{" << endl;
                mFile << "\t\t\t\t\"" <<  "gene_name" << "\":" << "\"" << fusion.mLeftGene.mName << "\"," << endl;
                mFile << "\t\t\t\t\"" <<  "gene_chr" << "\":" << "\"" << fusion.mLeftGene.mChr << "\"," << endl;
                mFile << "\t\t\t\t\"" <<  "position" << "\":" << fusion.mLeftGP.position << "," << endl;
                mFile << "\t\t\t\t\"" <<  "reference" << "\":" << "\"" << fusion.mLeftRef << "\"," << endl;
                mFile << "\t\t\t\t\"" <<  "ref_ext" << "\":" << "\"" << fusion.mLeftRefExt << "\"," << endl;
                mFile << "\t\t\t\t\"" <<  "pos_str" << "\":" << "\"" << fusion.mLeftPos << "\"," << endl;
                mFile << "\t\t\t\t\"" <<  "exon_or_intron" << "\":" << "\"" << (fusion.mLeftIsExon?"exon":"intron") << "\"," << endl;
                mFile << "\t\t\t\t\"" <<  "exon_or_intron_id" << "\":" << fusion.mLeftExonOrIntronID << "," << endl;
                mFile << "\t\t\t\t\"" <<  "strand" << "\":" << "\"" << (fusion.isLeftProteinForward()?"forward":"reversed") << "\"" << endl;
            mFile << "\t\t\t}, " << endl;
            mFile << "\t\t\t\"" <<  "right" << "\":{" << endl;
                mFile << "\t\t\t\t\"" <<  "gene_name" << "\":" << "\"" << fusion.mRightGene.mName << "\"," << endl;
                mFile << "\t\t\t\t\"" <<  "gene_chr" << "\":" << "\"" << fusion.mRightGene.mChr << "\"," << endl;
                mFile << "\t\t\t\t\"" <<  "position" << "\":" << fusion.mRightGP.position << "," << endl;
                mFile << "\t\t\t\t\"" <<  "reference" << "\":" << "\"" << fusion.mRightRef << "\"," << endl;
                mFile << "\t\t\t\t\"" <<  "ref_ext" << "\":" << "\"" << fusion.mRightRefExt << "\"," << endl;
                mFile << "\t\t\t\t\"" <<  "pos_str" << "\":" << "\"" << fusion.mRightPos << "\"," << endl;
                mFile << "\t\t\t\t\"" <<  "exon_or_intron" << "\":" << "\"" << (fusion.mRightIsExon?"exon":"intron") << "\"," << endl;
                mFile << "\t\t\t\t\"" <<  "exon_or_intron_id" << "\":" << fusion.mRightExonOrIntronID << "," << endl;
                mFile << "\t\t\t\t\"" <<  "strand" << "\":" << "\"" << (fusion.isRightProteinForward()?"forward":"reversed") << "\"" << endl;
            mFile << "\t\t\t}, " << endl;

        mFile << "\t\t\t\"" <<  "unique" << "\":" << fusion.mUnique << "," << endl;
        mFile << "\t\t\t\"" <<  "reads" << "\":[" << endl;
        for(int m=0; m<matches.size(); m++){
            mFile << "\t\t\t\t{" << endl;
            mFile << "\t\t\t\t\t\"break\":" << matches[m]->mReadBreak << "," << endl;
            mFile << "\t\t\t\t\t\"" << "strand" << "\":" << "\"" << (matches[m]->mReversed?"reversed":"forward") << "\"," << endl;
            matches[m]->printReadToJson(mFile, "\t\t\t\t\t");
            mFile << "\t\t\t\t}";
            if(m!=matches.size()-1)
                mFile << ",";
            mFile << endl;
        }
        mFile << "\t\t\t]" << endl;
        mFile << "\t\t}";
    }

    mFile << endl;
    mFile << "\t}" << endl;
    mFile << "}" << endl;
}